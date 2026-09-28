"""Score delivered MISATO ligand poses with the existing WSL plb_bench scorer.

This is a separate, conservative path for static-receptor docking. It does not
alter source SDFs, reference CIFs, or the old benchmark. Only a unique native
ligand with a successful graph/coordinate check is primary-score eligible.
"""
from __future__ import annotations

import argparse
import csv
import logging
import re
import signal
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


FIELDS = [
    "target_id", "method", "pose_path", "reference_path", "reference_residue",
    "reference_selection", "preaudit_status", "status", "bisy_rmsd_angstrom",
    "lddt_pli", "selected_ref_ost_chain", "selected_ref_ost_atom_count",
    "elapsed_seconds", "error",
]


def resolve_path(value: str, repo_root: Path) -> Path:
    """Accept the Windows paths stored in the raw manifests from inside WSL."""
    if match := re.match(r"^([A-Za-z]):[\\/](.*)$", value):
        return Path(f"/mnt/{match.group(1).lower()}/" + match.group(2).replace("\\", "/"))
    path = Path(value.replace("\\", "/"))
    return path if path.is_absolute() else repo_root / path


def read_audit(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def raise_timeout(_signum, _frame) -> None:
    raise TimeoutError("OST scoring exceeded per-stage runtime limit")


def initialize_worker() -> None:
    """Avoid a warning per residue in OST 2.11 on large crystal structures."""
    import ost

    ost.PushVerbosityLevel(0)
    logging.getLogger("plb_bench.scoring").setLevel(logging.ERROR)
    signal.signal(signal.SIGALRM, raise_timeout)


def initial_row(audit: dict[str, str]) -> dict[str, object]:
    eligible = (audit["reference_selection"] == "single_candidate"
                and audit["reference_status"] == "provisional_pose_rmsd"
                and bool(audit["reference_residues"]))
    return {
        "target_id": audit["target_id"], "method": audit["method"],
        "pose_path": audit["pose_path"], "reference_path": audit["reference_path"],
        "reference_residue": audit["reference_residues"],
        "reference_selection": audit["reference_selection"],
        "preaudit_status": audit["reference_status"],
        "status": "pending" if eligible else "excluded_reference_or_graph_check",
        "bisy_rmsd_angstrom": "", "lddt_pli": "",
        "selected_ref_ost_chain": "", "selected_ref_ost_atom_count": "",
        "elapsed_seconds": "", "error": "",
    }


def strip_to_polymer(reference: Path, output_pdb: Path) -> int:
    """Remove all bound small molecules/water before adding a predicted pose."""
    import gemmi

    structure = gemmi.read_structure(str(reference))
    if not structure or not structure[0]:
        raise ValueError("empty_reference_structure")
    model = structure[0]
    for chain in model:
        for index in range(len(chain) - 1, -1, -1):
            residue = chain[index]
            if residue.het_flag == "H" or residue.is_water():
                del chain[index]
    n_polymer = sum(len(chain) for chain in model)
    if n_polymer == 0:
        raise ValueError("no_polymer_residues_after_stripping")
    output_pdb.write_text(structure.make_pdb_string(), encoding="utf-8")
    return n_polymer


def score_target(target_id: str, input_rows: list[dict[str, object]],
                 repo_root_text: str, timeout_seconds: int) -> list[dict[str, object]]:
    """One WSL worker per target; reuse prepared receptor for both methods."""
    root = Path(repo_root_text)
    start = time.perf_counter()
    output = []
    # The old package is imported inside the WSL worker. Docker is not used.
    from plb_bench.sanitizer import sanitize
    from plb_bench.scoring import (_load_mmcif_text,
                                   _wire_lig_bonds_into_entity)
    from selected_ligand_metrics import score_selected, write_heavy_only_pose
    from ost import mol as ost_mol

    refs = {str(row["reference_path"]) for row in input_rows}
    selected = {str(row["reference_residue"]) for row in input_rows}
    if len(refs) != 1 or len(selected) != 1:
        return [{**row, "status": "failed_inconsistent_target_reference",
                 "error": "methods have different reference paths or residues"}
                for row in input_rows]
    reference = resolve_path(next(iter(refs)), root)
    try:
        signal.alarm(timeout_seconds)
        if not reference.is_file():
            raise FileNotFoundError(reference)
        with tempfile.TemporaryDirectory(prefix=f"misato_{target_id}_") as temp:
            receptor = Path(temp) / "receptor_polymer.pdb"
            strip_to_polymer(reference, receptor)
            reference_ent, _, _ = _load_mmcif_text(reference.read_text(encoding="utf-8"))
            signal.alarm(0)
            for row in input_rows:
                this = dict(row)
                pose = resolve_path(str(row["pose_path"]), root)
                try:
                    signal.alarm(timeout_seconds)
                    if not pose.is_file():
                        raise FileNotFoundError(pose)
                    heavy_pose = Path(temp) / f"{row['method']}_heavy_only.sdf"
                    write_heavy_only_pose(pose, heavy_pose)
                    model = sanitize(receptor, target_id, ligand_path=heavy_pose)
                    model_ent, _, _ = _load_mmcif_text(model.text)
                    if not hasattr(ost_mol, "STANDARD_EDIT") and hasattr(ost_mol, "UNBUFFERED_EDIT"):
                        # OST 2.11 renamed the edit mode; the legacy helper
                        # otherwise leaves LIG disconnected and full-graph
                        # matching fails for many valid SDFs.
                        ost_mol.STANDARD_EDIT = ost_mol.UNBUFFERED_EDIT
                    _wire_lig_bonds_into_entity(model_ent, model.text)
                    expected_name = str(row["reference_residue"]).split(":")[1]
                    metrics = score_selected(model_ent, reference_ent, pose,
                                             expected_residue_name=expected_name)
                    this["bisy_rmsd_angstrom"] = ("" if metrics["bisy_rmsd"] is None
                                                   else metrics["bisy_rmsd"])
                    this["lddt_pli"] = "" if metrics["lddt_pli"] is None else metrics["lddt_pli"]
                    this["selected_ref_ost_chain"] = metrics["selected_ref_ost_chain"]
                    this["selected_ref_ost_atom_count"] = metrics["selected_ref_ost_atom_count"]
                    this["status"] = ("scored_both" if metrics["bisy_rmsd"] is not None
                                      and metrics["lddt_pli"] is not None else "partial_or_no_metrics")
                    if this["status"] != "scored_both":
                        this["error"] = "OST exact-selected-ligand matching returned one or no metrics"
                except Exception as error:  # preserve a row for every delivered pose
                    this["status"] = ("failed_runtime_timeout" if isinstance(error, TimeoutError)
                                      else "failed_scoring")
                    this["error"] = f"{type(error).__name__}: {error}"
                finally:
                    signal.alarm(0)
                this["elapsed_seconds"] = round(time.perf_counter() - start, 3)
                output.append(this)
    except Exception as error:
        output = [{**row, "status": ("failed_runtime_timeout" if isinstance(error, TimeoutError)
                                      else "failed_receptor_preparation"),
                   "elapsed_seconds": round(time.perf_counter() - start, 3),
                   "error": f"{type(error).__name__}: {error}"}
                  for row in input_rows]
    finally:
        signal.alarm(0)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crystal-audit", type=Path, default=Path("misato_output/crystal_pose_all.csv"))
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=int, default=90,
                        help="Per receptor/pose stage; 90 seconds by default")
    parser.add_argument("--limit-targets", type=int, default=0,
                        help="Smoke-test first N eligible IDs; zero means all")
    parser.add_argument("--target-id", action="append", default=[],
                        help="Score one ID; may be repeated")
    parser.add_argument("--skip-target-id", action="append", default=[],
                        help="Record a known pathological target as excluded; may be repeated")
    parser.add_argument("--resume", action="store_true", help="Append only IDs absent from an existing CSV")
    args = parser.parse_args()
    if args.workers < 1 or args.limit_targets < 0 or args.timeout_seconds < 1:
        parser.error("--workers/--timeout-seconds must be positive and --limit-targets nonnegative")
    if args.out_csv.exists() and not args.resume:
        parser.error("Output exists; select a new file or pass --resume")
    if args.resume and not args.out_csv.exists():
        parser.error("--resume requires an existing output CSV")
    try:
        import ost  # noqa: F401
    except ImportError:
        parser.error("OpenStructure is not importable: run inside WSL Ubuntu with conda env plb")

    root = Path(__file__).resolve().parents[2]
    package_root = root / "scoring" / "plb_bench"
    sys.path.insert(0, str(package_root))
    audit = read_audit(args.crystal_audit)
    target_filter = {item.upper() for item in args.target_id}
    if target_filter:
        audit = [row for row in audit if row["target_id"].upper() in target_filter]
        missing = target_filter - {row["target_id"].upper() for row in audit}
        if missing:
            parser.error(f"target IDs absent from crystal audit: {sorted(missing)}")
    elif args.limit_targets:
        ids = sorted({row["target_id"] for row in audit
                      if initial_row(row)["status"] == "pending"})[:args.limit_targets]
        audit = [row for row in audit if row["target_id"] in ids]
    seen_keys = set()
    if args.resume:
        for row in read_audit(args.out_csv):
            seen_keys.add((row["target_id"], row["method"]))
    rows = [initial_row(row) for row in audit
            if (row["target_id"], row["method"]) not in seen_keys]
    skip_ids = {target_id.upper() for target_id in args.skip_target_id}
    for row in rows:
        if row["target_id"] in skip_ids and row["status"] == "pending":
            row["status"] = "excluded_runtime_guard"
            row["error"] = "Predeclared hard-case exclusion after a prior >20-minute OST stall"
    if len(rows) != len({(row["target_id"], row["method"]) for row in rows}):
        parser.error("Duplicate target/method rows in crystal audit")
    excluded = [row for row in rows if row["status"] != "pending"]
    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        if row["status"] == "pending":
            groups[str(row["target_id"])].append(row)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if args.resume else "x"
    counts: Counter[str] = Counter()
    with args.out_csv.open(mode, newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if not args.resume:
            writer.writeheader()
        writer.writerows(excluded)
        counts.update(str(row["status"]) for row in excluded)
        handle.flush()
        logging.basicConfig(level=logging.WARNING)
        with ProcessPoolExecutor(max_workers=args.workers,
                                 initializer=initialize_worker) as pool:
            futures = {pool.submit(score_target, target_id, group, str(root),
                                   args.timeout_seconds): target_id
                       for target_id, group in groups.items()}
            for future in as_completed(futures):
                target_id = futures[future]
                try:
                    completed = future.result()
                except Exception as error:
                    completed = [{**row, "status": "failed_worker",
                                  "error": f"{type(error).__name__}: {error}"}
                                 for row in groups[target_id]]
                writer.writerows(completed)
                counts.update(str(row["status"]) for row in completed)
                handle.flush()
    print(f"Wrote {len(rows)} rows to {args.out_csv}; statuses: {dict(counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
