"""Rescore delivered MISATO SDF poses with gnina --score_only; never redock.

The crystal receptor retains polymer chains, except graph-verified short
polymer native ligands, and removes non-polymer/branched species and water.
The output is one resumable status row per delivered pose. CNNscore is *not*
the method-native DiffDock confidence logit.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import gemmi

try:
    from .audit_crystal_poses import heavy_mol, signature
    from .audit_multiresidue_recovery import multiresidue_candidates
    from .score_multiresidue_recovery_wsl import graph_valid_copies
    from .score_static_docking_wsl import resolve_path
except ImportError:
    from audit_crystal_poses import heavy_mol, signature
    from audit_multiresidue_recovery import multiresidue_candidates
    from score_multiresidue_recovery_wsl import graph_valid_copies
    from score_static_docking_wsl import resolve_path


FIELDS = ["target_id", "method", "status", "cnnscore", "gnina_version",
          "receptor_mode", "removed_ligand_subchains", "error"]
CNN_SCORE = re.compile(r"^CNNscore:\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*$",
                       re.MULTILINE)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def parse_cnnscores(output: str, expected: int) -> list[float]:
    values = CNN_SCORE.findall(output)
    if len(values) != expected:
        raise ValueError(f"expected_{expected}_CNNscores_got_{len(values)}")
    scores = [float(value) for value in values]
    if any(not 0 <= value <= 1 for value in scores):
        raise ValueError("CNNscore_out_of_range")
    return scores


def gnina_environment(binary: Path, cuda_libs: Path | None) -> dict[str, str]:
    env = os.environ.copy()
    if cuda_libs:
        locations = sorted((cuda_libs / "nvidia").glob("*/lib"))
        if not locations:
            raise FileNotFoundError(f"No NVIDIA libraries under {cuda_libs}")
        env["LD_LIBRARY_PATH"] = ":".join(map(str, locations)) + \
            (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    return env


def verified_ligand_subchains(reference: Path,
                              graph_valid_pose_paths: list[Path]) -> set[str]:
    """Remove every graph-identical short-polymer native copy, not the receptor."""
    removed: set[str] = set()
    signatures_seen = set()
    for path in graph_valid_pose_paths:
        pose = heavy_mol(path)
        wanted = signature([atom.GetAtomicNum() for atom in pose.GetAtoms()])
        fingerprint = tuple(sorted(wanted.items()))
        if fingerprint in signatures_seen:
            continue
        signatures_seen.add(fingerprint)
        matches = graph_valid_copies(pose, multiresidue_candidates(reference, wanted))
        for candidate in matches:
            removed.update(candidate["subchains"])
    return removed


def prepare_receptor(reference: Path, destination: Path,
                     ligand_subchains: set[str]) -> None:
    structure = gemmi.read_structure(str(reference))
    if not structure or not structure[0]:
        raise ValueError("empty_crystal_structure")
    model = structure[0]
    for chain in model:
        for index in range(len(chain) - 1, -1, -1):
            residue = chain[index]
            if residue.entity_type != gemmi.EntityType.Polymer or \
                    residue.subchain in ligand_subchains:
                del chain[index]
    structure.remove_hydrogens()
    if not any(len(chain) for chain in model):
        raise ValueError("no_receptor_polymer_after_stripping")
    destination.write_text(structure.make_pdb_string(), encoding="utf-8")


def score_target(target_id: str, rows: list[dict[str, str]], all_target_rows: list[dict[str, str]], root: Path,
                 crystal_dir: Path, binary: Path, env: dict[str, str],
                 multi_eligible: set[tuple[str, str]], timeout_seconds: int,
                 no_gpu: bool, version: str, batch: bool) -> list[dict[str, str]]:
    result = []
    try:
        reference = crystal_dir / f"{target_id[:4]}.cif"
        graph_paths = [resolve_path(row["pose_path"], root) for row in all_target_rows
                       if (target_id, row["method"]) in multi_eligible]
        removed = verified_ligand_subchains(reference, graph_paths) if graph_paths else set()
        with tempfile.TemporaryDirectory(prefix=f"misato_gnina_{target_id}_") as directory:
            receptor = Path(directory) / "receptor.pdb"
            prepare_receptor(reference, receptor, removed)
            if batch and len(rows) > 1:
                pieces = [resolve_path(row["pose_path"], root).read_bytes() for row in rows]
                if all(piece.count(b"$$$$") == 1 for piece in pieces):
                    combined = Path(directory) / "delivered_poses.sdf"
                    combined.write_bytes(b"\n".join(pieces))
                    command = [str(binary), "--receptor", str(receptor),
                               "--ligand", str(combined), "--score_only", "--cpu", "1"]
                    if no_gpu:
                        command.append("--no_gpu")
                    completed = subprocess.run(command, capture_output=True, text=True,
                                               timeout=timeout_seconds * len(rows), env=env,
                                               check=False)
                    output = completed.stdout + "\n" + completed.stderr
                    if completed.returncode:
                        raise RuntimeError(f"gnina_batch_exit_{completed.returncode}: " +
                                           output[-500:].replace("\n", " "))
                    scores = parse_cnnscores(output, len(rows))
                    return [{"target_id": target_id, "method": row["method"],
                             "status": "scored", "cnnscore": score,
                             "gnina_version": version,
                             "receptor_mode": "static_crystal_polymer",
                             "removed_ligand_subchains": ";".join(sorted(removed)),
                             "error": ""}
                            for row, score in zip(rows, scores)]
            for row in rows:
                item = {field: "" for field in FIELDS}
                item.update(target_id=target_id, method=row["method"], status="pending",
                            gnina_version=version, receptor_mode="static_crystal_polymer",
                            removed_ligand_subchains=";".join(sorted(removed)))
                try:
                    pose = resolve_path(row["pose_path"], root)
                    command = [str(binary), "--receptor", str(receptor),
                               "--ligand", str(pose), "--score_only", "--cpu", "1"]
                    if no_gpu:
                        command.append("--no_gpu")
                    completed = subprocess.run(command, capture_output=True, text=True,
                                               timeout=timeout_seconds, env=env, check=False)
                    output = completed.stdout + "\n" + completed.stderr
                    if completed.returncode:
                        raise RuntimeError(f"gnina_exit_{completed.returncode}: " +
                                           output[-500:].replace("\n", " "))
                    item["cnnscore"] = parse_cnnscores(output, 1)[0]
                    item["status"] = "scored"
                except Exception as error:  # noqa: BLE001 - retain per-pose failure provenance
                    item["status"] = "failed"
                    item["error"] = f"{type(error).__name__}: {error}"[:800]
                result.append(item)
    except Exception as error:  # noqa: BLE001 - retain per-target failure provenance
        result = [{"target_id": target_id, "method": row["method"],
                   "status": "receptor_failed", "error": f"{type(error).__name__}: {error}"[:800]}
                  for row in rows]
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allcopies", type=Path,
                        default=Path("misato_output/misato_score_allcopies_v1.csv"))
    parser.add_argument("--recovery-audit", type=Path,
                        default=Path("misato_output/multiresidue_recovery_all_v2.csv"))
    parser.add_argument("--crystal-dir", type=Path,
                        default=Path("misato_output/rcsb_asymmetric_unit"))
    parser.add_argument("--gnina", type=Path,
                        default=Path("misato_output/tools/gnina.1.3.3.cuda12.8.static"))
    parser.add_argument("--cuda-libs", type=Path,
                        default=Path("misato_output/tools/gnina_cuda_libs"))
    parser.add_argument("--no-gpu", action="store_true")
    parser.add_argument("--batch", action="store_true",
                        help="Score the two delivered SDFs per target in one gnina invocation")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout-seconds", type=int, default=180)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--target-id", action="append", default=[])
    parser.add_argument("--limit-targets", type=int, default=0)
    args = parser.parse_args()
    if args.workers < 1 or args.timeout_seconds < 1:
        parser.error("workers and timeout must be positive")
    if args.out_csv.exists() and not args.resume:
        parser.error("Output exists; choose another file or --resume")
    if args.resume and not args.out_csv.exists():
        parser.error("--resume requires an existing output")
    root = Path(__file__).resolve().parents[2]
    binary = args.gnina.resolve()
    if not binary.is_file():
        parser.error(f"gnina binary missing: {binary}")
    env = gnina_environment(binary, None if args.no_gpu else args.cuda_libs.resolve())
    version_check = subprocess.run([str(binary), "--version"], capture_output=True,
                                   text=True, env=env, check=False)
    if version_check.returncode:
        parser.error("gnina is not executable: " + version_check.stderr[-500:])
    version = version_check.stdout.strip().splitlines()[0]
    multi_eligible = {(row["target_id"], row["method"]) for row in
                      read_rows(args.recovery_audit) if row["graph_status"] == "exact_full_graph"}
    wanted = {target.upper() for target in args.target_id}
    by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_rows(args.allcopies):
        if not wanted or row["target_id"] in wanted:
            by_id[row["target_id"]].append(row)
    if args.limit_targets:
        by_id = {key: by_id[key] for key in sorted(by_id)[:args.limit_targets]}
    all_target_rows = {key: list(rows) for key, rows in by_id.items()}
    prior = {(row["target_id"], row["method"]) for row in read_rows(args.out_csv)} \
        if args.resume else set()
    for target_id in list(by_id):
        by_id[target_id] = [row for row in by_id[target_id]
                            if (target_id, row["method"]) not in prior]
        if not by_id[target_id]:
            del by_id[target_id]
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    with args.out_csv.open("a" if args.resume else "x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if not args.resume:
            writer.writeheader()
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(score_target, target_id, rows, all_target_rows[target_id], root,
                                   args.crystal_dir.resolve(), binary, env, multi_eligible,
                                   args.timeout_seconds, args.no_gpu, version, args.batch): target_id
                       for target_id, rows in by_id.items()}
            for future in as_completed(futures):
                target_id = futures[future]
                try:
                    completed = future.result()
                except Exception as error:  # noqa: BLE001 - retain worker failure provenance
                    completed = [{"target_id": target_id, "method": row["method"],
                                  "status": "worker_failed",
                                  "error": f"{type(error).__name__}: {error}"[:800]}
                                 for row in by_id[target_id]]
                writer.writerows(completed)
                counts.update(str(row["status"]) for row in completed)
                handle.flush()
    print(f"Wrote {sum(counts.values())} gnina rows: {dict(counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
