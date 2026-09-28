"""Score MISATO static-docking poses against *all* matching crystal ligand copies.

Differences from ``score_static_docking_wsl.py`` (v3):

* Multi-copy references (homo-oligomers) are scored the same way as the
  co-folding benchmark: every crystal copy with the pose's heavy-element
  composition is passed to OST, which assigns the best-matching copy. v3
  excluded these whenever the copy had to be chosen by pose proximity.
* Template transfer. The pose is written into a copy of the crystal entity by
  moving the atoms of one native ligand copy onto the pose coordinates, using
  the pose->native heavy-atom graph mapping. Model and target ligands therefore
  carry identical residue names, atom names, and bonds, so OST's ligand
  ``identity`` check cannot fail on bond-perception differences between the
  QM-derived SDF and the crystal compound. The receptor is the unchanged crystal
  polymer (all non-polymer chains removed), exactly as in the static docking.
* The reference is reduced to heavy atoms before scoring, so native ligands
  deposited with hydrogens (NMR, neutron, riding-H models) match the
  heavy-atom pose.

Composite ligands, peptide ligands, and poses without a heavy-element-matched
crystal residue are carried over as excluded rows; every delivered pose keeps a
row.
"""
from __future__ import annotations

import argparse
import csv
import logging
import signal
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from score_static_docking_wsl import initialize_worker, read_audit, resolve_path

FIELDS = [
    "target_id", "method", "pose_path", "reference_path", "reference_resname",
    "preaudit_status", "preaudit_selection", "status", "bisy_rmsd_angstrom",
    "lddt_pli", "n_target_copies", "assigned_ref_ost_chain", "elapsed_seconds", "error",
]
# Pre-audit outcomes that still have a single-component native ligand in the
# crystal (possibly several copies). Everything else is carried as excluded.
ATTEMPT = {"provisional_pose_rmsd", "ambiguous_crystal_copies"}


def initial_row(audit: dict[str, str]) -> dict[str, object]:
    return {
        "target_id": audit["target_id"], "method": audit["method"],
        "pose_path": audit["pose_path"], "reference_path": audit["reference_path"],
        "reference_resname": "", "preaudit_status": audit["reference_status"],
        "preaudit_selection": audit["reference_selection"],
        "status": "pending" if audit["reference_status"] in ATTEMPT and audit["reference_path"]
        else "excluded_reference_or_graph_check",
        "bisy_rmsd_angstrom": "", "lddt_pli": "", "n_target_copies": "",
        "assigned_ref_ost_chain": "", "elapsed_seconds": "", "error": "",
    }


def _edit_mode():
    import ost.mol as mol

    mode = getattr(mol, "UNBUFFERED_EDIT", None)
    return mode if mode is not None else mol._ost_mol.EditMode.UNBUFFERED_EDIT


def _heavy_signature(residue) -> Counter[str]:
    return Counter(str(atom.element).upper() for atom in residue.atoms
                   if str(atom.element).upper() not in {"H", "D"})


def pose_to_native_names(pose_path: Path, reference: Path):
    """Map pose heavy-atom coordinates onto native atom names.

    Returns (resname, {native_atom_name: (x, y, z)}, element signature).
    """
    from audit_crystal_poses import (heavy_mol, native_graph, reference_residues,
                                     signature, untyped_graph)

    pose = heavy_mol(pose_path)
    wanted = signature([atom.GetAtomicNum() for atom in pose.GetAtoms()])
    copies = [res for res in reference_residues(reference) if res["signature"] == wanted]
    if not copies:
        raise ValueError("no_single_component_native_copy")
    pose_graph = untyped_graph(pose)
    # Connectivity is perceived from each copy's coordinates; a poorly resolved
    # copy can fail, so use the first copy whose graph maps onto the pose.
    template, match = None, None
    for candidate in copies:
        try:
            native = native_graph(candidate)
        except (ValueError, RuntimeError):
            continue
        if native.GetNumAtoms() != pose.GetNumAtoms():
            continue
        match = pose_graph.GetSubstructMatch(untyped_graph(native), useChirality=False)
        if match:
            template = candidate
            break
    if template is None:
        raise ValueError("no_graph_isomorphism_any_copy")
    coords = pose.GetConformer().GetPositions()
    names = {template["atoms"][i].name.strip(): tuple(float(v) for v in coords[match[i]])
             for i in range(len(match))}
    if len(names) != len(match):
        raise ValueError("duplicate_native_atom_names")
    elements = Counter(template["atoms"][i].element.name.upper() for i in range(len(match)))
    return template["label"].split(":")[1].upper(), names, elements


def build_model(ref_ent, carrier, name_to_xyz):
    """Copy the crystal entity, keep polymer chains plus one carrier ligand whose
    atoms are moved onto the pose. Returns (model_entity, model_carrier_residue)."""
    from ost import geom

    model = ref_ent.Copy()
    target_chain, target_num = carrier.chain.name, carrier.number
    editor = model.EditXCS(_edit_mode())
    for chain in list(model.chains):
        kind = str(chain.chain_type)
        if "POLY" in kind and "NON_POLY" not in kind:
            continue
        if chain.name != target_chain:
            editor.DeleteChain(chain)
            continue
        for residue in list(chain.residues):
            if residue.number != target_num:
                editor.DeleteResidue(residue)
    residue = model.FindResidue(target_chain, target_num)
    if not residue.IsValid():
        raise ValueError("carrier_residue_lost")
    present = {atom.name.strip() for atom in residue.atoms}
    missing = set(name_to_xyz) - present
    if missing:
        raise ValueError(f"carrier_missing_atoms_{len(missing)}")
    for atom in list(residue.atoms):
        name = atom.name.strip()
        if name in name_to_xyz:
            editor.SetAtomPos(atom, geom.Vec3(*name_to_xyz[name]))
        else:  # hydrogens / extra altloc atoms absent from the heavy pose
            editor.DeleteAtom(atom)
    del editor
    return model, model.FindResidue(target_chain, target_num)


def _assigned(scorer):
    import math

    if not scorer.assignment:
        return None, None
    trg_i, mdl_i = scorer.assignment[0]
    value = scorer.score_matrix[trg_i, mdl_i]
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None, None
    return (None, None) if math.isnan(value) or math.isinf(value) else (value, trg_i)


def score_target(target_id: str, rows: list[dict[str, object]], repo_root_text: str,
                 timeout_seconds: int) -> list[dict[str, object]]:
    root = Path(repo_root_text)
    start = time.perf_counter()
    from plb_bench.scoring import _get_ligands, _load_mmcif_text
    from ost.mol.alg import ligand_scoring_lddtpli, ligand_scoring_scrmsd

    out = []
    refs = {str(row["reference_path"]) for row in rows}
    if len(refs) != 1:
        return [{**row, "status": "failed_inconsistent_target_reference",
                 "error": "methods have different reference paths"} for row in rows]
    reference = resolve_path(next(iter(refs)), root)
    try:
        signal.alarm(timeout_seconds)
        loaded, _, _ = _load_mmcif_text(reference.read_text(encoding="utf-8"))
        # Heavy atoms only. NMR / neutron / riding-H depositions otherwise give
        # the native ligand hydrogens the heavy-only pose lacks, and OST's
        # full-graph identity check fails.
        from ost import mol as ost_mol
        ref_ent = ost_mol.CreateEntityFromView(loaded.Select("ele!=H and ele!=D"), True)
        ref_ligands = _get_ligands(ref_ent)
        signal.alarm(0)
    except Exception as error:
        signal.alarm(0)
        return [{**row, "status": "failed_receptor_preparation",
                 "elapsed_seconds": round(time.perf_counter() - start, 3),
                 "error": f"{type(error).__name__}: {error}"} for row in rows]

    for row in rows:
        this = dict(row)
        try:
            signal.alarm(timeout_seconds)
            resname, name_to_xyz, elements = pose_to_native_names(
                resolve_path(str(row["pose_path"]), root), reference)
            this["reference_resname"] = resname
            copies = [res for res in ref_ligands
                      if res.name.upper() == resname and _heavy_signature(res) == elements]
            if not copies:
                raise ValueError("no_ost_native_copy")
            this["n_target_copies"] = len(copies)
            model_ent, model_lig = build_model(ref_ent, copies[0].handle
                                               if hasattr(copies[0], "handle") else copies[0],
                                               name_to_xyz)
            results = {}
            for key, scorer_type in (("bisy_rmsd_angstrom", ligand_scoring_scrmsd.SCRMSDScorer),
                                     ("lddt_pli", ligand_scoring_lddtpli.LDDTPLIScorer)):
                scorer = scorer_type(model=model_ent, target=ref_ent, model_ligands=[model_lig],
                                     target_ligands=copies, substructure_match=False)
                results[key] = _assigned(scorer)
            this["bisy_rmsd_angstrom"] = "" if results["bisy_rmsd_angstrom"][0] is None \
                else results["bisy_rmsd_angstrom"][0]
            this["lddt_pli"] = "" if results["lddt_pli"][0] is None else results["lddt_pli"][0]
            trg_i = results["bisy_rmsd_angstrom"][1]
            this["assigned_ref_ost_chain"] = "" if trg_i is None else copies[trg_i].chain.name
            both = this["bisy_rmsd_angstrom"] != "" and this["lddt_pli"] != ""
            this["status"] = "scored_both" if both else "partial_or_no_metrics"
            if not both:
                this["error"] = "OST returned one or no metrics after template transfer"
        except TimeoutError as error:
            this["status"], this["error"] = "failed_runtime_timeout", str(error)
        except ValueError as error:
            this["status"], this["error"] = "failed_mapping", str(error)
        except Exception as error:  # keep a row for every delivered pose
            this["status"], this["error"] = "failed_scoring", f"{type(error).__name__}: {error}"
        finally:
            signal.alarm(0)
        this["elapsed_seconds"] = round(time.perf_counter() - start, 3)
        out.append(this)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--crystal-audit", type=Path, default=Path("misato_output/crystal_pose_all.csv"))
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=int, default=90)
    parser.add_argument("--limit-targets", type=int, default=0, help="Smoke-test first N IDs")
    parser.add_argument("--target-id", action="append", default=[], help="Score one ID; repeatable")
    parser.add_argument("--skip-target-id", action="append", default=[],
                        help="Record a known pathological target as excluded; repeatable")
    parser.add_argument("--resume", action="store_true", help="Append only IDs absent from --out-csv")
    args = parser.parse_args()
    if args.workers < 1 or args.limit_targets < 0 or args.timeout_seconds < 1:
        parser.error("--workers/--timeout-seconds must be positive, --limit-targets nonnegative")
    if args.out_csv.exists() and not args.resume:
        parser.error("Output exists; select a new file or pass --resume")
    if args.resume and not args.out_csv.exists():
        parser.error("--resume requires an existing output CSV")
    try:
        import ost  # noqa: F401
    except ImportError:
        parser.error("OpenStructure is not importable: run inside WSL Ubuntu with conda env plb")

    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root / "scoring" / "plb_bench"))
    audit = read_audit(args.crystal_audit)
    wanted = {item.upper() for item in args.target_id}
    if wanted:
        audit = [row for row in audit if row["target_id"].upper() in wanted]
        if wanted - {row["target_id"].upper() for row in audit}:
            parser.error(f"IDs absent from crystal audit: {sorted(wanted - {r['target_id'].upper() for r in audit})}")
    elif args.limit_targets:
        ids = sorted({row["target_id"] for row in audit
                      if initial_row(row)["status"] == "pending"})[:args.limit_targets]
        audit = [row for row in audit if row["target_id"] in ids]
    seen = set()
    if args.resume:
        seen = {(row["target_id"], row["method"]) for row in read_audit(args.out_csv)}
    rows = [initial_row(row) for row in audit if (row["target_id"], row["method"]) not in seen]
    skip = {item.upper() for item in args.skip_target_id}
    for row in rows:
        if row["target_id"] in skip and row["status"] == "pending":
            row["status"] = "excluded_runtime_guard"
            row["error"] = "Predeclared hard-case exclusion after a prior >20-minute OST stall"
    if len(rows) != len({(row["target_id"], row["method"]) for row in rows}):
        parser.error("Duplicate target/method rows in crystal audit")

    excluded = [row for row in rows if row["status"] != "pending"]
    groups: dict[str, list] = defaultdict(list)
    for row in rows:
        if row["status"] == "pending":
            groups[str(row["target_id"])].append(row)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    with args.out_csv.open("a" if args.resume else "x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if not args.resume:
            writer.writeheader()
        writer.writerows(excluded)
        counts.update(str(row["status"]) for row in excluded)
        handle.flush()
        logging.basicConfig(level=logging.WARNING)
        done = 0
        with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize_worker) as pool:
            futures = {pool.submit(score_target, tid, group, str(root), args.timeout_seconds): tid
                       for tid, group in groups.items()}
            for future in as_completed(futures):
                tid = futures[future]
                try:
                    completed = future.result()
                except Exception as error:
                    completed = [{**row, "status": "failed_worker",
                                  "error": f"{type(error).__name__}: {error}"} for row in groups[tid]]
                writer.writerows(completed)
                counts.update(str(row["status"]) for row in completed)
                handle.flush()
                done += 1
                if done % 250 == 0:
                    print(f"  {done}/{len(groups)} targets; {dict(counts)}", flush=True)
    print(f"Wrote {len(rows)} rows to {args.out_csv}; statuses: {dict(counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
