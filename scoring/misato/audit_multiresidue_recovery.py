"""Read-only eligibility audit for MISATO peptide/branched/composite references.

Only exact heavy-element and full untyped-graph matches are called recoverable.
This script does not merge residues or score with OST.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import gemmi
import numpy as np

try:
    from .audit_crystal_poses import (
        direct_pose_rmsd,
        heavy_mol,
        native_graph,
        reference_residues,
        signature,
    )
except ImportError:  # direct script execution
    from audit_crystal_poses import (
        direct_pose_rmsd,
        heavy_mol,
        native_graph,
        reference_residues,
        signature,
    )


FIELDS = ["target_id", "method", "old_status", "candidate_kind", "candidate_label",
          "candidate_residue_count", "exact_element_candidates", "graph_status",
          "direct_rmsd_angstrom", "candidate_heavy_atoms", "pose_heavy_atoms"]


def selected_atoms(residues) -> list:
    """One heavy conformer per residue/atom name, retaining residue order."""
    atoms = []
    for residue in residues:
        by_name = {}
        for atom in residue:
            if atom.element.atomic_number <= 1:
                continue
            name = atom.name.strip()
            altloc = str(atom.altloc).strip("\x00 ")
            if name not in by_name or altloc in {"", "A"}:
                by_name[name] = atom
        atoms.extend(by_name.values())
    return atoms


def multiresidue_candidates(path: Path, wanted: Counter[int]) -> list[dict]:
    structure = gemmi.read_structure(str(path))
    if not structure or not structure[0]:
        return []
    candidates = []
    for chain in structure[0]:
        by_subchain = defaultdict(list)
        for residue in chain:
            if residue.entity_type in (gemmi.EntityType.Polymer, gemmi.EntityType.Branched):
                by_subchain[residue.subchain].append(residue)
        for subchain, residues in by_subchain.items():
            if not (2 <= len(residues) <= 30):
                continue
            atoms = selected_atoms(residues)
            if len(atoms) < 2:
                continue
            kind = ("short_polymer" if residues[0].entity_type == gemmi.EntityType.Polymer
                    else "branched_glycan")
            candidate_signature = signature([a.element.atomic_number for a in atoms])
            if candidate_signature != wanted:
                continue
            candidates.append({"kind": kind, "label": f"{chain.name}:{subchain}",
                               "residue_count": len(residues), "atoms": atoms,
                               "subchains": [subchain],
                               "signature": candidate_signature})
    # The old audit already found pairs of non-polymer residues by element
    # signature. Rebuild the same pairs here to test full graph connectivity.
    singles = reference_residues(path)
    for a, b in combinations(singles, 2):
        if a["signature"] + b["signature"] != wanted:
            continue
        # Elemental equality alone can pair different receptor chains whose
        # components are many angstroms apart. Require a plausible covalent
        # contact before testing the assembled molecular graph.
        minimum_distance = float(np.min(np.linalg.norm(
            a["coords"][:, None, :] - b["coords"][None, :, :], axis=2)))
        if minimum_distance > 2.3:
            continue
        atoms = a["atoms"] + b["atoms"]
        candidates.append({"kind": "nonpolymer_pair", "label": a["label"] + "+" + b["label"],
                           "residue_count": 2, "atoms": atoms, "subchains": [],
                           "signature": a["signature"] + b["signature"]})
    return candidates


def audit_pose(row: dict[str, str], candidates: list[dict]) -> dict[str, object]:
    pose = heavy_mol(Path(row["pose_path"]))
    wanted = signature([a.GetAtomicNum() for a in pose.GetAtoms()])
    exact = [candidate for candidate in candidates if candidate["signature"] == wanted]
    result: dict[str, object] = {
        "target_id": row["target_id"], "method": row["method"], "old_status": row["status"],
        "candidate_kind": "", "candidate_label": "", "candidate_residue_count": "",
        "exact_element_candidates": len(exact), "graph_status": "no_exact_element_candidate",
        "direct_rmsd_angstrom": "", "candidate_heavy_atoms": "",
        "pose_heavy_atoms": pose.GetNumAtoms(),
    }
    outcomes = []
    for candidate in exact:
        try:
            native = native_graph(candidate)
            if native.GetNumBonds() != pose.GetNumBonds():
                outcomes.append((candidate, None, "bond_count_mismatch"))
                continue
            rmsd, _, issue = direct_pose_rmsd(pose, native)
            outcomes.append((candidate, rmsd, issue or "exact_full_graph"))
        except (ValueError, RuntimeError) as error:
            outcomes.append((candidate, None, type(error).__name__))
    matched = [item for item in outcomes if item[1] is not None]
    if matched:
        candidate, rmsd, issue = min(matched, key=lambda item: item[1])
    elif outcomes:
        candidate, rmsd, issue = outcomes[0]
    else:
        return result
    result.update({"candidate_kind": candidate["kind"], "candidate_label": candidate["label"],
                   "candidate_residue_count": candidate["residue_count"],
                   "candidate_heavy_atoms": len(candidate["atoms"]),
                   "graph_status": issue,
                   "direct_rmsd_angstrom": "" if rmsd is None else round(rmsd, 4)})
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=Path("misato_output/misato_score_allcopies_v1.csv"))
    parser.add_argument("--crystal-dir", type=Path, default=Path("misato_output/rcsb_asymmetric_unit"))
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--limit-targets", type=int, default=0)
    args = parser.parse_args()
    if args.out_csv.exists() or args.limit_targets < 0:
        parser.error("Choose a new output file and nonnegative limit")
    with args.scores.open(newline="", encoding="utf-8-sig") as handle:
        score_rows = [row for row in csv.DictReader(handle)
                      if row["preaudit_selection"] in
                      {"no_native_element_match", "composite_reference_review"}]
    by_id = defaultdict(list)
    for row in score_rows:
        by_id[row["target_id"]].append(row)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    with args.out_csv.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for index, (target_id, rows) in enumerate(sorted(by_id.items())):
            if args.limit_targets and index >= args.limit_targets:
                break
            path = args.crystal_dir / f"{target_id[:4]}.cif"
            try:
                first_pose = heavy_mol(Path(rows[0]["pose_path"]))
                wanted = signature([a.GetAtomicNum() for a in first_pose.GetAtoms()])
                candidates = multiresidue_candidates(path, wanted)
                for row in rows:
                    outcome = audit_pose(row, candidates)
                    writer.writerow(outcome)
                    counts[str(outcome["graph_status"])] += 1
            except Exception as error:  # noqa: BLE001 - retain per-pose audit failures
                for row in rows:
                    writer.writerow({"target_id": target_id, "method": row["method"],
                                     "old_status": row["status"], "graph_status":
                                     f"audit_error:{type(error).__name__}:{error}"})
                    counts["audit_error"] += 1
            handle.flush()
    print(f"Audited {min(len(by_id), args.limit_targets or len(by_id))} targets: {dict(counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
