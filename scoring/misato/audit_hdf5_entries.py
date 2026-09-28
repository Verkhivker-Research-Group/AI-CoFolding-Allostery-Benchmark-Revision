"""Audit requested MISATO HDF5 target groups without bulk trajectory export.

This reports group-level dimensions and the public MISATO ligand-selection
metadata. It deliberately does not claim which frame Lucas supplied to either
docking model; that must come from his scripts/logs or input structures.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import h5py
import numpy as np


FIELDS = [
    "target_id", "hdf5_status", "n_atoms", "n_frames", "coordinate_shape",
    "molecule_boundaries", "last_molecule_start", "last_molecule_atom_count",
    "residue_zero_atom_count", "ligand_selection_rule", "selected_ligand_atom_count",
    "selected_ligand_element_signature", "note",
]
REQUIRED = {
    "trajectory_coordinates", "atoms_number", "atoms_residue",
    "molecules_begin_atom_index",
}


def element_signature(numbers: np.ndarray) -> str:
    return ".".join(str(int(number)) for number in sorted(numbers[numbers > 1]))


def audit_group(target_id: str, group) -> dict[str, object]:
    missing = sorted(REQUIRED - set(group.keys()))
    if missing:
        return {"target_id": target_id, "hdf5_status": "missing_required_fields", "n_atoms": "", "n_frames": "", "coordinate_shape": "", "molecule_boundaries": "", "last_molecule_start": "", "last_molecule_atom_count": "", "residue_zero_atom_count": "", "ligand_selection_rule": "", "selected_ligand_atom_count": "", "selected_ligand_element_signature": "", "note": ", ".join(missing)}
    atoms_number = group["atoms_number"][:]
    atoms_residue = group["atoms_residue"][:]
    boundaries = group["molecules_begin_atom_index"][:]
    coordinates = group["trajectory_coordinates"]
    zero_indices = np.flatnonzero(atoms_residue == 0)
    last_start = int(boundaries[-1]) if len(boundaries) else 0
    if len(zero_indices):
        selected = zero_indices
        rule = "atoms_residue_eq_0"
    else:
        selected = np.arange(last_start, len(atoms_number))
        rule = "fallback_last_molecule_for_peptide_ligand"
    return {
        "target_id": target_id, "hdf5_status": "available", "n_atoms": len(atoms_number),
        "n_frames": coordinates.shape[0], "coordinate_shape": "x".join(str(value) for value in coordinates.shape),
        "molecule_boundaries": ";".join(str(int(value)) for value in boundaries),
        "last_molecule_start": last_start, "last_molecule_atom_count": len(atoms_number) - last_start,
        "residue_zero_atom_count": len(zero_indices), "ligand_selection_rule": rule,
        "selected_ligand_atom_count": len(selected),
        "selected_ligand_element_signature": element_signature(atoms_number[selected]), "note": "",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", type=Path, required=True)
    parser.add_argument("--target-manifest", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    args = parser.parse_args()
    with args.target_manifest.open(newline="", encoding="utf-8") as handle:
        target_ids = [row["target_id"].upper() for row in csv.DictReader(handle) if row["in_valid_primary_cohort"].lower() == "true"]
    rows = []
    with h5py.File(args.h5, "r") as dataset:
        available = {entry.upper(): entry for entry in dataset.keys()}
        for target_id in target_ids:
            if target_id not in available:
                rows.append({"target_id": target_id, "hdf5_status": "target_not_in_hdf5", "n_atoms": "", "n_frames": "", "coordinate_shape": "", "molecule_boundaries": "", "last_molecule_start": "", "last_molecule_atom_count": "", "residue_zero_atom_count": "", "ligand_selection_rule": "", "selected_ligand_atom_count": "", "selected_ligand_element_signature": "", "note": ""})
                continue
            rows.append(audit_group(target_id, dataset[available[target_id]]))
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"rows": len(rows), "statuses": dict(sorted(Counter(row["hdf5_status"] for row in rows).items()))}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
