"""Validate ligand graph agreement between primary MISATO docking outputs.

This is a pre-scoring safeguard. It compares EquiBind and DiffDock output
ligands for every target in the primary cohort; it does not compare either
prediction with the experimental reference ligand.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


FIELDS = [
    "target_id", "equibind_path", "diffdock_path", "equibind_status",
    "diffdock_status", "equibind_heavy_atoms", "diffdock_heavy_atoms",
    "equibind_canonical_smiles", "diffdock_canonical_smiles",
    "equibind_connectivity_smiles", "diffdock_connectivity_smiles",
    "equibind_element_signature", "diffdock_element_signature", "graph_match",
    "connectivity_match", "element_signature_match", "disposition", "note",
]


def descriptor(path_text: str) -> tuple[str, int | str, str, str, str, str]:
    if not path_text:
        return "missing_path", "", "", "", "", ""
    path = Path(path_text)
    try:
        mol = Chem.MolFromMolFile(str(path), sanitize=True, removeHs=False)
    except Exception as error:
        return "rdkit_error", "", "", "", "", type(error).__name__
    if mol is None:
        # Some valid organophosphates cannot pass RDKit's standard valence
        # sanitizer. Keep them auditable rather than categorically failing.
        mol = Chem.MolFromMolFile(str(path), sanitize=False, removeHs=False)
        if mol is None:
            return "unreadable_sdf", "", "", "", "", "RDKit returned None"
        status = "unsanitized_only"
    else:
        status = "ok"
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    elements = ".".join(str(atom.GetAtomicNum()) for atom in sorted(heavy_atoms, key=lambda atom: atom.GetAtomicNum()))
    try:
        no_hydrogens = Chem.RemoveHs(mol)
        try:
            canonical = Chem.MolToSmiles(no_hydrogens, canonical=True, isomericSmiles=True)
            connectivity = Chem.MolToSmiles(no_hydrogens, canonical=True, isomericSmiles=False)
        except Exception:
            canonical, connectivity = "", ""
        return (
            status,
            len(heavy_atoms),
            canonical,
            connectivity,
            elements,
            "",
        )
    except Exception as error:
        # The coordinate-bearing SDF still parsed. This is normally an RDKit
        # hydrogen/bond-normalization edge case, not evidence of a failed run.
        return "graph_descriptor_unavailable", len(heavy_atoms), "", "", elements, type(error).__name__


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    args = parser.parse_args()

    with args.target_manifest.open(newline="", encoding="utf-8") as handle:
        targets = list(csv.DictReader(handle))

    rows: list[dict[str, object]] = []
    for target in targets:
        if target.get("in_valid_primary_cohort", "").lower() != "true":
            continue
        eq_status, eq_atoms, eq_smiles, eq_connectivity, eq_elements, eq_note = descriptor(target["equibind_primary_path"])
        dd_status, dd_atoms, dd_smiles, dd_connectivity, dd_elements, dd_note = descriptor(target["diffdock_primary_path"])
        readable = eq_status in {"ok", "unsanitized_only", "graph_descriptor_unavailable"} and dd_status in {"ok", "unsanitized_only", "graph_descriptor_unavailable"}
        graph_match = readable and bool(eq_smiles) and eq_smiles == dd_smiles
        connectivity_match = readable and bool(eq_connectivity) and eq_connectivity == dd_connectivity
        element_match = readable and bool(eq_elements) and eq_elements == dd_elements
        if graph_match:
            disposition, note = "eligible_exact_graph", ""
        elif connectivity_match:
            disposition, note = "eligible_stereochemistry_review", "Connectivity agrees; stereochemistry or bond representation differs"
        elif element_match:
            disposition, note = "eligible_representation_review", "Heavy-element composition agrees; verify reference ligand graph before scoring"
        elif not readable:
            disposition, note = "exclude_pending_unreadable_prediction", "; ".join(x for x in (eq_note, dd_note) if x)
        else:
            disposition, note = "exclude_pending_ligand_mismatch", "Predicted ligand heavy-element compositions differ"
        rows.append({
            "target_id": target["target_id"], "equibind_path": target["equibind_primary_path"],
            "diffdock_path": target["diffdock_primary_path"], "equibind_status": eq_status,
            "diffdock_status": dd_status, "equibind_heavy_atoms": eq_atoms,
            "diffdock_heavy_atoms": dd_atoms, "equibind_canonical_smiles": eq_smiles,
            "diffdock_canonical_smiles": dd_smiles,
            "equibind_connectivity_smiles": eq_connectivity,
            "diffdock_connectivity_smiles": dd_connectivity,
            "equibind_element_signature": eq_elements,
            "diffdock_element_signature": dd_elements,
            "graph_match": graph_match, "connectivity_match": connectivity_match,
            "element_signature_match": element_match, "disposition": disposition,
            "note": note,
        })

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "rows": len(rows),
        "dispositions": dict(sorted(Counter(row["disposition"] for row in rows).items())),
        "exact_graph_match": sum(bool(row["graph_match"]) for row in rows),
        "connectivity_match": sum(bool(row["connectivity_match"]) for row in rows),
        "element_signature_match": sum(bool(row["element_signature_match"]) for row in rows),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
