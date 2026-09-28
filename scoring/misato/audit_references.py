"""Audit experimental PDB references for the primary MISATO ligand cohort.

The audit is deliberately staged. It validates that PDB IDs resolve and finds
candidate experimental ligand residues by heavy-element signature. It does not
claim that the downloaded PDB is the exact receptor conformation used by the
docking runs, and it does not calculate RMSD.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import tempfile
from collections import Counter
from itertools import combinations_with_replacement
from pathlib import Path

import gemmi
import requests


FIELDS = [
    "target_id", "prediction_disposition", "prediction_element_signature",
    "reference_path", "reference_source", "reference_status",
    "candidate_component_ids", "candidate_residue_instances",
    "candidate_component_type_count", "is_composite_candidate", "audit_disposition", "note",
]
SOLVENT_AND_IONS = {
    "HOH", "WAT", "DOD", "EDO", "GOL", "PEG", "PO4", "SO4", "ACT",
    "MG", "ZN", "CA", "NA", "CL", "K", "MN", "FE", "CU", "CO", "NI",
    "CD", "HG", "PB", "AL", "CS", "RB", "IN", "TL", "SR", "BA", "LI",
    "BR", "F", "I", "NH4", "NO3", "IOD", "FLC",
}


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix=".tmp_", suffix=path.suffix)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
        os.replace(temporary, path)
    except Exception:
        Path(temporary).unlink(missing_ok=True)
        raise


def reference(target_id: str, cache_dir: Path) -> tuple[Path, str]:
    for stem in (f"{target_id.upper()}-assembly1.cif.gz", f"{target_id.upper()}.cif.gz"):
        cached = cache_dir / stem
        if cached.is_file():
            return cached, "cache"
    errors = []
    for stem in (f"{target_id.upper()}-assembly1.cif.gz", f"{target_id.upper()}.cif.gz"):
        url = f"https://files.rcsb.org/download/{stem}"
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                output = cache_dir / stem
                atomic_write(output, response.content)
                return output, "rcsb"
            errors.append(f"{stem}: HTTP {response.status_code}")
        except requests.RequestException as error:
            errors.append(f"{stem}: {type(error).__name__}")
    raise RuntimeError("; ".join(errors))


def signature(residue) -> str:
    return ".".join(str(atom.element.atomic_number) for atom in sorted(residue, key=lambda atom: atom.element.atomic_number) if atom.element.atomic_number > 1)


def combined_signature(signatures: list[str]) -> str:
    elements = []
    for value in signatures:
        elements.extend(int(item) for item in value.split(".") if item)
    return ".".join(str(item) for item in sorted(elements))


def candidates(path: Path, expected_signature: str, max_components: int) -> tuple[list[str], list[str], bool]:
    structure = gemmi.read_structure(str(path))
    if not structure:
        raise ValueError("empty reference structure")
    components: dict[str, dict[str, object]] = {}
    for chain in structure[0]:
        for residue in chain:
            if residue.het_flag != "H" or residue.name.upper() in SOLVENT_AND_IONS:
                continue
            data = components.setdefault(residue.name, {"signature": signature(residue), "instances": []})
            data["instances"].append(f"{chain.name}:{residue.name}:{residue.seqid.num}")
    labels, instances = [], []
    component_ids = sorted(components)
    for count in range(1, max_components + 1):
        for group in combinations_with_replacement(component_ids, count):
            if combined_signature([components[item]["signature"] for item in group]) != expected_signature:
                continue
            label = "+".join(group)
            labels.append(label)
            instances.append("+".join("|".join(components[item]["instances"]) for item in group))
    return labels, instances, any("+" in label for label in labels)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--identity-csv", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0, help="Maximum eligible IDs; 0 means all")
    parser.add_argument("--max-component-combination", type=int, default=3, help="Largest non-polymer component combination to consider")
    args = parser.parse_args()

    with args.identity_csv.open(newline="", encoding="utf-8") as handle:
        inputs = [row for row in csv.DictReader(handle) if row["disposition"].startswith("eligible_")]
    if args.limit:
        inputs = inputs[:args.limit]
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for row in inputs:
        target_id = row["target_id"].upper()
        output = {
            "target_id": target_id, "prediction_disposition": row["disposition"],
            "prediction_element_signature": row["equibind_element_signature"],
            "reference_path": "", "reference_source": "", "reference_status": "",
            "candidate_component_ids": "", "candidate_residue_instances": "",
            "candidate_component_type_count": 0, "is_composite_candidate": False,
            "audit_disposition": "", "note": "",
        }
        try:
            path, source = reference(target_id, args.cache_dir)
            component_ids, instances, composite = candidates(path, row["equibind_element_signature"], args.max_component_combination)
            output.update({"reference_path": str(path), "reference_source": source, "reference_status": "available", "candidate_component_ids": ";".join(component_ids), "candidate_residue_instances": ";".join(instances), "candidate_component_type_count": len(component_ids), "is_composite_candidate": composite})
            if len(component_ids) == 1:
                output["audit_disposition"] = "composite_reference_ligand_requires_graph_confirmation" if composite else "candidate_reference_ligand_identified"
            elif len(component_ids) > 1:
                output["audit_disposition"] = "ambiguous_reference_ligand"
            else:
                output["audit_disposition"] = "no_signature_matched_reference_ligand"
        except Exception as error:
            output.update({"reference_status": "unavailable_or_unreadable", "audit_disposition": "reference_review_required", "note": f"{type(error).__name__}: {error}"})
        rows.append(output)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"rows": len(rows), "dispositions": dict(sorted(Counter(row["audit_disposition"] for row in rows).items()))}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
