"""Conservatively audit MISATO SDF poses against static RCSB crystal ligands.

This is a *provisional* direct-coordinate pose audit, not the final OST
BiSyRMSD/lDDT-PLI benchmark. It never modifies source SDFs or crystal CIFs.
Ambiguous and composite native ligands are reported rather than force-scored.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter
from itertools import combinations
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds


RDLogger.DisableLog("rdApp.*")
EXCLUDE = {"HOH", "WAT", "DOD", "EDO", "GOL", "PEG", "PO4", "SO4", "ACT",
           "MG", "ZN", "CA", "NA", "CL", "K", "MN", "FE", "CU", "CO", "NI",
           "CD", "HG", "PB", "AL", "CS", "RB", "IN", "TL", "SR", "BA", "LI",
           "BR", "F", "I", "NH4", "NO3", "IOD", "FLC"}
FIELDS = ["target_id", "pdb_id", "method", "pose_path", "reference_path",
          "reference_status", "reference_candidate_count", "reference_selection",
          "reference_residues", "pose_heavy_atoms", "native_heavy_atoms",
          "pose_to_native_centroid_angstrom", "pose_rmsd_angstrom_provisional",
          "graph_mapping_count", "note"]


def heavy_mol(path: Path):
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    mol = next((item for item in supplier if item is not None and item.GetNumConformers()), None)
    if mol is None:
        raise ValueError("unreadable_sdf")
    heavy = Chem.RemoveHs(mol, sanitize=False)
    if heavy.GetNumAtoms() < 2:
        raise ValueError("too_few_heavy_atoms")
    return heavy


def signature(elements: list[int]) -> Counter[int]:
    return Counter(number for number in elements if number > 1)


def reference_residues(cif_path: Path) -> list[dict]:
    structure = gemmi.read_structure(str(cif_path))
    if not structure:
        raise ValueError("empty_crystal_structure")
    out = []
    for chain in structure[0]:
        for residue in chain:
            if residue.het_flag != "H" or residue.name.upper() in EXCLUDE:
                continue
            # One conformer per atom name; prefer blank/A altloc where present.
            atoms_by_name = {}
            for atom in residue:
                if atom.element.atomic_number <= 1:
                    continue
                name = atom.name.strip()
                altloc = str(atom.altloc).strip("\x00 ")
                if name not in atoms_by_name or altloc in {"", "A"}:
                    atoms_by_name[name] = atom
            atoms = list(atoms_by_name.values())
            if len(atoms) < 2:
                continue
            coords = np.array([[atom.pos.x, atom.pos.y, atom.pos.z] for atom in atoms], dtype=float)
            label = f"{chain.name}:{residue.name}:{residue.seqid.num}{residue.seqid.icode.strip()}"
            out.append({"label": label, "atoms": atoms, "coords": coords,
                        "signature": signature([atom.element.atomic_number for atom in atoms])})
    return out


def native_graph(candidate: dict):
    atoms = candidate["atoms"]
    lines = [str(len(atoms)), "native ligand"]
    lines.extend(f"{atom.element.name} {atom.pos.x:.5f} {atom.pos.y:.5f} {atom.pos.z:.5f}" for atom in atoms)
    mol = Chem.MolFromXYZBlock("\n".join(lines) + "\n")
    if mol is None:
        raise ValueError("native_xyz_unreadable")
    rdDetermineBonds.DetermineConnectivity(mol)
    return mol


def untyped_graph(mol):
    graph = Chem.RWMol(mol)
    for atom in graph.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetIsAromatic(False)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in graph.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetStereo(Chem.BondStereo.STEREONONE)
    result = graph.GetMol()
    result.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(result)
    return result


def direct_pose_rmsd(pose, native, max_matches: int = 4096) -> tuple[float | None, int, str]:
    if pose.GetNumAtoms() != native.GetNumAtoms():
        return None, 0, "atom_count_mismatch"
    pose_graph = untyped_graph(pose)
    native_graph = untyped_graph(native)
    matches = pose_graph.GetSubstructMatches(native_graph, uniquify=False,
                                             useChirality=False, maxMatches=max_matches)
    if len(matches) >= max_matches:
        return None, len(matches), "graph_mapping_limit_reached"
    if not matches:
        return None, 0, "no_graph_isomorphism"
    pose_coords = np.asarray(pose.GetConformer().GetPositions())
    native_coords = np.asarray(native.GetConformer().GetPositions())
    rmsd = min(math.sqrt(np.mean(np.sum((pose_coords[list(mapping)] - native_coords) ** 2, axis=1)))
               for mapping in matches)
    return rmsd, len(matches), ""


def crystal_candidates(residues: list[dict], wanted: Counter[int]) -> tuple[list[dict], bool]:
    singles = [residue for residue in residues if residue["signature"] == wanted]
    if singles:
        return singles, False
    eligible = [residue for residue in residues if all(count <= wanted[number]
                for number, count in residue["signature"].items())]
    # Two-component assemblies are identified but not force-scored here.
    pairs = [{"label": f"{a['label']}+{b['label']}"} for a, b in combinations(eligible, 2)
             if a["signature"] + b["signature"] == wanted]
    return pairs, bool(pairs)


def audit_target(target: dict, crystal_dir: Path) -> list[dict[str, object]]:
    target_id = target["target_id"].strip().upper()
    match = re.match(r"^[A-Z0-9]{4}", target_id)
    if not match:
        raise ValueError(f"Invalid target ID: {target_id}")
    pdb_id = match.group(0)
    crystal_path = crystal_dir / f"{pdb_id}.cif"
    methods = [(method, target[f"{method}_primary_path"])
               for method in ("equibind", "diffdock") if target[f"{method}_status"] == "candidate"]
    rows = []
    poses = {}
    for method, path_text in methods:
        try:
            poses[method] = heavy_mol(Path(path_text))
        except (OSError, ValueError) as error:
            poses[method] = error
    if crystal_path.is_file():
        try:
            residues = reference_residues(crystal_path)
            crystal_error = ""
        except (OSError, ValueError, RuntimeError) as error:
            residues, crystal_error = [], f"{type(error).__name__}: {error}"
    else:
        residues, crystal_error = [], "crystal_cif_not_cached"

    # Use one native copy for both methods. DiffDock is the anchor when present;
    # selecting a different copy separately for each method would bias scores.
    anchor = poses.get("diffdock")
    if not hasattr(anchor, "GetAtoms"):
        anchor = poses.get("equibind")
    wanted = signature([atom.GetAtomicNum() for atom in anchor.GetAtoms()]) if hasattr(anchor, "GetAtoms") else Counter()
    candidates, composite = crystal_candidates(residues, wanted) if wanted else ([], False)
    chosen = None
    selection = ""
    if candidates and not composite:
        if len(candidates) == 1:
            chosen, selection = candidates[0], "single_candidate"
        elif hasattr(anchor, "GetConformer"):
            anchor_center = np.asarray(anchor.GetConformer().GetPositions()).mean(axis=0)
            ranked = sorted(((float(np.linalg.norm(anchor_center - item["coords"].mean(axis=0))), item)
                             for item in candidates), key=lambda pair: pair[0])
            if ranked[1][0] - ranked[0][0] >= 5.0:
                chosen, selection = ranked[0][1], "copy_by_diffdock_proximity" if "diffdock" in poses else "copy_by_equibind_proximity"
            else:
                selection = "ambiguous_crystal_copies"
    elif composite:
        selection = "composite_reference_review"
    elif crystal_error:
        selection = "crystal_unavailable"
    else:
        selection = "no_native_element_match"
    native = None
    if chosen:
        try:
            native = native_graph(chosen)
        except (ValueError, RuntimeError) as error:
            selection = "native_graph_failed"
            crystal_error = f"{type(error).__name__}: {error}"
    for method, path_text in methods:
        pose = poses[method]
        row = {"target_id": target_id, "pdb_id": pdb_id, "method": method,
               "pose_path": path_text, "reference_path": str(crystal_path) if crystal_path.is_file() else "",
               "reference_status": selection, "reference_candidate_count": len(candidates),
               "reference_selection": selection, "reference_residues": chosen["label"] if chosen else "",
               "pose_heavy_atoms": pose.GetNumAtoms() if hasattr(pose, "GetNumAtoms") else "",
               "native_heavy_atoms": native.GetNumAtoms() if native else "",
               "pose_to_native_centroid_angstrom": "", "pose_rmsd_angstrom_provisional": "",
               "graph_mapping_count": "", "note": crystal_error}
        if not hasattr(pose, "GetConformer"):
            row["reference_status"] = "pose_unreadable"
            row["note"] = str(pose)
        elif native:
            row["pose_to_native_centroid_angstrom"] = round(float(np.linalg.norm(
                np.asarray(pose.GetConformer().GetPositions()).mean(axis=0) - chosen["coords"].mean(axis=0))), 4)
            rmsd, mappings, issue = direct_pose_rmsd(pose, native)
            row["graph_mapping_count"] = mappings
            if rmsd is not None:
                row["pose_rmsd_angstrom_provisional"] = round(rmsd, 4)
                row["reference_status"] = "provisional_pose_rmsd"
            else:
                row["reference_status"] = issue
                row["note"] = issue
        rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", type=Path, default=Path("misato_output/inventory/target_manifest.csv"))
    parser.add_argument("--crystal-dir", type=Path, default=Path("misato_output/rcsb_asymmetric_unit"))
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0, help="First N targets; 0 means all")
    args = parser.parse_args()
    if args.limit < 0 or args.out_csv.exists():
        parser.error("--limit must be nonnegative and --out-csv must not already exist")
    with args.target_manifest.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        targets = [row for row in reader if row["equibind_status"] == "candidate"
                   or row["diffdock_status"] == "candidate"]
    if args.limit:
        targets = targets[:args.limit]
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    statuses: Counter[str] = Counter()
    with args.out_csv.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for target in targets:
            for row in audit_target(target, args.crystal_dir):
                writer.writerow(row)
                statuses[str(row["reference_status"])] += 1
            handle.flush()
    print(f"Targets: {len(targets)}; rows: {sum(statuses.values())}; statuses: {dict(statuses)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
