"""OST ligand metrics constrained to one independently selected native residue.

The generic scorer optimizes over *all* nonpolymer ligands in a crystal CIF.
That is unsafe for MISATO static docking: partial copies and cofactors can
yield an attractive but incorrect score. This adapter uses only the one
pre-audited native ligand and disallows substructure-only matches.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path


def load_pose(path: Path):
    from rdkit import Chem

    for sanitize in (True, False):
        supplier = Chem.SDMolSupplier(str(path), removeHs=True, sanitize=sanitize)
        mol = next((item for item in supplier if item is not None and item.GetNumConformers()), None)
        if mol is not None:
            return mol
    raise ValueError("pose_sdf_unreadable_by_rdkit")


def write_heavy_only_pose(source: Path, destination: Path) -> int:
    """Write a temporary heavy-atom SDF for OST's full-graph matching."""
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(source), removeHs=False, sanitize=False)
    mol = next((item for item in supplier if item is not None and item.GetNumConformers()), None)
    if mol is None:
        raise ValueError("pose_sdf_unreadable_by_rdkit")
    mol.UpdatePropertyCache(strict=False)
    heavy = Chem.RemoveHs(mol, sanitize=False)
    if heavy.GetNumAtoms() < 2:
        raise ValueError("too_few_heavy_atoms")
    writer = Chem.SDWriter(str(destination))
    try:
        writer.write(heavy)
    finally:
        writer.close()
    return heavy.GetNumAtoms()


def choose_native(candidates, expected_name: str, pose_mol):
    wanted = Counter(atom.GetSymbol().upper() for atom in pose_mol.GetAtoms()
                     if atom.GetAtomicNum() > 1)
    matches = []
    for residue in candidates:
        if residue.name.upper() != expected_name.upper():
            continue
        elements = Counter(str(atom.element).upper() for atom in residue.atoms
                           if str(atom.element).upper() not in {"H", "D"})
        if elements == wanted:
            matches.append(residue)
    if len(matches) != 1:
        raise ValueError(f"selected_native_residue_count_{len(matches)}")
    return matches[0]


def _extract(value, metric: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, dict):
        for key in (("rmsd", "score", "bisy_rmsd") if metric == "rmsd"
                    else ("lddt_pli", "lddt", "score")):
            if key in value and value[key] is not None:
                try:
                    return float(value[key])
                except (TypeError, ValueError):
                    pass
        for item in value.values():
            try:
                return float(item)
            except (TypeError, ValueError):
                pass
    for key in (("rmsd", "bisy_rmsd", "score", "value") if metric == "rmsd"
                else ("lddt_pli", "lddt", "score", "value")):
        try:
            item = getattr(value, key, None)
            if item is not None:
                return float(item)
        except (TypeError, ValueError):
            pass
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def score_selected(model_ent, ref_ent, pose_path: Path,
                   expected_residue_name: str) -> dict[str, object]:
    from ost.mol.alg import ligand_scoring_lddtpli, ligand_scoring_scrmsd
    from plb_bench.scoring import _get_ligands

    pose_mol = load_pose(pose_path)
    native = choose_native(_get_ligands(ref_ent), expected_residue_name, pose_mol)
    result: dict[str, object] = {
        "bisy_rmsd": None, "lddt_pli": None,
        "selected_ref_ost_chain": native.chain.name,
        "selected_ref_ost_atom_count": len(native.atoms),
    }
    model_candidates = [[pose_mol]]
    handles = _get_ligands(model_ent)
    if handles:
        model_candidates.append(handles)
    for metric, scorer_type in (("bisy_rmsd", ligand_scoring_scrmsd.SCRMSDScorer),
                                ("lddt_pli", ligand_scoring_lddtpli.LDDTPLIScorer)):
        for model_ligands in model_candidates:
            try:
                scorer = scorer_type(model=model_ent, target=ref_ent,
                                     model_ligands=model_ligands,
                                     target_ligands=[native],
                                     substructure_match=False)
                values = [_extract(value, "rmsd" if metric == "bisy_rmsd" else "lddt")
                          for value in (scorer.score or {}).values()]
                values = [value for value in values if value is not None]
                if values:
                    result[metric] = min(values) if metric == "bisy_rmsd" else max(values)
                    break
            except Exception:
                continue
    return result
