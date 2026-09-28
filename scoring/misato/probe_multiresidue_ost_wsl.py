"""Small WSL validation probe for a graph-verified multi-residue native ligand.

This does not write scores or modify the existing v3/all-copies results.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from audit_crystal_poses import heavy_mol, native_graph, signature, untyped_graph
from audit_multiresidue_recovery import multiresidue_candidates


def build_entity(receptor, candidate, graph, xyz):
    from ost import geom, mol

    entity = receptor.Copy()
    editor = entity.EditXCS(mol.BUFFERED_EDIT)
    chain = editor.InsertChain("MISATO_LIG")
    editor.SetChainType(chain, mol.ChainType.CHAINTYPE_NON_POLY)
    ligand = editor.AppendResidue(chain, "LIG")
    inserted = []
    for index, (atom, position) in enumerate(zip(candidate["atoms"], xyz)):
        inserted.append(editor.InsertAtom(ligand, f"A{index:04d}",
                                          geom.Vec3(*[float(v) for v in position]),
                                          atom.element.name.upper()))
    for bond in graph.GetBonds():
        editor.Connect(inserted[bond.GetBeginAtomIdx()], inserted[bond.GetEndAtomIdx()])
    editor.UpdateICS()
    return entity, ligand


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-id", required=True)
    parser.add_argument("--method", choices=("diffdock", "equibind"), default="diffdock")
    parser.add_argument("--candidate-kind", default="")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root / "scoring" / "plb_bench"))
    import csv

    from ost import mol
    from ost.mol.alg import ligand_scoring_lddtpli, ligand_scoring_scrmsd
    from plb_bench.scoring import _load_mmcif_text
    from score_static_docking_wsl import resolve_path

    with (root / "misato_output" / "misato_score_allcopies_v1.csv").open(newline="") as handle:
        score = next(row for row in csv.DictReader(handle)
                     if row["target_id"] == args.target_id.upper() and row["method"] == args.method)
    pose = heavy_mol(resolve_path(score["pose_path"], root))
    reference = resolve_path(score["reference_path"], root)
    wanted = signature([a.GetAtomicNum() for a in pose.GetAtoms()])
    candidates = [candidate for candidate in multiresidue_candidates(reference, wanted)
                  if not args.candidate_kind or candidate["kind"] == args.candidate_kind]
    for candidate in candidates:
        graph = native_graph(candidate)
        matches = untyped_graph(pose).GetSubstructMatches(untyped_graph(graph),
                                                           uniquify=False, maxMatches=32)
        if matches:
            break
    else:
        raise ValueError("No full-graph candidate in requested category")
    match = matches[0]
    loaded, _, _ = _load_mmcif_text(reference.read_text(encoding="utf-8"))
    receptor = mol.CreateEntityFromView(loaded.Select("ele!=H and ele!=D"), True)
    editor = receptor.EditXCS(mol.BUFFERED_EDIT)
    for chain in list(receptor.chains):
        if ("POLY_PEPTIDE" in str(chain.chain_type) or "POLY_NUCLEOTIDE" in str(chain.chain_type)) \
                and len(chain.residues) > 30:
            continue
        editor.DeleteChain(chain)
    editor.UpdateICS()
    native_xyz = [[atom.pos.x, atom.pos.y, atom.pos.z] for atom in candidate["atoms"]]
    pose_xyz = pose.GetConformer().GetPositions()[list(match)]
    ref_ent, ref_lig = build_entity(receptor, candidate, graph, native_xyz)
    mdl_ent, mdl_lig = build_entity(receptor, candidate, graph, pose_xyz)
    print("candidate", candidate["kind"], candidate["label"], "atoms", len(native_xyz),
          "bonds", graph.GetNumBonds(), "receptor_chains", len(receptor.chains), flush=True)
    for label, scorer_type in (("bisy", ligand_scoring_scrmsd.SCRMSDScorer),
                               ("lddt", ligand_scoring_lddtpli.LDDTPLIScorer)):
        control = scorer_type(model=ref_ent, target=ref_ent, model_ligands=[ref_lig],
                              target_ligands=[ref_lig], substructure_match=False)
        print(label, "self_control", control.score, flush=True)
        scorer = scorer_type(model=mdl_ent, target=ref_ent, model_ligands=[mdl_lig],
                             target_ligands=[ref_lig], substructure_match=False)
        print(label, scorer.score, "assignment", scorer.assignment, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
