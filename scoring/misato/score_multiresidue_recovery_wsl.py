"""Score graph-verified MISATO multi-residue ligands without changing v3.

Each native crystal assembly becomes one generic pseudo-residue. Atom names
are keyed to the QM pose atom order (not the duplicated original residue atom
names), and bonds follow the verified full heavy-atom graph. All graph-valid
crystal copies are offered to OST for assignment.
"""
from __future__ import annotations

import argparse
import csv
import math
import signal
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

try:
    from .audit_crystal_poses import heavy_mol, native_graph, signature, untyped_graph
    from .audit_multiresidue_recovery import multiresidue_candidates
    from .score_static_docking_wsl import initialize_worker, resolve_path
except ImportError:  # direct WSL script execution
    from audit_crystal_poses import heavy_mol, native_graph, signature, untyped_graph
    from audit_multiresidue_recovery import multiresidue_candidates
    from score_static_docking_wsl import initialize_worker, resolve_path


FIELDS = ["target_id", "method", "status", "bisy_rmsd_angstrom", "lddt_pli",
          "candidate_kind", "n_target_copies", "assigned_copy_label",
          "elapsed_seconds", "error"]


def _assigned(scorer):
    if not scorer.assignment:
        return None, None
    target_index, model_index = scorer.assignment[0]
    try:
        value = float(scorer.score_matrix[target_index, model_index])
    except (TypeError, ValueError):
        return None, None
    return (None, None) if not math.isfinite(value) else (value, target_index)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def graph_valid_copies(pose, candidates: list[dict]) -> list[dict]:
    """Return chemically matching native assemblies, deduplicated by atoms."""
    pose_graph = untyped_graph(pose)
    valid = []
    seen_atoms = set()
    for candidate in candidates:
        try:
            native = native_graph(candidate)
            if native.GetNumBonds() != pose.GetNumBonds():
                continue
            mapping = pose_graph.GetSubstructMatch(untyped_graph(native), useChirality=False)
        except (ValueError, RuntimeError):
            continue
        if len(mapping) != pose.GetNumAtoms():
            continue
        coordinates = tuple(sorted((atom.element.name.upper(), round(atom.pos.x, 3),
                                    round(atom.pos.y, 3), round(atom.pos.z, 3))
                                   for atom in candidate["atoms"]))
        if coordinates in seen_atoms:
            continue
        seen_atoms.add(coordinates)
        native_by_pose = [None] * len(mapping)
        for native_index, pose_index in enumerate(mapping):
            atom = candidate["atoms"][native_index]
            native_by_pose[pose_index] = (atom.pos.x, atom.pos.y, atom.pos.z)
        valid.append({**candidate, "native_xyz_by_pose": native_by_pose})
    return valid


def same_indexed_heavy_graph(first, second) -> bool:
    if first.GetNumAtoms() != second.GetNumAtoms():
        return False
    if [atom.GetAtomicNum() for atom in first.GetAtoms()] != [
            atom.GetAtomicNum() for atom in second.GetAtoms()]:
        return False
    def edges(mol):
        return {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
                for bond in mol.GetBonds()}
    return edges(first) == edges(second)


def receptor_without_ligands(reference_text: str, ligand_subchains: set[str]):
    from ost import mol
    from plb_bench.scoring import _load_mmcif_text

    loaded, _, _ = _load_mmcif_text(reference_text)
    receptor = mol.CreateEntityFromView(loaded.Select("ele!=H and ele!=D"), True)
    editor = receptor.EditXCS(mol.BUFFERED_EDIT)
    for chain in list(receptor.chains):
        kind = str(chain.chain_type)
        is_receptor_polymer = "POLY" in kind and "NON_POLY" not in kind \
            and "OLIGOSACCHARIDE" not in kind
        if chain.name in ligand_subchains or not is_receptor_polymer:
            editor.DeleteChain(chain)
    editor.UpdateICS()
    if not receptor.chains:
        raise ValueError("no_receptor_polymer_after_ligand_removal")
    return receptor


def add_pseudo(editor, chain_name: str, elements: list[str], xyz, graph):
    from ost import geom, mol

    chain = editor.InsertChain(chain_name)
    editor.SetChainType(chain, mol.ChainType.CHAINTYPE_NON_POLY)
    residue = editor.AppendResidue(chain, "LIG")
    atoms = [editor.InsertAtom(residue, f"A{index:04d}",
                               geom.Vec3(*[float(v) for v in point]), element)
             for index, (element, point) in enumerate(zip(elements, xyz))]
    for bond in graph.GetBonds():
        editor.Connect(atoms[bond.GetBeginAtomIdx()], atoms[bond.GetEndAtomIdx()])
    return residue


def build_complexes(receptor, pose, copies: list[dict]):
    from ost import mol

    graph = untyped_graph(pose)
    elements = [atom.GetSymbol().upper() for atom in pose.GetAtoms()]
    target = receptor.Copy()
    editor = target.EditXCS(mol.BUFFERED_EDIT)
    labels = []
    for index, candidate in enumerate(copies):
        add_pseudo(editor, f"MISATO_REF_{index}", elements,
                   candidate["native_xyz_by_pose"], graph)
        labels.append(candidate["label"])
    editor.UpdateICS()
    target_ligands = [target.FindChain(f"MISATO_REF_{index}").residues[0]
                      for index in range(len(copies))]
    model = receptor.Copy()
    editor = model.EditXCS(mol.BUFFERED_EDIT)
    add_pseudo(editor, "MISATO_MODEL", elements,
               pose.GetConformer().GetPositions(), graph)
    editor.UpdateICS()
    model_ligand = model.FindChain("MISATO_MODEL").residues[0]
    return model, target, model_ligand, target_ligands, labels


def score_target(target_id: str, rows: list[dict[str, str]], root_text: str,
                 crystal_dir_text: str, timeout_seconds: int) -> list[dict[str, object]]:
    root = Path(root_text)
    start = time.perf_counter()
    from ost.mol.alg import ligand_scoring_lddtpli, ligand_scoring_scrmsd

    reference = Path(crystal_dir_text) / f"{target_id[:4]}.cif"
    result = []
    try:
        signal.alarm(timeout_seconds)
        anchor = heavy_mol(resolve_path(rows[0]["pose_path"], root))
        wanted = signature([atom.GetAtomicNum() for atom in anchor.GetAtoms()])
        all_candidates = multiresidue_candidates(reference, wanted)
        copies = graph_valid_copies(anchor, all_candidates)
        if not copies:
            raise ValueError("no_full_graph_native_assembly")
        if len(copies) > 32:
            raise ValueError(f"too_many_native_copies_{len(copies)}")
        ligand_subchains = {name for candidate in all_candidates
                            for name in candidate["subchains"]}
        receptor = receptor_without_ligands(reference.read_text(encoding="utf-8"),
                                            ligand_subchains)
        signal.alarm(0)
    except Exception as error:  # noqa: BLE001 - retain per-target failure provenance
        signal.alarm(0)
        return [{"target_id": target_id, "method": row["method"],
                 "status": "failed_reference_preparation", "error": f"{type(error).__name__}: {error}",
                 "elapsed_seconds": round(time.perf_counter() - start, 3)} for row in rows]
    for row in rows:
        item: dict[str, object] = {"target_id": target_id, "method": row["method"],
                                   "status": "pending", "bisy_rmsd_angstrom": "",
                                   "lddt_pli": "", "candidate_kind": ";".join(sorted({c["kind"] for c in copies})),
                                   "n_target_copies": len(copies), "assigned_copy_label": "",
                                   "elapsed_seconds": "", "error": ""}
        try:
            signal.alarm(timeout_seconds)
            pose = heavy_mol(resolve_path(row["pose_path"], root))
            if not same_indexed_heavy_graph(pose, anchor):
                raise ValueError("method_pose_atom_order_or_graph_differs_from_anchor")
            model, target, model_ligand, target_ligands, labels = build_complexes(
                receptor, pose, copies)
            assignments = {}
            for key, scorer_type in (("bisy_rmsd_angstrom", ligand_scoring_scrmsd.SCRMSDScorer),
                                     ("lddt_pli", ligand_scoring_lddtpli.LDDTPLIScorer)):
                scorer = scorer_type(model=model, target=target,
                                     model_ligands=[model_ligand], target_ligands=target_ligands,
                                     substructure_match=False)
                assignments[key] = _assigned(scorer)
            for key, (value, _) in assignments.items():
                item[key] = "" if value is None else value
            assigned_index = assignments["bisy_rmsd_angstrom"][1]
            if assigned_index is not None:
                item["assigned_copy_label"] = labels[assigned_index]
            item["status"] = ("scored_both" if item["bisy_rmsd_angstrom"] != ""
                              and item["lddt_pli"] != "" else "partial_or_no_metrics")
            if item["status"] != "scored_both":
                item["error"] = "OST returned one or no metrics"
        except Exception as error:  # noqa: BLE001 - retain per-pose failure provenance
            item["status"] = "failed_runtime_timeout" if isinstance(error, TimeoutError) else "failed_scoring"
            item["error"] = f"{type(error).__name__}: {error}"
        finally:
            signal.alarm(0)
        item["elapsed_seconds"] = round(time.perf_counter() - start, 3)
        result.append(item)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path,
                        default=Path("misato_output/multiresidue_recovery_all_v1.csv"))
    parser.add_argument("--allcopies", type=Path,
                        default=Path("misato_output/misato_score_allcopies_v1.csv"))
    parser.add_argument("--crystal-dir", type=Path,
                        default=Path("misato_output/rcsb_asymmetric_unit"))
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=int, default=120)
    parser.add_argument("--target-id", action="append", default=[])
    parser.add_argument("--skip-target-id", action="append", default=[],
                        help="Record a known pathological target as a runtime exclusion")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.timeout_seconds < 1:
        parser.error("workers and timeout must be positive")
    if args.out_csv.exists() and not args.resume:
        parser.error("Output already exists; choose a new file or --resume")
    if args.resume and not args.out_csv.exists():
        parser.error("--resume requires an existing CSV")
    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root / "scoring" / "plb_bench"))
    eligible = {(row["target_id"], row["method"]) for row in read_rows(args.audit)
                if row["graph_status"] == "exact_full_graph"}
    requested = {target.upper() for target in args.target_id}
    if requested:
        eligible = {key for key in eligible if key[0] in requested}
    prior = {(row["target_id"], row["method"]) for row in read_rows(args.out_csv)} \
        if args.resume else set()
    skip = {target.upper() for target in args.skip_target_id}
    by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    skipped = []
    for row in read_rows(args.allcopies):
        key = (row["target_id"], row["method"])
        if key in eligible and key not in prior:
            if row["target_id"] in skip:
                skipped.append({"target_id": row["target_id"], "method": row["method"],
                                "status": "excluded_runtime_guard", "bisy_rmsd_angstrom": "",
                                "lddt_pli": "", "error":
                                "Predeclared hard-case exclusion after >8-minute OST stall"})
            else:
                by_id[row["target_id"]].append(row)
    if sum(map(len, by_id.values())) + len(skipped) != len(eligible - prior):
        parser.error("Eligible audit rows absent from all-copies score input")
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    with args.out_csv.open("a" if args.resume else "x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if not args.resume:
            writer.writeheader()
        writer.writerows(skipped)
        counts.update(str(row["status"]) for row in skipped)
        handle.flush()
        with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize_worker) as pool:
            futures = {pool.submit(score_target, target_id, rows, str(root),
                                   str(args.crystal_dir.resolve()), args.timeout_seconds): target_id
                       for target_id, rows in by_id.items()}
            for future in as_completed(futures):
                target_id = futures[future]
                try:
                    completed = future.result()
                except Exception as error:  # noqa: BLE001 - retain worker failure provenance
                    completed = [{"target_id": target_id, "method": row["method"],
                                  "status": "failed_worker", "error": f"{type(error).__name__}: {error}"}
                                 for row in by_id[target_id]]
                writer.writerows(completed)
                counts.update(str(row["status"]) for row in completed)
                handle.flush()
    print(f"Wrote {sum(counts.values())} recovery rows: {dict(counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
