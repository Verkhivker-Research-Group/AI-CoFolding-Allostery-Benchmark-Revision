"""Static-receptor, ligand-inclusive contact QS for delivered MISATO poses.

Unlike co-folding QS, this compares two ligand placements in the *same*
crystal-polymer receptor. It uses the contact-distance weighting from
``scoring/fix_qs_dynamicbind.py`` but neither receptor superposition nor
nearest-residue remapping is appropriate here.
"""
from __future__ import annotations

import argparse
import csv
import math
import signal
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

try:
    from .audit_crystal_poses import heavy_mol, signature
    from .audit_multiresidue_recovery import multiresidue_candidates
    from .score_static_docking_wsl import initialize_worker, resolve_path
except ImportError:
    from audit_crystal_poses import heavy_mol, signature
    from audit_multiresidue_recovery import multiresidue_candidates
    from score_static_docking_wsl import initialize_worker, resolve_path


CONTACT_D = 12.0  # Angstrom, matching fix_qs_dynamicbind.py
FIELDS = ["target_id", "method", "status", "static_receptor_qs", "contact_cutoff_angstrom",
          "native_contact_residues", "pose_contact_residues", "assigned_crystal_copy", "error"]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def residue_min_distances(receptor_xyz: np.ndarray, residue_indices: np.ndarray,
                          n_residues: int, ligand_xyz: np.ndarray,
                          contact_d: float = CONTACT_D) -> np.ndarray:
    """Minimum heavy-atom distance from each receptor residue to a ligand."""
    if not len(receptor_xyz) or not len(ligand_xyz):
        raise ValueError("empty_receptor_or_ligand")
    nearest, _ = cKDTree(ligand_xyz).query(receptor_xyz, distance_upper_bound=contact_d)
    distances = np.full(n_residues, math.inf, dtype=float)
    np.minimum.at(distances, residue_indices, nearest)
    return distances


def contact_qs(ref_distances: np.ndarray, pose_distances: np.ndarray,
               contact_d: float = CONTACT_D) -> float | None:
    """The ligand-inclusive QS formula from fix_qs_dynamicbind.py."""
    if ref_distances.shape != pose_distances.shape or contact_d <= 0:
        raise ValueError("incompatible_distance_vectors_or_cutoff")
    ref_contact = ref_distances <= contact_d
    pose_contact = pose_distances <= contact_d
    union = ref_contact | pose_contact
    if not union.any():
        return None
    shared = ref_contact & pose_contact
    shared_w = float(np.maximum(0.0, 1.0 -
                     np.abs(ref_distances[shared] - pose_distances[shared]) / contact_d).sum())
    nonshared = int(np.count_nonzero(union & ~shared))
    denominator = shared_w + nonshared
    return shared_w / denominator if denominator > 0 else None


def receptor_atoms(reference, ligand_subchains: set[str]) -> tuple[np.ndarray, np.ndarray, int]:
    """Use the crystal polymer, excluding the identified multi-residue ligand."""
    xyz = []
    indices = []
    n_residues = 0
    for chain in reference.chains:
        kind = str(chain.chain_type)
        keep = ("POLY" in kind and "NON_POLY" not in kind
                and "OLIGOSACCHARIDE" not in kind and chain.name not in ligand_subchains)
        if not keep:
            continue
        for residue in chain.residues:
            atoms = [atom for atom in residue.atoms
                     if str(atom.element).upper() not in {"H", "D"}]
            if not atoms:
                continue
            for atom in atoms:
                xyz.append((atom.pos.x, atom.pos.y, atom.pos.z))
                indices.append(n_residues)
            n_residues += 1
    if not n_residues:
        raise ValueError("no_receptor_polymer_residues")
    return np.asarray(xyz, dtype=float), np.asarray(indices, dtype=int), n_residues


def single_native_xyz(reference, pose, row: dict[str, str],
                      old: dict[str, str], broad: dict[str, str]) -> np.ndarray:
    from plb_bench.scoring import _get_ligands

    wanted = Counter(atom.GetSymbol().upper() for atom in pose.GetAtoms())
    source = row["metric_source"]
    if source == "v3_strict_original":
        name = old["reference_residue"].split(":")[1].upper()
        chain_name = old["selected_ref_ost_chain"]
    else:
        name = broad["reference_resname"].upper()
        chain_name = (broad["assigned_ref_ost_chain"] if source == "allcopies_v1"
                      else row["assigned_crystal_copy"])
    candidates = []
    for residue in _get_ligands(reference):
        elements = Counter(str(atom.element).upper() for atom in residue.atoms
                           if str(atom.element).upper() not in {"H", "D"})
        if residue.name.upper() == name and residue.chain.name == chain_name and elements == wanted:
            candidates.append(residue)
    if len(candidates) != 1:
        raise ValueError(f"assigned_single_native_copy_count_{len(candidates)}")
    return np.asarray([(atom.pos.x, atom.pos.y, atom.pos.z)
                       for atom in candidates[0].atoms
                       if str(atom.element).upper() not in {"H", "D"}], dtype=float)


def score_target(target_id: str, rows: list[dict[str, str]], root_text: str,
                 crystal_dir_text: str, old_rows: dict[str, dict[str, str]],
                 broad_rows: dict[str, dict[str, str]], contact_d: float,
                 timeout_seconds: int) -> list[dict[str, object]]:
    from ost import mol
    from plb_bench.scoring import _load_mmcif_text

    root = Path(root_text)
    result = []
    try:
        signal.alarm(timeout_seconds)
        cif = Path(crystal_dir_text) / f"{target_id[:4]}.cif"
        reference_text = cif.read_text(encoding="utf-8")
        loaded, _, _ = _load_mmcif_text(reference_text)
        reference = mol.CreateEntityFromView(loaded.Select("ele!=H and ele!=D"), True)
        signal.alarm(0)
    except Exception as error:  # noqa: BLE001 - preserve a failure row for each delivered pose
        signal.alarm(0)
        return [{"target_id": target_id, "method": row["method"],
                 "status": "reference_failed", "error": f"{type(error).__name__}: {error}"}
                for row in rows]
    for row in rows:
        item: dict[str, object] = {field: "" for field in FIELDS}
        item.update(target_id=target_id, method=row["method"], status="pending",
                    contact_cutoff_angstrom=contact_d,
                    assigned_crystal_copy=row["assigned_crystal_copy"])
        try:
            signal.alarm(timeout_seconds)
            method = row["method"]
            pose_path = resolve_path(broad_rows[method]["pose_path"], root)
            pose = heavy_mol(pose_path)
            pose_xyz = np.asarray(pose.GetConformer().GetPositions(), dtype=float)
            source = row["metric_source"]
            ligand_subchains: set[str] = set()
            if source.startswith("ost_pseudo_multiresidue"):
                wanted = signature([atom.GetAtomicNum() for atom in pose.GetAtoms()])
                candidates = multiresidue_candidates(cif, wanted)
                matches = [c for c in candidates if c["label"] == row["assigned_crystal_copy"]]
                if len(matches) != 1:
                    raise ValueError(f"assigned_multiresidue_copy_count_{len(matches)}")
                ligand_subchains = {subchain for c in candidates for subchain in c["subchains"]}
                native_xyz = np.asarray([(atom.pos.x, atom.pos.y, atom.pos.z)
                                         for atom in matches[0]["atoms"]], dtype=float)
            else:
                native_xyz = single_native_xyz(reference, pose, row, old_rows[method],
                                               broad_rows[method])
            receptor_xyz, residue_indices, n_residues = receptor_atoms(reference,
                                                                       ligand_subchains)
            native_dist = residue_min_distances(receptor_xyz, residue_indices, n_residues,
                                                 native_xyz, contact_d)
            pose_dist = residue_min_distances(receptor_xyz, residue_indices, n_residues,
                                               pose_xyz, contact_d)
            control = contact_qs(native_dist, native_dist, contact_d)
            if control is None or abs(control - 1.0) > 1e-12:
                raise ValueError("native_self_control_not_one")
            outside = residue_min_distances(receptor_xyz, residue_indices, n_residues,
                                            native_xyz + 10000.0, contact_d)
            if contact_qs(native_dist, outside, contact_d) != 0.0:
                raise ValueError("out_of_pocket_control_not_zero")
            value = contact_qs(native_dist, pose_dist, contact_d)
            if value is None:
                raise ValueError("no_native_or_pose_contacts")
            item["static_receptor_qs"] = value
            item["native_contact_residues"] = int(np.count_nonzero(native_dist <= contact_d))
            item["pose_contact_residues"] = int(np.count_nonzero(pose_dist <= contact_d))
            item["status"] = "scored"
        except Exception as error:  # noqa: BLE001 - retain per-pose failure provenance
            item["status"] = "failed"
            item["error"] = f"{type(error).__name__}: {error}"
        finally:
            signal.alarm(0)
        result.append(item)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path,
                        default=Path("misato_output/misato_score_recovered_v1.csv"))
    parser.add_argument("--v3", type=Path,
                        default=Path("misato_output/misato_score_all_heavy_selected_v3.csv"))
    parser.add_argument("--allcopies", type=Path,
                        default=Path("misato_output/misato_score_allcopies_v1.csv"))
    parser.add_argument("--crystal-dir", type=Path,
                        default=Path("misato_output/rcsb_asymmetric_unit"))
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--timeout-seconds", type=int, default=120)
    parser.add_argument("--contact-d", type=float, default=CONTACT_D)
    parser.add_argument("--limit-targets", type=int, default=0)
    parser.add_argument("--target-id", action="append", default=[])
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.timeout_seconds < 1 or args.contact_d <= 0:
        parser.error("workers, timeout, and contact distance must be positive")
    if args.out_csv.exists() and not args.resume:
        parser.error("Output exists; choose a new name or --resume")
    if args.resume and not args.out_csv.exists():
        parser.error("--resume requires an existing CSV")
    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root / "scoring" / "plb_bench"))
    by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    wanted_ids = {target.upper() for target in args.target_id}
    for row in read_rows(args.scores):
        if row["score_status"] == "scored_both" and (not wanted_ids or row["target_id"] in wanted_ids):
            by_id[row["target_id"]].append(row)
    if args.limit_targets:
        by_id = {key: by_id[key] for key in sorted(by_id)[:args.limit_targets]}
    prior = {(row["target_id"], row["method"]) for row in read_rows(args.out_csv)} \
        if args.resume else set()
    for target_id in list(by_id):
        by_id[target_id] = [row for row in by_id[target_id]
                            if (target_id, row["method"]) not in prior]
        if not by_id[target_id]:
            del by_id[target_id]
    old = defaultdict(dict)
    broad = defaultdict(dict)
    for row in read_rows(args.v3):
        if row["target_id"] in by_id:
            old[row["target_id"]][row["method"]] = row
    for row in read_rows(args.allcopies):
        if row["target_id"] in by_id:
            broad[row["target_id"]][row["method"]] = row
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    with args.out_csv.open("a" if args.resume else "x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        if not args.resume:
            writer.writeheader()
        with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize_worker) as pool:
            futures = {pool.submit(score_target, target_id, rows, str(root),
                                   str(args.crystal_dir.resolve()), old[target_id],
                                   broad[target_id], args.contact_d,
                                   args.timeout_seconds): target_id
                       for target_id, rows in by_id.items()}
            for future in as_completed(futures):
                target_id = futures[future]
                try:
                    completed = future.result()
                except Exception as error:  # noqa: BLE001 - retain worker failure provenance
                    completed = [{"target_id": target_id, "method": row["method"],
                                  "status": "worker_failed",
                                  "error": f"{type(error).__name__}: {error}"}
                                 for row in by_id[target_id]]
                writer.writerows(completed)
                counts.update(str(row["status"]) for row in completed)
                handle.flush()
    print(f"Wrote {sum(counts.values())} static-QS rows: {dict(counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
