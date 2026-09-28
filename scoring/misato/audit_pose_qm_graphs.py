"""Compare every delivered pose's heavy-atom connectivity with QM.hdf5."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import networkx as nx
from rdkit import Chem, RDLogger


RDLogger.DisableLog("rdApp.*")
FIELDS = ["target_id", "method", "pose_path", "graph_status", "qm_heavy_atoms",
          "pose_heavy_atoms", "qm_bonds", "pose_bonds", "error"]


def pose_graph(path: Path) -> tuple[list[int], set[tuple[int, int]]]:
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
    mol = next((item for item in supplier if item is not None and item.GetNumConformers()), None)
    if mol is None:
        raise ValueError("unreadable_sdf")
    heavy = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    remap = {old: new for new, old in enumerate(heavy)}
    elements = [mol.GetAtomWithIdx(index).GetAtomicNum() for index in heavy]
    edges = {tuple(sorted((remap[bond.GetBeginAtomIdx()], remap[bond.GetEndAtomIdx()])))
             for bond in mol.GetBonds()
             if bond.GetBeginAtomIdx() in remap and bond.GetEndAtomIdx() in remap}
    return elements, edges


def graph_isomorphic(qm_atoms: list[int], qm_edges: set[tuple[int, int]],
                     pose_atoms: list[int], pose_edges: set[tuple[int, int]]) -> bool:
    if len(qm_atoms) != len(pose_atoms) or len(qm_edges) != len(pose_edges):
        return False
    if sorted(qm_atoms) != sorted(pose_atoms):
        return False
    if qm_atoms == pose_atoms and qm_edges == pose_edges:
        return True
    def build(atoms, edges):
        graph = nx.Graph()
        graph.add_nodes_from((index, {"element": number}) for index, number in enumerate(atoms))
        graph.add_edges_from(edges)
        return graph
    return nx.is_isomorphic(build(qm_atoms, qm_edges), build(pose_atoms, pose_edges),
                           node_match=lambda a, b: a["element"] == b["element"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qm-graphs", type=Path,
                        default=Path("misato_output/inventory/qm_heavy_graphs_v1.jsonl"))
    parser.add_argument("--manifest", type=Path,
                        default=Path("misato_output/inventory/target_manifest.csv"))
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.exists():
        parser.error("--out-dir exists; choose a new versioned directory")
    qm = {row["target_id"]: row for row in
          (json.loads(line) for line in args.qm_graphs.read_text(encoding="utf-8").splitlines())}
    with args.manifest.open(newline="", encoding="utf-8-sig") as handle:
        manifest = list(csv.DictReader(handle))
    args.out_dir.mkdir(parents=True)
    statuses = {method: Counter() for method in ("diffdock", "equibind")}
    with (args.out_dir / "pose_qm_graph_audit.csv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for target in manifest:
            target_id = target["target_id"].upper()
            source = qm.get(target_id)
            for method in statuses:
                if target[f"{method}_status"] != "candidate":
                    continue
                path = Path(target[f"{method}_primary_path"])
                row = {"target_id": target_id, "method": method, "pose_path": str(path),
                       "graph_status": "", "qm_heavy_atoms": "", "pose_heavy_atoms": "",
                       "qm_bonds": "", "pose_bonds": "", "error": ""}
                try:
                    if not source or source["status"] != "ok":
                        raise ValueError("qm_graph_unavailable")
                    qm_atoms = source["atomic_numbers"]
                    qm_edges = {tuple(bond[:2]) for bond in source["bonds"]}
                    elements, edges = pose_graph(path)
                    row["qm_heavy_atoms"] = len(qm_atoms)
                    row["pose_heavy_atoms"] = len(elements)
                    row["qm_bonds"] = len(qm_edges)
                    row["pose_bonds"] = len(edges)
                    if qm_atoms == elements and qm_edges == edges:
                        row["graph_status"] = "same_indexed_graph"
                    elif graph_isomorphic(qm_atoms, qm_edges, elements, edges):
                        row["graph_status"] = "isomorphic_reordered_graph"
                    else:
                        row["graph_status"] = "graph_mismatch"
                except Exception as error:
                    row["graph_status"] = "audit_failed"
                    row["error"] = f"{type(error).__name__}: {error}"
                writer.writerow(row)
                statuses[method][row["graph_status"]] += 1
    summary = {
        "input_qm_graphs": str(args.qm_graphs),
        "by_method": {method: dict(sorted(counts.items())) for method, counts in statuses.items()},
        "interpretation": "Heavy-atom topology only; atom mapping, bond order, protonation, stereochemistry, and initial coordinates require separate checks.",
    }
    (args.out_dir / "pose_qm_graph_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
