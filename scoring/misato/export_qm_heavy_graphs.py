"""Export only delivered QM ligand heavy-atom graphs for local identity audit."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import h5py

from audit_qm_coverage import PUBLISHED_QM_MD5, md5_file


def graph_from_group(target_id: str, group) -> dict:
    atoms = [int(item.decode("ascii")) for item in group["atom_properties"]["atom_names"][:]]
    heavy = [index for index, number in enumerate(atoms) if number > 1]
    old_to_new = {old: new for new, old in enumerate(heavy)}
    edges = {}
    for row in group["atom_properties"]["bonds"][:]:
        left, right = int(row[0]), int(row[1])
        if left in old_to_new and right in old_to_new and left != right:
            edge = tuple(sorted((old_to_new[left], old_to_new[right])))
            edges[edge] = float(row[2])
    return {"target_id": target_id, "atomic_numbers": [atoms[index] for index in heavy],
            "bonds": [[left, right, order] for (left, right), order in sorted(edges.items())]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qm", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=Path("misato_output/inventory/target_manifest.csv"))
    parser.add_argument("--out-jsonl", type=Path, required=True)
    args = parser.parse_args()
    if args.out_jsonl.exists():
        parser.error("--out-jsonl exists; choose a new versioned filename")
    if md5_file(args.qm) != PUBLISHED_QM_MD5:
        parser.error("QM.hdf5 does not match Zenodo MD5")
    with args.manifest.open(newline="", encoding="utf-8-sig") as handle:
        ids = sorted({row["target_id"].upper() for row in csv.DictReader(handle)})
    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    statuses = Counter()
    with h5py.File(args.qm, "r") as source, args.out_jsonl.open("x", encoding="utf-8") as output:
        for target_id in ids:
            try:
                row = graph_from_group(target_id, source[target_id])
                row["status"] = "ok"
            except (KeyError, ValueError, TypeError, UnicodeDecodeError) as error:
                row = {"target_id": target_id, "status": "failed",
                       "error": f"{type(error).__name__}: {error}"}
            statuses[row["status"]] += 1
            output.write(json.dumps(row, separators=(",", ":")) + "\n")
    print(f"QM graphs: {len(ids)}; statuses: {dict(statuses)}; output: {args.out_jsonl}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
