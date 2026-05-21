# 4_score/score_pocket_rmsd.py
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from ost import io
from ost.mol.alg import superpose

POCKET_RADIUS_A = 8.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--pred-cif", required=True)
    p.add_argument("--ref-cif", required=True)
    p.add_argument("--pocket-json", required=True)
    p.add_argument("--out-csv", required=True)
    return p.parse_args()


def res_sel(res: dict) -> str:
    return f"cname={res['chain']} and rnum={int(res['resnum'])} and aname=CA"


def main() -> int:
    args = parse_args()
    pred = io.LoadEntity(str(args.pred_cif), format="auto")
    ref = io.LoadEntity(str(args.ref_cif), format="auto")

    pocket = json.loads(Path(args.pocket_json).read_text()).get("pocket_residues", [])
    if not pocket:
        raise RuntimeError("Pocket residue set is empty.")

    sel_expr = " or ".join([res_sel(r) for r in pocket])
    pred_view = pred.Select(sel_expr)
    ref_view = ref.Select(sel_expr)

    sp = superpose.SuperposeSVD(pred_view, ref_view)
    rmsd = float(sp.rmsd)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pocket_CA_RMSD", "n_atoms", "pocket_radius_A"])
        w.writeheader()
        w.writerow({"pocket_CA_RMS
