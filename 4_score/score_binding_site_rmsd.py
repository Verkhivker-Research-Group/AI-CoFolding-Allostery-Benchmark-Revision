# 4_score/score_binding_site_rmsd.py
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from ost import io
from ost.mol.alg import superpose


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--pred-cif", required=True)
    p.add_argument("--ref-cif", required=True)
    p.add_argument("--binding-site-json", required=True)
    p.add_argument("--out-csv", required=True)
    return p.parse_args()


def res_sel(res: dict) -> str:
    return f"cname={res['chain']} and rnum={int(res['resnum'])} and aname=CA"


def main() -> int:
    args = parse_args()
    pred = io.LoadEntity(str(args.pred_cif), format="auto")
    ref = io.LoadEntity(str(args.ref_cif), format="auto")

    bs = json.loads(Path(args.bindi
