# 4_score/score_lddt_pli.py
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from ost.mol.alg.ligand_scoring_lddtpli import LDDTPLIScorer

from ost_utils import filter_ligands, load_complex_mmcif


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--pred-cif", required=True)
    p.add_argument("--ref-cif", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--substructure-match", action="store_true")
    p.add_argument("--exclude-resnames", default="HOH,WAT,DOD")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    exclude = {x.strip().upper() for x in args.exclude_resnames.split(",") if x.strip()}

    mdl = load_complex_mmcif(args.pred_cif, extract_nonpoly=True)
    trg = load_complex_mmcif(args.ref_cif, extract_nonpoly=True)

    mdl_ligs = filter_ligands(mdl.ligands, exclude)
    trg_ligs = filter_ligands(trg.ligands, exclude)

    sc = LDDTPLIScorer(
        mdl.ent,
        trg.ent,
        [r.Select("ele != H") for r in mdl_ligs],
        [r.Select("ele != H") for r in trg_ligs],
        substructure_match=bool(args.substructure_match),
    )

    rows = []
    for (trg_i, mdl_i) in sc.assignment:
        rows.append(
            {
                "target_idx": trg_i,
                "model_idx": mdl_i,
                "target_resname": trg_ligs[trg_i].GetName(),
                "model_resname": mdl_ligs[mdl_i].GetName(),
                "LDDT_PLI": float(sc.score_matrix[trg_i, mdl_i]),
            }
        )

    if not rows:
        raise RuntimeError("No ligand assignment produced by scorer.")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote -> {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
