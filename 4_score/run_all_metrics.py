# 4_score/run_all_metrics.py
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from ost import io
from ost.mol.alg import qsscore, superpose
from ost.mol.alg.ligand_scoring_lddtpli import LDDTPLIScorer
from ost.mol.alg.ligand_scoring_scrmsd import SCRMSDScorer

from ost_utils import filter_ligands, load_complex_mmcif


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--pred-cif", required=True)
    p.add_argument("--ref-cif", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--exclude-resnames", default="HOH,WAT,DOD")
    p.add_argument("--substructure-match", action="store_true")
    p.add_argument("--pocket-json", default="")
    p.add_argument("--binding-site-json", default="")
    return p.parse_args()


def ca_sel(res: dict) -> str:
    return f"cname={res['chain']} and rnum={int(res['resnum'])} and aname=CA"


def ca_rmsd(pred_ent, ref_ent, residues: list[dict]) -> float:
    sel_expr = " or ".join([ca_sel(r) for r in residues])
    pv = pred_ent.Select(sel_expr)
    rv = ref_ent.Select(sel_expr)
    sp = superpose.SuperposeSVD(pv, rv)
    return float(sp.rmsd)


def main() -> int:
    args = parse_args()
    exclude = {x.strip().upper() for x in args.exclude_resnames.split(",") if x.strip()}

    mdl = load_complex_mmcif(args.pred_cif, extract_nonpoly=True)
    trg = load_complex_mmcif(args.ref_cif, extract_nonpoly=True)

    mdl_ligs = filter_ligands(mdl.ligands, exclude)
    trg_ligs = filter_ligands(trg.ligands, exclude)

    scr = SCRMSDScorer(
        mdl.ent,
        trg.ent,
        [r.Select("ele != H") for r in mdl_ligs],
        [r.Select("ele != H") for r in trg_ligs],
        substructure_match=bool(args.substructure_match),
    )
    trg_i, mdl_i = scr.assignment[0]
    bisyrmsd = float(scr.score_matrix[trg_i, mdl_i])

    pli = LDDTPLIScorer(
        mdl.ent,
        trg.ent,
        [r.Select("ele != H") for r in mdl_ligs],
        [r.Select("ele != H") for r in trg_ligs],
        substructure_match=bool(args.substructure_match),
    )
    trg_i2, mdl_i2 = pli.assignment[0]
    lddt_pli = float(pli.score_matrix[trg_i2, mdl_i2])

    mdl_q = qsscore.QSEntity(io.LoadEntity(str(args.pred_cif), format="auto"))
    ref_q = qsscore.QSEntity(io.LoadEntity(str(args.ref_cif), format="auto"))
    qsres = qsscore.QSScorer(ref_q, mdl_q).Score()
    qs_global = getattr(qsres, "qsματο_global", getattr(qsres, "qs_global", ""))

    row = {"BiSyRMSD": bisyrmsd, "LDDT_PLI": lddt_pli, "QS_global": qs_global}

    if args.pocket_json:
        pocket = json.loads(Path(args.pocket_json).read_text()).get("pocket_residues", [])
        row["pocket_CA_RMSD"] = ca_rmsd(mdl.ent, trg.ent, pocket) if pocket else ""

    if args.binding_site_json:
        bs = json.loads(Path(args.binding_site_json).read_text()).get("binding_site_residues", [])
        row["binding_site_CA_RMSD"] = ca_rmsd(mdl.ent, trg.ent, bs) if bs else ""

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        w.writeheader()
        w.writerow(row)

    print(f"Wrote -> {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
