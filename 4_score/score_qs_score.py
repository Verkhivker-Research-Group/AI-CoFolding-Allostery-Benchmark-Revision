# 4_score/score_qs_score.py
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from ost import io
from ost.mol.alg import qsscore


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--pred-cif", required=True)
    p.add_argument("--ref-cif", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--contact-d", type=float, default=12.0)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    mdl = io.LoadEntity(str(args.pred_cif), format="auto")
    ref = io.LoadEntity(str(args.ref_cif), format="auto")

    mdl_q = qsscore.QSEntity(mdl, contact_d=float(args.contact_d))
    ref_q = qsscore.QSEntity(ref, contact_d=float(args.contact_d))

    res = qsscore.QSScorer(ref_q, mdl_q).Score()

    row = {
        "QS_global": getattr(res, "qs_global", ""),
        "QS_best": getattr(res, "qs_best", ""),
        "ICS": getattr(res, "ics", ""),
        "IPS": getattr(res, "ips", ""),
        "n_contacts_ref": getattr(res, "n_contacts_ref", ""),
        "n_contacts_mdl": getattr(res, "n_contacts_mdl", ""),
    }

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
