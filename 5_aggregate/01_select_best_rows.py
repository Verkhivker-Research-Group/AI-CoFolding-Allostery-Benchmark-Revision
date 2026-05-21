"""Select best (lowest RMSD or highest score) rows from raw model outputs.

Replaces notebook cells like:
- AF3 allo: group by Complex_Seed prefix and keep min RMSD
- Protenix allo: compute max across columns, keep per-protein result
"""

from __future__ import annotations

import argparse
from pathlib import Path

import sys

# Allow running as a standalone script without installing as a package.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


import pandas as pd

from scripts.utils.csv_ops import best_row_per_group, write_csv
from scripts.utils.ids import extract_prefix


def af3_best_by_complex_seed(in_csv: str, out_csv: str, seed_col: str = "Complex_Seed", rmsd_col: str = "RMSD") -> None:
    df = pd.read_csv(in_csv)
    df["group_code"] = df[seed_col].astype(str).apply(lambda x: extract_prefix(x, 5)).str.upper()
    best = best_row_per_group(df, "group_code", rmsd_col, keep="min")
    write_csv(best.drop(columns=["group_code"], errors="ignore"), out_csv)


def protenix_rowwise_max_to_rmsd(in_csv: str, out_csv: str, id_col: str = "Protein ID") -> None:
    df = pd.read_csv(in_csv)
    rmsd_cols = [c for c in df.columns if c != id_col]
    df["RMSD"] = df[rmsd_cols].apply(pd.to_numeric, errors="coerce").max(axis=1)
    out = df[[id_col, "RMSD"]].copy()
    write_csv(out, out_csv)


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)

    p1 = sub.add_parser("af3_best")
    p1.add_argument("--in_csv", required=True)
    p1.add_argument("--out_csv", required=True)
    p1.add_argument("--seed_col", default="Complex_Seed")
    p1.add_argument("--rmsd_col", default="RMSD")

    p2 = sub.add_parser("protenix_max")
    p2.add_argument("--in_csv", required=True)
    p2.add_argument("--out_csv", required=True)
    p2.add_argument("--id_col", default="Protein ID")

    args = ap.parse_args()

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)

    if args.mode == "af3_best":
        af3_best_by_complex_seed(args.in_csv, args.out_csv, args.seed_col, args.rmsd_col)
    elif args.mode == "protenix_max":
        protenix_rowwise_max_to_rmsd(args.in_csv, args.out_csv, args.id_col)


if __name__ == "__main__":
    main()
