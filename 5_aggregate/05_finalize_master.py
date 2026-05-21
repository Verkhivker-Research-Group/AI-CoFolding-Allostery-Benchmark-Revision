import sys

# Allow running as a standalone script without installing as a package.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

"""Finalize a master results table.

Replaces the last notebook steps:
- rename columns to final schema
- optional row filtering
- drop rows with missing values
"""

from __future__ import annotations

import argparse

import pandas as pd

from scripts.utils.csv_ops import write_csv


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--dropna", action="store_true")
    ap.add_argument("--exclude_source", default=None)
    ap.add_argument("--exclude_dataset", default=None)
    ap.add_argument("--rename_to_final", action="store_true")
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)

    if args.rename_to_final:
        rename_map = {
            "RMSD": "ligand_rmsd",
            "model_rmsd": "rmsd",
            "model_tm_score": "tm_score",
            "model_plddt": "plddt",
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    if args.exclude_source and args.exclude_dataset:
        df = df[~((df["Source"].astype(str).str.lower() == args.exclude_source.lower()) &
                  (df["dataset"].astype(str).str.lower() == args.exclude_dataset.lower()))]

    if args.dropna:
        df = df.dropna()

    write_csv(df, args.out_csv)
    print(f"[ok] wrote {args.out_csv} (rows={len(df)})")


if __name__ == "__main__":
    main()
