import sys

# Allow running as a standalone script without installing as a package.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

"""Normalize a combined spreadsheet into the canonical master schema.

Replaces notebook cell 7.

Input columns expected:
- Source
- dataset (optional)
- PDB_ID
- RMSD

This script will:
- normalize PDB_ID (seed handling, underscore removal)
- normalize Source + dataset inferred from Source if missing
- output columns ordered as: Source, dataset, PDB_ID, RMSD
"""

from __future__ import annotations

import argparse

import pandas as pd

from scripts.utils.csv_ops import write_csv
from scripts.utils.ids import normalize_pdb_id, parse_source_and_dataset


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--pdb_col", default="PDB_ID")
    ap.add_argument("--source_col", default="Source")
    ap.add_argument("--dataset_col", default="dataset")
    ap.add_argument("--rmsd_col", default="RMSD")
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)

    # Ensure columns exist
    for c in [args.pdb_col, args.source_col, args.rmsd_col]:
        if c not in df.columns:
            raise SystemExit(f"Missing required column: {c}")

    # Normalize PDB_ID
    df[args.pdb_col] = df[args.pdb_col].apply(normalize_pdb_id)

    # Normalize Source/dataset
    if args.dataset_col not in df.columns:
        df[args.dataset_col] = ""

    parsed = df[args.source_col].apply(parse_source_and_dataset)
    df["Source"] = parsed.apply(lambda x: x.source)
    # Only overwrite dataset if empty
    inferred_dataset = parsed.apply(lambda x: x.dataset)
    df["dataset"] = df[args.dataset_col].where(df[args.dataset_col].astype(str).str.strip() != "", inferred_dataset)

    out = df[["Source", "dataset", args.pdb_col, args.rmsd_col]].copy()
    out = out.rename(columns={args.pdb_col: "PDB_ID", args.rmsd_col: "RMSD"})

    write_csv(out, args.out_csv)
    print(f"[ok] wrote {args.out_csv}")


if __name__ == "__main__":
    main()
