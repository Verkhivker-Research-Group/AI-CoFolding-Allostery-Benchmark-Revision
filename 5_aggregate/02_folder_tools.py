"""Folder-level utilities.

Replaces notebook cells:
- filter each CSV to leftmost col + RMSD
- combine CSVs with a Source column from filename
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import sys

# Allow running as a standalone script without installing as a package.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


import pandas as pd

from scripts.utils.csv_ops import write_csv


def filter_to_leftmost_and_rmsd(folder: str, pattern: str = "*.csv", rmsd_col: str = "RMSD") -> None:
    for file_path in glob.glob(os.path.join(folder, pattern)):
        df = pd.read_csv(file_path)
        leftmost = df.columns[0]
        if rmsd_col not in df.columns:
            print(f"[skip] {file_path} has no '{rmsd_col}' column")
            continue
        out = df[[leftmost, rmsd_col]].copy()
        out_path = str(file_path).replace(".csv", "_filtered.csv")
        write_csv(out, out_path)
        print(f"[ok] wrote {out_path}")


def combine_csvs_with_source(folder: str, out_csv: str, pattern: str = "*.csv") -> None:
    dfs = []
    for file_path in glob.glob(os.path.join(folder, pattern)):
        df = pd.read_csv(file_path)
        name = os.path.splitext(os.path.basename(file_path))[0]
        df.insert(0, "Source", name)
        dfs.append(df)
    if not dfs:
        raise SystemExit(f"No CSVs matched {pattern} in {folder}")
    combined = pd.concat(dfs, ignore_index=True)
    write_csv(combined, out_csv)
    print(f"[ok] wrote {out_csv}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("filter")
    p1.add_argument("--folder", required=True)
    p1.add_argument("--pattern", default="*.csv")
    p1.add_argument("--rmsd_col", default="RMSD")

    p2 = sub.add_parser("combine")
    p2.add_argument("--folder", required=True)
    p2.add_argument("--pattern", default="*.csv")
    p2.add_argument("--out_csv", required=True)

    args = ap.parse_args()

    Path(args.folder).mkdir(parents=True, exist_ok=True)

    if args.cmd == "filter":
        filter_to_leftmost_and_rmsd(args.folder, args.pattern, args.rmsd_col)
    else:
        combine_csvs_with_source(args.folder, args.out_csv, args.pattern)


if __name__ == "__main__":
    main()
