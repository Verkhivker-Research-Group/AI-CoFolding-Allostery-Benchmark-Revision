from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.utils.rcsb import download_structure


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument('--excel', required=True, help='Path to Excel file (e.g., data.xlsx)')
    p.add_argument('--column', default='entryName', help='Column containing PDB identifiers')
    p.add_argument('--out-dir', default='data/structures/cif', help='Output directory for downloaded .cif files')
    p.add_argument('--overwrite', action='store_true', help='Re-download even if file exists')
    p.add_argument('--sleep', type=float, default=0.0, help='Seconds to sleep between requests')
    return p.parse_args()


def extract_pdb_ids(series: pd.Series) -> list[str]:
    return (
        series.astype(str)
        .str.split('_')
        .str[0]
        .str.strip()
        .str.upper()
        .dropna()
        .drop_duplicates()
        .tolist()
    )


def main() -> int:
    args = parse_args()
    df = pd.read_excel(args.excel)
    if args.column not in df.columns:
        raise KeyError(f"Column '{args.column}' not found in {args.excel}. Available: {list(df.columns)}")

    pdb_ids = extract_pdb_ids(df[args.column])
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ok = 0
    bad = 0
    for pdb_id in pdb_ids:
        res = download_structure(pdb_id, 'cif', out_dir, overwrite=args.overwrite, sleep_s=args.sleep)
        if res.ok:
            ok += 1
        else:
            bad += 1
            print(f"[WARN] {pdb_id}: {res.error} ({res.status_code})")

    print(f"Done. Downloaded OK: {ok}  Failed: {bad}  Out: {out_dir}")
    return 0 if bad == 0 else 2


if __name__ == '__main__':
    raise SystemExit(main())
