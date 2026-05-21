from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Allow `python scripts/...` to import `scripts.utils.*`
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.utils.rcsb import download_structure  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--excel", required=True, help="Path to Excel file (e.g., data.xlsx)")
    p.add_argument("--column", default="entryName", help="Column containing entry names (default: entryName)")
    p.add_argument("--fmt", choices=["pdb", "cif"], required=True, help="Download format")
    p.add_argument("--out_dir", required=True, help="Output directory")
    p.add_argument("--overwrite", action="store_true", help="Re-download even if file exists")
    p.add_argument("--sleep_s", type=float, default=0.0, help="Sleep between requests")
    p.add_argument("--limit", type=int, default=None, help="Optional limit number of IDs")
    return p.parse_args()


def extract_pdb_id(entry: str) -> str:
    entry = str(entry)
    return entry.split("_")[0].strip()


def main() -> None:
    args = parse_args()
    df = pd.read_excel(args.excel)
    if args.column not in df.columns:
        raise ValueError(f"Column '{args.column}' not found. Columns: {list(df.columns)}")

    pdb_ids = (
        df[args.column]
        .dropna()
        .map(extract_pdb_id)
        .drop_duplicates()
        .tolist()
    )

    if args.limit is not None:
        pdb_ids = pdb_ids[: args.limit]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ok = 0
    bad = 0
    for pid in pdb_ids:
        res = download_structure(pid, args.fmt, out_dir, overwrite=args.overwrite, sleep_s=args.sleep_s)
        if res.ok:
            ok += 1
            print(f"OK  {res.pdb_id}.{res.fmt} -> {res.path}")
        else:
            bad += 1
            print(f"ERR {res.pdb_id}.{res.fmt} ({res.status_code}) {res.error}")

    print(f"Done. Success={ok} Failed={bad} Total={len(pdb_ids)}")


if __name__ == "__main__":
    main()
