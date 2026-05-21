# 2_run_models/validate_predictions.py
from __future__ import annotations

import argparse
import csv
from pathlib import Path

ALLOWED_EXTS = {".cif", ".mmcif", ".pdb"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--pred-root", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--min-bytes", type=int, default=5000)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    pred_root = Path(args.pred_root)

    rows = []
    for pdb_dir in sorted([p for p in pred_root.iterdir() if p.is_dir()]):
        models = sorted([f for f in pdb_dir.iterdir() if f.is_file() and f.suffix.lower() in ALLOWED_EXTS])
        status = "OK"
        msg = ""
        if models:
            tiny = [f.name for f in models if f.stat().st_size < args.min_bytes]
            if tiny:
                status = "WARN"
                msg = f"tiny files: {', '.join(tiny[:5])}" + ("..." if len(tiny) > 5 else "")
        else:
            status = "FAIL"
            msg = "no model files found"

        rows.append({"PDB": pdb_dir.name, "n_models": len(models), "status": status, "message": msg})

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["PDB", "n_models", "status", "message"])
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
