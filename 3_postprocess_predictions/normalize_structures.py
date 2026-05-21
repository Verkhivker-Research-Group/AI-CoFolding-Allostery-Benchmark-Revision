# 3_postprocess_predictions/normalize_structures.py
from __future__ import annotations

import argparse
from pathlib import Path

from ost import io


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--in-struct", required=True)
    p.add_argument("--out-cif", required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    ent = io.LoadEntity(str(args.in_struct), format="auto")
    v = ent.Select("ele != H")
    out = Path(args.out_cif)
    out.parent.mkdir(parents=True, exist_ok=True)
    io.SaveMMCIF(v, str(out))
    print(f"Saved -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
