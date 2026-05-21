# 3_postprocess_predictions/extract_entities_from_cif.py
from __future__ import annotations

import argparse
from pathlib import Path

from ost import io
from ost.mol.alg.scoring_base import MMCIFPrep


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--in-cif", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--write-receptor", action="store_true")
    p.add_argument("--write-ligands", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    in_cif = Path(args.in_cif)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ent, ligs = MMCIFPrep(str(in_cif), extract_nonpoly=True)

    if args.write_receptor:
        rec = ent.Select("polymer")
        io.SaveMMCIF(rec, str(out_dir / f"{in_cif.stem}_receptor.cif"))

    if args.write_ligands:
        for i, lig in enumerate(ligs):
            lv = lig.Select("ele != H")
            io.SaveMMCIF(lv, str(out_dir / f"{in_cif.stem}_lig_{i:02d}_{lig.GetName()}.cif"))

    print(f"Done -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
