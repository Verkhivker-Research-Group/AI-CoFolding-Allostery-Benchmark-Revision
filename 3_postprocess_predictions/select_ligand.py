# 3_postprocess_predictions/select_ligand.py
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ost.mol.alg.scoring_base import MMCIFPrep


@dataclass(frozen=True)
class LigandPick:
    chain: str
    resnum: int
    ins_code: str
    resname: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--ref-cif", required=True)
    p.add_argument("--out-json", required=True)
    p.add_argument("--allow-multiple", action="store_true")
    p.add_argument("--exclude-resnames", default="HOH,WAT,DOD")
    return p.parse_args()


def residue_to_pick(res) -> LigandPick:
    num = res.GetNumber()
    return LigandPick(
        chain=res.GetChain().GetName(),
        resnum=int(num.GetNum()),
        ins_code=str(num.GetInsCode() or ""),
        resname=str(res.GetName()),
    )


def main() -> int:
    args = parse_args()
    exclude = {x.strip().upper() for x in args.exclude_resnames.split(",") if x.strip()}

    ent, ligs = MMCIFPrep(str(args.ref_cif), extract_nonpoly=True)
    picks = [residue_to_pick(r) for r in ligs if r.GetName().upper() not in exclude]

    if not picks:
        raise RuntimeError("No ligands found after exclusions.")

    if not args.allow_multiple:
        def atom_count(p: LigandPick) -> int:
            sel = ent.Select(f"cname={p.chain} and rnum={p.resnum}")
            return sel.GetAtomCount()

        picks = [max(picks, key=atom_count)]

    out: dict[str, Any] = {
        "ref_cif": str(Path(args.ref_cif)),
        "ligands": [
            {
                "chain": p.chain,
                "resnum": p.resnum,
                "ins_code": p.ins_code,
                "resname": p.resname,
                "ost_residue_sel": f"cname={p.chain} and rnum={p.resnum}",
            }
            for p in picks
        ],
    }

    out_path = Path(args.out_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    print(f"Wrote -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
