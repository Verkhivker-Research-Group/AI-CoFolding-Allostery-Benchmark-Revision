# 3_postprocess_predictions/build_pocket_residue_set.py
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ost.mol.alg.scoring_base import MMCIFPrep


@dataclass(frozen=True)
class ResID:
    chain: str
    resnum: int
    ins: str
    resname: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--ref-cif", required=True)
    p.add_argument("--ligand-json", required=True)
    p.add_argument("--radius", type=float, default=8.0)
    p.add_argument("--out-json", required=True)
    return p.parse_args()


def load_ligand_selections(lig_json: Path) -> list[str]:
    data = json.loads(lig_json.read_text())
    return [x["ost_residue_sel"] for x in data["ligands"]]


def unique_residues(view) -> list[ResID]:
    out: dict[tuple[str, int, str], ResID] = {}
    for res in view.GetResidues():
        num = res.GetNumber()
        key = (res.GetChain().GetName(), int(num.GetNum()), str(num.GetInsCode() or ""))
        out[key] = ResID(chain=key[0], resnum=key[1], ins=key[2], resname=res.GetName())
    return list(out.values())


def main() -> int:
    args = parse_args()
    ref_cif = Path(args.ref_cif)
    lig_json = Path(args.ligand_json)

    ent, _ = MMCIFPrep(str(ref_cif), extract_nonpoly=False)
    lig_sels = load_ligand_selections(lig_json)

    prot = ent.Select("protein and ele != H")
    near_atoms = prot.Select(f"within {args.radius} of ({' or '.join(lig_sels)})")
    pocket_res = unique_residues(near_atoms)

    out: dict[str, Any] = {
        "ref_cif": str(ref_cif),
        "ligand_json": str(lig_json),
        "radius": float(args.radius),
        "pocket_residues": [{"chain": r.chain, "resnum": r.resnum, "ins": r.ins, "resname": r.resname} for r in pocket_res],
    }

    out_path = Path(args.out_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    print(f"Wrote -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
