# 4_score/ost_utils.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import ost
from ost.mol.alg.scoring_base import MMCIFPrep


@dataclass(frozen=True)
class LoadedComplex:
    ent: ost.mol.EntityHandle
    ligands: list


def load_complex_mmcif(path: str | Path, extract_nonpoly: bool = True) -> LoadedComplex:
    ent, ligs = MMCIFPrep(str(path), extract_nonpoly=extract_nonpoly)
    return LoadedComplex(ent=ent, ligands=list(ligs))


def ligand_residue_to_id(res) -> dict[str, Any]:
    num = res.GetNumber()
    return {
        "chain": res.GetChain().GetName(),
        "resnum": int(num.GetNum()),
        "ins": str(num.GetInsCode() or ""),
        "resname": res.GetName(),
        "n_atoms": res.Select("ele != H").GetAtomCount(),
    }


def filter_ligands(ligs: list, exclude_resnames: set[str]) -> list:
    out = []
    for r in ligs:
        if r.GetName().upper() in exclude_resnames:
            continue
        if r.Select("ele != H").GetAtomCount() == 0:
            continue
        out.append(r)
    return out


def pick_largest_ligand(ligs: list):
    return max(ligs, key=lambda r: r.Select("ele != H").GetAtomCount())
