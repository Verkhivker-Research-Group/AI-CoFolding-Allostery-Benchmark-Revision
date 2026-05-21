# 2_run_models/collect_boltz2_outputs.py
from __future__ import annotations

import argparse
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

STRUCT_EXTS = {".cif", ".mmcif", ".pdb"}


@dataclass(frozen=True)
class Candidate:
    path: Path
    pdb_id: str
    score_hint: float | None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--boltz-root", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--max-per-target", type=int, default=5)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def infer_pdb_id(path: Path) -> str:
    tokens = re.split(r"[^A-Za-z0-9]+", (path.stem + "_" + path.parent.name).upper())
    for t in tokens:
        if len(t) == 4 and t.isalnum():
            return t
    s = path.stem.upper()
    return (s[:4] if len(s) >= 4 else s).ljust(4, "X")


def score_hint_from_name(path: Path) -> float | None:
    name = path.name.lower()

    m = re.search(r"(?:rank|sample)[_-]?(\d+)", name)
    if m:
        return float(m.group(1))

    m = re.search(r"(?:ptm|plddt|confidence)[_-]?(\d+(?:\.\d+)?)", name)
    if m:
        return -float(m.group(1))

    return None


def iter_structure_files(root: Path) -> Iterable[Path]:
    for p in root.rglob("*"):
        if p.is_file() and p.suffix.lower() in STRUCT_EXTS:
            yield p


def choose_top(cands: list[Candidate], k: int) -> list[Candidate]:
    def key(c: Candidate) -> tuple[int, float, str]:
        has = 0 if c.score_hint is not None else 1
        hint = c.score_hint if c.score_hint is not None else 1e18
        return (has, hint, str(c.path))

    return sorted(cands, key=key)[:k]


def main() -> int:
    args = parse_args()
    root = Path(args.boltz_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    by_pdb: dict[str, list[Candidate]] = {}
    for f in iter_structure_files(root):
        pdb_id = infer_pdb_id(f)
        by_pdb.setdefault(pdb_id, []).append(Candidate(f, pdb_id, score_hint_from_name(f)))

    wrote = 0
    for pdb_id, cands in sorted(by_pdb.items()):
        top = choose_top(cands, args.max_per_target)
        target_dir = out_dir / pdb_id
        target_dir.mkdir(parents=True, exist_ok=True)

        for i, cand in enumerate(top):
            dst = target_dir / f"model_{i:03d}{cand.path.suffix.lower()}"
            if dst.exists() and not args.overwrite:
                continue
            shutil.copy2(cand.path, dst)
            wrote += 1

    print(f"Done. Wrote {wrote} files into {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
