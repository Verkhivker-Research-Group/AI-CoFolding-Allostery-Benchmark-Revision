"""CLI entry points.

Three subcommands:

* ``plb-bench normalize``  — raw producer outputs → ``data/<producer>/<pdb_id>/model_NNN.cif``
* ``plb-bench score``      — score an already-normalized tree with OST
* ``plb-bench run``        — legacy: raw tree → scored table in one shot
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .pipeline import PipelineConfig, run as legacy_run
from .normalize import DISCOVERERS, normalize_producer
from .score_tree import ScoreConfig, run as score_run
from .registry import PARSERS


def _add_normalize(sub):
    n = sub.add_parser("normalize",
                        help="Turn a producer's raw outputs into data/<producer>/<pdb_id>/model_NNN.cif")
    n.add_argument("--producer", required=True, choices=sorted(DISCOVERERS),
                    help="Which producer this raw tree is for")
    n.add_argument("--raw", required=True, type=Path,
                    help="Producer's raw root directory (contains <pdb_id>/ subdirs)")
    n.add_argument("--out", required=True, type=Path,
                    help="Normalized data root (output goes to <out>/<producer>/<pdb_id>/...)")
    n.add_argument("--no-overwrite", action="store_true",
                    help="Leave stale model_*.cif from previous runs in place")
    n.add_argument("-v", "--verbose", action="store_true")


def _add_score(sub):
    s = sub.add_parser("score", help="Score a normalized data tree with OST")
    s.add_argument("--data", required=True, type=Path,
                    help="Normalized data root (expects <producer>/<pdb_id>/model_NNN.cif layout)")
    s.add_argument("--refs", required=True, type=Path)
    s.add_argument("--out", required=True, type=Path)
    s.add_argument("--producers", nargs="+", help="Restrict to a subset")
    s.add_argument("--pdb-ids", nargs="+", help="Restrict to a subset of PDB IDs")
    s.add_argument("--workers", type=int, default=None)
    s.add_argument("--ray", action="store_true")
    s.add_argument("--no-download", action="store_true")
    s.add_argument("--format", choices=("parquet", "csv"), default="parquet")
    s.add_argument("-v", "--verbose", action="store_true")


def _add_run(sub):
    r = sub.add_parser("run", help="Legacy: raw tree (<pdb_id>/<producer>/) → scored table")
    r.add_argument("--input", required=True, type=Path)
    r.add_argument("--refs", required=True, type=Path)
    r.add_argument("--out", required=True, type=Path)
    r.add_argument("--producers", nargs="+", choices=sorted(PARSERS))
    r.add_argument("--workers", type=int, default=None)
    r.add_argument("--ray", action="store_true")
    r.add_argument("--no-download", action="store_true")
    r.add_argument("--no-canonical", action="store_true")
    r.add_argument("--format", choices=("parquet", "csv"), default="parquet")
    r.add_argument("-v", "--verbose", action="store_true")


def main() -> int:
    p = argparse.ArgumentParser("plb-bench", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    _add_normalize(sub)
    _add_score(sub)
    _add_run(sub)
    args = p.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if getattr(args, "verbose", False) else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.cmd == "normalize":
        entries = normalize_producer(
            producer=args.producer,
            raw_root=args.raw,
            output_root=args.out,
            overwrite=not args.no_overwrite,
        )
        n_ok = sum(e.n_models for e in entries)
        n_fail = sum(e.failed for e in entries)
        print(f"Normalized {len(entries)} (producer, pdb_id) groups: "
              f"{n_ok} models written, {n_fail} failed")
        print(f"Output root: {args.out / args.producer}")
        return 0

    if args.cmd == "score":
        cfg = ScoreConfig(
            normalized_root=args.data,
            refs_dir=args.refs,
            output_dir=args.out,
            producers=args.producers,
            pdb_ids=args.pdb_ids,
            n_workers=args.workers,
            use_ray=args.ray,
            download_missing_refs=not args.no_download,
            output_format=args.format,
        )
        out = score_run(cfg)
        print(f"Done. Output: {out}")
        return 0

    if args.cmd == "run":
        cfg = PipelineConfig(
            input_root=args.input,
            refs_dir=args.refs,
            output_dir=args.out,
            producers=args.producers,
            n_workers=args.workers,
            use_ray=args.ray,
            download_missing_refs=not args.no_download,
            write_canonical_models=not args.no_canonical,
            output_format=args.format,
        )
        out = legacy_run(cfg)
        print(f"Done. Output: {out}")
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
