"""Phase 1: Input Harmonization & Confidence Scoring.

Each generative model dumps its outputs differently. We crawl a directory
organized as::

    root/
      <pdb_id>/
        af3/           <- AF3 output dir (JSON + CIF)
        proteinx/      <- ProteinX output dir
        chai/          <- Chai-1 output dir (NPZ + CIF)
        dynamicbind/   <- DynamicBind output (PDB + SDF)

Any missing producer is silently skipped. Adding a new producer = one function
decorated with ``@register_parser``; no other code changes.
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Callable, Iterator

import numpy as np

from .schema import ModelRecord

log = logging.getLogger(__name__)

# Registry: producer_name -> parser callable(pdb_dir: Path, pdb_id: str) -> list[ModelRecord]
PARSERS: dict[str, Callable[[Path, str], list[ModelRecord]]] = {}


def register_parser(name: str):
    """Decorator to register a new producer parser."""
    def wrap(fn):
        PARSERS[name] = fn
        return fn
    return wrap


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _safe_json(path: Path) -> dict | None:
    try:
        with path.open() as fh:
            return json.load(fh)
    except Exception as e:
        log.warning("bad json %s: %s", path, e)
        return None


def _safe_npz(path: Path) -> dict | None:
    try:
        with np.load(path, allow_pickle=True) as arch:
            return {k: arch[k] for k in arch.files}
    except Exception as e:
        log.warning("bad npz %s: %s", path, e)
        return None


def _mean(x) -> float | None:
    try:
        arr = np.asarray(x, dtype=float)
        if arr.size == 0:
            return None
        return float(np.nanmean(arr))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# AlphaFold 3
# ---------------------------------------------------------------------------
# AF3 writes <job>/seed-*_sample-*/model.cif alongside summary_confidences.json
# and confidences.json. The canonical per-model score is `ranking_score`
# (combined pTM/ipTM/clash heuristic).
@register_parser("af3")
def parse_af3(pdb_dir: Path, pdb_id: str) -> list[ModelRecord]:
    root = pdb_dir / "af3"
    if not root.is_dir():
        return []
    out: list[ModelRecord] = []
    # AF3 sometimes nests one level (seed_*/sample_*), sometimes flat
    for cif in sorted(root.rglob("*.cif")):
        summary = cif.with_name("summary_confidences.json")
        if not summary.exists():
            # fall back to a sibling by stem
            cands = list(cif.parent.glob("*summary_confidences*.json"))
            summary = cands[0] if cands else None
        meta = _safe_json(summary) if summary else {}
        meta = meta or {}
        raw = {
            k: float(meta[k])
            for k in ("ranking_score", "iptm", "ptm", "fraction_disordered",
                     "has_clash", "num_recycles")
            if k in meta and isinstance(meta[k], (int, float))
        }
        # AF3 ranking_score is the best proxy for "pick this pose"
        confidence = raw.get("ranking_score", raw.get("iptm", raw.get("ptm")))
        out.append(ModelRecord(
            pdb_id=pdb_id, producer="af3", structure_path=cif,
            raw_scores=raw, confidence=confidence,
        ))
    return out


# ---------------------------------------------------------------------------
# ProteinX  (Bytedance / ByteDance Protenix — "proteinx" covers both spellings)
# ---------------------------------------------------------------------------
# Protenix produces predictions/<job>/seed_*/predictions/ with *.cif and a
# confidence json (confidence_sample_*.json). The headline metric is
# "complex_plddt" / "ranking_confidence".
@register_parser("proteinx")
def parse_proteinx(pdb_dir: Path, pdb_id: str) -> list[ModelRecord]:
    for alias in ("proteinx", "protenix"):
        root = pdb_dir / alias
        if root.is_dir():
            break
    else:
        return []
    out: list[ModelRecord] = []
    for cif in sorted(root.rglob("*.cif")):
        # Real Protenix file naming:
        #   <name>_<seed>_sample_N.cif
        #   <name>_<seed>_summary_confidence_sample_N.json
        # Pair by sample index, preferring the summary_confidence file.
        m = re.search(r"sample[_-]?(\d+)", cif.stem)
        sample_idx = m.group(1) if m else None
        conf_json = None
        if sample_idx is not None:
            # Prefer the summary_confidence file; fall back to any matching sample json
            preferred = list(cif.parent.glob(
                f"*summary_confidence*sample*{sample_idx}*.json"))
            if preferred:
                conf_json = preferred[0]
            else:
                cands = list(cif.parent.glob(f"*sample*{sample_idx}*.json"))
                if cands:
                    conf_json = cands[0]
        if conf_json is None:
            cands = list(cif.parent.glob("*summary_confidence*.json")) \
                    or list(cif.parent.glob("*confidence*.json"))
            conf_json = cands[0] if cands else None
        meta = _safe_json(conf_json) if conf_json else {}
        meta = meta or {}
        # Real Protenix summary_confidence fields (per bytedance/Protenix docs):
        #   ranking_score, plddt, gpde, ptm, iptm, has_clash, disorder,
        #   chain_ptm[], chain_iptm[], chain_plddt[], chain_pair_iptm[][], ...
        # The headline ranker is `ranking_score` (same convention as AF3).
        raw: dict[str, float] = {}
        for k in ("ranking_score", "plddt", "gpde", "ptm", "iptm",
                   "has_clash", "disorder",
                   # legacy / pipeline-wrapper names that some forks emit
                   "ranking_confidence", "complex_plddt"):
            v = meta.get(k)
            if isinstance(v, bool):
                raw[k] = float(v)
            elif isinstance(v, (int, float)):
                raw[k] = float(v)
            elif isinstance(v, list):
                mv = _mean(v)
                if mv is not None:
                    raw[k] = mv
        confidence = (raw.get("ranking_score")
                      or raw.get("ranking_confidence")
                      or raw.get("iptm")
                      or raw.get("ptm")
                      or raw.get("plddt"))
        out.append(ModelRecord(
            pdb_id=pdb_id, producer="proteinx", structure_path=cif,
            raw_scores=raw, confidence=confidence,
        ))
    return out


# ---------------------------------------------------------------------------
# Chai-1
# ---------------------------------------------------------------------------
# Chai emits pred.model_idx_<k>.cif + scores.model_idx_<k>.npz with keys
# `aggregate_score`, `ptm`, `iptm`, `per_chain_ptm`, `per_chain_pair_iptm`, etc.
@register_parser("chai")
def parse_chai(pdb_dir: Path, pdb_id: str) -> list[ModelRecord]:
    root = pdb_dir / "chai"
    if not root.is_dir():
        return []
    out: list[ModelRecord] = []
    for cif in sorted(root.rglob("pred.model_idx_*.cif")):
        m = re.search(r"model_idx_(\d+)", cif.stem)
        if not m:
            continue
        k = m.group(1)
        npz = cif.parent / f"scores.model_idx_{k}.npz"
        data = _safe_npz(npz) if npz.exists() else {}
        data = data or {}
        raw = {}
        for key in ("aggregate_score", "ptm", "iptm",
                    "has_inter_chain_clashes", "chain_chain_clashes"):
            if key in data:
                raw[key] = _mean(data[key]) if np.ndim(data[key]) > 0 else float(data[key])
        confidence = raw.get("aggregate_score", raw.get("iptm", raw.get("ptm")))
        out.append(ModelRecord(
            pdb_id=pdb_id, producer="chai", structure_path=cif,
            raw_scores=raw, confidence=confidence,
        ))
    return out


# ---------------------------------------------------------------------------
# DynamicBind
# ---------------------------------------------------------------------------
# DynamicBind docks a ligand into a flexible receptor. For each sample it writes
# rank<k>_receptor_*.pdb and rank<k>_ligand_*.sdf plus a complex_confidence.csv
# / affinity_prediction.csv. The rank number IS the ordering (rank1 = best).
@register_parser("dynamicbind")
def parse_dynamicbind(pdb_dir: Path, pdb_id: str) -> list[ModelRecord]:
    root = pdb_dir / "dynamicbind"
    if not root.is_dir():
        return []

    # Confidence csv (one row per rank)
    conf_map: dict[int, dict[str, float]] = {}
    for csv_name in ("complex_confidence.csv", "affinity_prediction.csv"):
        csv_path = root / csv_name
        if not csv_path.exists():
            continue
        try:
            import csv as _csv
            with csv_path.open() as fh:
                rdr = _csv.DictReader(fh)
                for row in rdr:
                    rank_str = row.get("rank") or row.get("model_rank") or row.get("index")
                    if rank_str is None:
                        continue
                    try:
                        rank = int(rank_str)
                    except ValueError:
                        continue
                    bucket = conf_map.setdefault(rank, {})
                    for k, v in row.items():
                        if k in ("rank", "model_rank", "index"):
                            continue
                        try:
                            bucket[k] = float(v)
                        except (TypeError, ValueError):
                            pass
        except Exception as e:
            log.warning("dynamicbind csv parse failed %s: %s", csv_path, e)

    out: list[ModelRecord] = []
    for pdb in sorted(root.rglob("rank*_receptor*.pdb")):
        m = re.search(r"rank(\d+)", pdb.name)
        if not m:
            continue
        rank = int(m.group(1))
        sdf_cands = list(pdb.parent.glob(f"rank{rank}_ligand*.sdf"))
        sdf = sdf_cands[0] if sdf_cands else None
        raw = conf_map.get(rank, {})
        # DynamicBind's "confidence" increases with rank-goodness; if we only have
        # the rank number, invert it so higher = better to stay consistent.
        confidence = raw.get("confidence", raw.get("affinity", -float(rank)))
        out.append(ModelRecord(
            pdb_id=pdb_id, producer="dynamicbind",
            structure_path=pdb, ligand_path=sdf,
            raw_scores=raw, confidence=confidence,
        ))
    return out


# ---------------------------------------------------------------------------
# Top-level crawling
# ---------------------------------------------------------------------------
def iter_pdb_dirs(root: Path) -> Iterator[tuple[str, Path]]:
    """Yield (pdb_id, pdb_dir) for every immediate subdirectory."""
    for p in sorted(root.iterdir()):
        if p.is_dir() and not p.name.startswith("."):
            yield p.name.lower(), p


def parse_all(pdb_dir: Path, pdb_id: str,
              producers: list[str] | None = None) -> list[ModelRecord]:
    """Run every registered parser (or a subset) on one pdb_id directory."""
    producers = producers or list(PARSERS.keys())
    records: list[ModelRecord] = []
    for name in producers:
        fn = PARSERS.get(name)
        if fn is None:
            continue
        try:
            records.extend(fn(pdb_dir, pdb_id))
        except Exception as e:
            log.exception("parser %s failed for %s: %s", name, pdb_id, e)
    return records


def rank_and_assign(records: list[ModelRecord]) -> list[ModelRecord]:
    """Within each producer, sort by confidence desc and assign model_idx 0..N-1."""
    by_producer: dict[str, list[ModelRecord]] = {}
    for r in records:
        by_producer.setdefault(r.producer, []).append(r)
    for prod, group in by_producer.items():
        group.sort(
            key=lambda r: (r.confidence if r.confidence is not None else -float("inf")),
            reverse=True,
        )
        for i, r in enumerate(group):
            r.model_idx = i
    return records
