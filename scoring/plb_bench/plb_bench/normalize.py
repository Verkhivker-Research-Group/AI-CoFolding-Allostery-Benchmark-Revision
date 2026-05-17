"""Normalization: raw producer outputs → unified ``data/<model>/<pdb_id>/model_NNN.cif`` tree.

Each generative model emits its own directory layout with its own naming and
confidence metadata format. The scoring pipeline wants a clean, uniform input:

    data/
      af3/
        1abc/
          model_000.cif        # highest ranking_score
          model_001.cif
          manifest.json        # which raw file each model came from, plus raw scores
      proteinx/
        1abc/
          model_000.cif        # highest ranking_score
          ...
      chai/
        6dql/
          ...
      dynamicbind/
        6dql/
          model_000.cif        # merged receptor + ligand
          ...

Different models can cover different PDB sets — that's the point of the
model-first layout. Each ``<model>/<pdb_id>/`` directory is self-contained and
carries a ``manifest.json`` documenting where the ranked models originated,
plus all raw per-sample scores for downstream analysis.

Normalization routes through the existing sanitizer so the emitted CIFs are
OST-ready (entity metadata rebuilt, SDF ligands merged as HETATM chains, etc.).
"""
from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator

import numpy as np

from .schema import ModelRecord
from .sanitizer import sanitize

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers (local copies — keeping normalizer independent of registry internals)
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
# Raw-record discovery (one per producer).
# Each discoverer takes the producer's own raw root directory and yields
# (pdb_id, list[ModelRecord]) — flexible about layout within that root.
# ---------------------------------------------------------------------------
Discoverer = Callable[[Path], Iterator[tuple[str, list[ModelRecord]]]]
DISCOVERERS: dict[str, Discoverer] = {}


def register_discoverer(name: str):
    def wrap(fn: Discoverer):
        DISCOVERERS[name] = fn
        return fn
    return wrap


def _iter_pdb_subdirs(root: Path) -> Iterator[tuple[str, Path]]:
    """Yield (pdb_id_lowercased, subdir) for each immediate subdir that looks
    like a PDB ID (4 alphanumeric chars). Skips hidden dirs."""
    for p in sorted(root.iterdir()):
        if not p.is_dir() or p.name.startswith("."):
            continue
        yield p.name.lower(), p


# ─── AF3 ────────────────────────────────────────────────────────────────────
# Expected raw layout:  <raw_root>/<pdb_id>/<any>/model.cif + summary_confidences.json
# (AF3's on-disk tree has seed_*/sample_* nesting; rglob handles any depth.)
@register_discoverer("af3")
def discover_af3(raw_root: Path) -> Iterator[tuple[str, list[ModelRecord]]]:
    for pdb_id, pdb_dir in _iter_pdb_subdirs(raw_root):
        records: list[ModelRecord] = []
        for cif in sorted(pdb_dir.rglob("*.cif")):
            summary = cif.with_name("summary_confidences.json")
            if not summary.exists():
                cands = list(cif.parent.glob("*summary_confidences*.json"))
                summary = cands[0] if cands else None
            meta = _safe_json(summary) if summary else {}
            meta = meta or {}
            raw = {k: float(meta[k]) for k in (
                "ranking_score", "iptm", "ptm", "fraction_disordered",
                "has_clash", "num_recycles",
            ) if k in meta and isinstance(meta[k], (int, float))}
            confidence = raw.get("ranking_score",
                                  raw.get("iptm", raw.get("ptm")))
            records.append(ModelRecord(
                pdb_id=pdb_id, producer="af3", structure_path=cif,
                raw_scores=raw, confidence=confidence,
            ))
        if records:
            yield pdb_id, records


# ─── Protenix / ProteinX ────────────────────────────────────────────────────
# Real Protenix layout per bytedance/Protenix docs:
#   <raw_root>/<pdb_id>/<seed>/<name>_<seed>_sample_N.cif
#   <raw_root>/<pdb_id>/<seed>/<name>_<seed>_summary_confidence_sample_N.json
# Some pipeline wrappers insert a ``predictions/`` folder; we handle both.
@register_discoverer("proteinx")
def discover_proteinx(raw_root: Path) -> Iterator[tuple[str, list[ModelRecord]]]:
    for pdb_id, pdb_dir in _iter_pdb_subdirs(raw_root):
        records: list[ModelRecord] = []
        for cif in sorted(pdb_dir.rglob("*.cif")):
            m = re.search(r"sample[_-]?(\d+)", cif.stem)
            sample_idx = m.group(1) if m else None
            conf_json: Path | None = None
            if sample_idx is not None:
                preferred = list(cif.parent.glob(
                    f"*summary_confidence*sample*{sample_idx}*.json"))
                if preferred:
                    conf_json = preferred[0]
                else:
                    cands = list(cif.parent.glob(f"*sample*{sample_idx}*.json"))
                    if cands:
                        conf_json = cands[0]
            if conf_json is None:
                cands = (list(cif.parent.glob("*summary_confidence*.json"))
                         or list(cif.parent.glob("*confidence*.json")))
                conf_json = cands[0] if cands else None
            meta = _safe_json(conf_json) if conf_json else {}
            meta = meta or {}
            raw: dict[str, float] = {}
            for k in ("ranking_score", "plddt", "gpde", "ptm", "iptm",
                       "has_clash", "disorder",
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
            records.append(ModelRecord(
                pdb_id=pdb_id, producer="proteinx", structure_path=cif,
                raw_scores=raw, confidence=confidence,
            ))
        if records:
            yield pdb_id, records


# ─── Chai-1 ─────────────────────────────────────────────────────────────────
# Expected raw layout:
#   <raw_root>/<pdb_id>/pred.model_idx_<k>.cif + scores.model_idx_<k>.npz
@register_discoverer("chai")
def discover_chai(raw_root: Path) -> Iterator[tuple[str, list[ModelRecord]]]:
    for pdb_id, pdb_dir in _iter_pdb_subdirs(raw_root):
        records: list[ModelRecord] = []
        for cif in sorted(pdb_dir.rglob("pred.model_idx_*.cif")):
            m = re.search(r"model_idx_(\d+)", cif.stem)
            if not m:
                continue
            k = m.group(1)
            npz = cif.parent / f"scores.model_idx_{k}.npz"
            data = _safe_npz(npz) if npz.exists() else {}
            data = data or {}
            raw: dict[str, float] = {}
            for key in ("aggregate_score", "ptm", "iptm",
                         "has_inter_chain_clashes", "chain_chain_clashes"):
                if key in data:
                    val = data[key]
                    raw[key] = _mean(val) if np.ndim(val) > 0 else float(val)
            confidence = raw.get("aggregate_score",
                                  raw.get("iptm", raw.get("ptm")))
            records.append(ModelRecord(
                pdb_id=pdb_id, producer="chai", structure_path=cif,
                raw_scores=raw, confidence=confidence,
            ))
        if records:
            yield pdb_id, records


# ─── DynamicBind ────────────────────────────────────────────────────────────
# Expected raw layout:
#   <raw_root>/<pdb_id>/rank<k>_receptor*.pdb + rank<k>_ligand*.sdf
#                         complex_confidence.csv  (or affinity_prediction.csv)
@register_discoverer("dynamicbind")
def discover_dynamicbind(raw_root: Path) -> Iterator[tuple[str, list[ModelRecord]]]:
    import csv as _csv
    for pdb_id, pdb_dir in _iter_pdb_subdirs(raw_root):
        # CSV with rank → confidence mapping
        conf_map: dict[int, dict[str, float]] = {}
        for csv_name in ("complex_confidence.csv", "affinity_prediction.csv"):
            csv_path = pdb_dir / csv_name
            if not csv_path.exists():
                continue
            try:
                with csv_path.open() as fh:
                    for row in _csv.DictReader(fh):
                        rank_str = (row.get("rank") or row.get("model_rank")
                                    or row.get("index"))
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
                log.warning("[dynamicbind %s] csv parse: %s", pdb_id, e)

        records: list[ModelRecord] = []
        for pdb in sorted(pdb_dir.rglob("rank*_receptor*.pdb")):
            m = re.search(r"rank(\d+)", pdb.name)
            if not m:
                continue
            rank = int(m.group(1))
            sdfs = list(pdb.parent.glob(f"rank{rank}_ligand*.sdf"))
            sdf = sdfs[0] if sdfs else None
            raw = conf_map.get(rank, {})
            confidence = raw.get("confidence",
                                  raw.get("affinity", -float(rank)))
            records.append(ModelRecord(
                pdb_id=pdb_id, producer="dynamicbind",
                structure_path=pdb, ligand_path=sdf,
                raw_scores=raw, confidence=confidence,
            ))
        if records:
            yield pdb_id, records


# ---------------------------------------------------------------------------
# Normalization driver
# ---------------------------------------------------------------------------
@dataclass
class NormalizedEntry:
    """One entry = one (model, pdb_id) group after ranking + materialization."""
    producer: str
    pdb_id: str
    output_dir: Path
    n_models: int
    failed: int = 0
    best_confidence: float | None = None


def _rank_records(records: list[ModelRecord]) -> list[ModelRecord]:
    records.sort(
        key=lambda r: (r.confidence if r.confidence is not None else -float("inf")),
        reverse=True,
    )
    for i, r in enumerate(records):
        r.model_idx = i
    return records


def normalize_one(
    producer: str,
    pdb_id: str,
    records: list[ModelRecord],
    output_root: Path,
    overwrite: bool = True,
) -> NormalizedEntry:
    """Rank records by confidence and materialize ``model_NNN.cif`` + manifest.

    Returns a summary. Individual sanitize failures are logged but don't abort
    the group — failed models simply don't get a file and appear in the manifest
    with ``status = "sanitize_failed"``.
    """
    records = _rank_records(records)
    out_dir = output_root / producer / pdb_id
    out_dir.mkdir(parents=True, exist_ok=True)

    if overwrite:
        # Clean stale model_*.cif so a rerun with fewer models doesn't leave orphans
        for stale in out_dir.glob("model_*.cif"):
            try:
                stale.unlink()
            except OSError:
                pass

    manifest_models: list[dict] = []
    n_ok, n_fail = 0, 0
    for r in records:
        entry = {
            "model_idx": r.model_idx,
            "confidence": r.confidence,
            "source_path": str(r.structure_path),
            "source_ligand_path": str(r.ligand_path) if r.ligand_path else None,
            "raw_scores": r.raw_scores,
            "status": "ok",
            "error": None,
            "canonical_path": None,
        }
        try:
            sanitized = sanitize(
                structure_path=r.structure_path,
                pdb_id=pdb_id,
                ligand_path=r.ligand_path,
            )
            dest = out_dir / f"model_{r.model_idx:03d}.cif"
            sanitized.write(dest)
            entry["canonical_path"] = str(dest)
            entry["n_chains"] = sanitized.n_chains
            entry["n_ligand_residues"] = sanitized.n_ligand_residues
            entry["ligand_chain_ids"] = list(sanitized.ligand_chain_ids)
            n_ok += 1
        except Exception as e:
            log.warning("[%s/%s/model_%03d] sanitize failed: %s",
                         producer, pdb_id, r.model_idx, e)
            entry["status"] = "sanitize_failed"
            entry["error"] = f"{type(e).__name__}: {e}"
            n_fail += 1
        manifest_models.append(entry)

    manifest = {
        "producer": producer,
        "pdb_id": pdb_id,
        "n_models": len(records),
        "n_ok": n_ok,
        "n_failed": n_fail,
        "models": manifest_models,
    }
    with (out_dir / "manifest.json").open("w") as fh:
        json.dump(manifest, fh, indent=2, default=str)

    best_conf = records[0].confidence if records else None
    return NormalizedEntry(
        producer=producer, pdb_id=pdb_id, output_dir=out_dir,
        n_models=n_ok, failed=n_fail, best_confidence=best_conf,
    )


def normalize_producer(
    producer: str,
    raw_root: Path,
    output_root: Path,
    overwrite: bool = True,
) -> list[NormalizedEntry]:
    """Normalize a whole producer: walk raw_root, write ``output_root/<producer>/<pdb>/...``.

    ``raw_root`` is the *producer's own* top-level directory — whatever layout
    the producer dumps. For Protenix that's typically ``raw/proteinx/``
    (containing per-PDB subdirectories).
    """
    raw_root = Path(raw_root)
    output_root = Path(output_root)
    if not raw_root.is_dir():
        raise FileNotFoundError(f"raw_root does not exist: {raw_root}")

    discover = DISCOVERERS.get(producer)
    if discover is None:
        raise KeyError(f"no discoverer registered for producer '{producer}'. "
                        f"Known: {sorted(DISCOVERERS)}")

    entries: list[NormalizedEntry] = []
    for pdb_id, records in discover(raw_root):
        try:
            entry = normalize_one(producer, pdb_id, records, output_root,
                                  overwrite=overwrite)
            entries.append(entry)
            log.info("normalized %s/%s: %d ok, %d failed",
                      producer, pdb_id, entry.n_models, entry.failed)
        except Exception as e:
            log.exception("normalize_one crashed for %s/%s: %s",
                           producer, pdb_id, e)
    return entries


def normalize_all(
    producer_roots: dict[str, Path],
    output_root: Path,
    overwrite: bool = True,
) -> dict[str, list[NormalizedEntry]]:
    """Normalize multiple producers in one call.

    ``producer_roots`` maps producer name → its raw root directory:

    .. code-block:: python

        normalize_all({
            "proteinx": "/data/raw/proteinx",
            "af3":      "/data/raw/af3",
        }, output_root="/data/normalized")
    """
    result: dict[str, list[NormalizedEntry]] = {}
    for producer, raw_root in producer_roots.items():
        result[producer] = normalize_producer(
            producer, Path(raw_root), output_root, overwrite=overwrite)
    return result
