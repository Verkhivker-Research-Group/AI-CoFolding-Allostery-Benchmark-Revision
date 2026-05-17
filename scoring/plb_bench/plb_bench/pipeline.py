"""Top-level orchestrator.

Scheduling unit
---------------
We parallelize at the **(pdb_id, producer)** level, not per-model. A single
PDB has up to ~5 models per producer; keeping them together lets a worker
amortize the reference-loading cost across models.

Error policy
------------
Every worker wraps its work in try/except. A failed model becomes a
``ScoreRecord`` with ``status != "ok"`` — it still appears in the output,
so downstream analysis can see *why* predictions were skipped. A failed
worker never poisons the pool.
"""
from __future__ import annotations

import logging
import os
import shutil
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .schema import ModelRecord, ScoreRecord
from .registry import PARSERS, iter_pdb_dirs, parse_all, rank_and_assign
from .references import get_reference
from .sanitizer import sanitize
from . import scoring as _scoring_mod

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
@dataclass
class PipelineConfig:
    input_root: Path
    refs_dir: Path
    output_dir: Path
    producers: list[str] | None = None       # None -> all registered
    n_workers: int | None = None             # None -> os.cpu_count()
    use_ray: bool = False
    download_missing_refs: bool = True
    write_canonical_models: bool = True      # emit ranked model_000.cif ... copies
    output_format: str = "parquet"           # "parquet" | "csv"


# ---------------------------------------------------------------------------
# Per-task worker
# ---------------------------------------------------------------------------
def _canonical_path(output_dir: Path, pdb_id: str,
                    producer: str, idx: int) -> Path:
    return output_dir / "models" / pdb_id.lower() / producer / f"model_{idx:03d}.cif"


def _process_group(args: tuple) -> list[dict]:
    """Worker entry point. Must be top-level / picklable.

    Processes all models of one (pdb_id, producer) pair:
      1. Sanitize each model (in-memory).
      2. Optionally write canonical model_NNN.cif.
      3. Fetch reference once.
      4. Score each model via OST.
    """
    (pdb_id, producer, records_serialized, refs_dir,
     output_dir, download_missing, write_canonical) = args

    refs_dir = Path(refs_dir)
    output_dir = Path(output_dir)
    logging.basicConfig(level=logging.INFO,
                         format="%(asctime)s [%(process)d] %(levelname)s %(message)s")

    # Rehydrate ModelRecord objects
    records = [ModelRecord(**r) for r in records_serialized]
    for r in records:
        r.structure_path = Path(r.structure_path)
        if r.ligand_path is not None:
            r.ligand_path = Path(r.ligand_path)
        if r.canonical_path is not None:
            r.canonical_path = Path(r.canonical_path)

    out: list[dict] = []

    # Reference: fetch once for the whole group
    ref_path: Path | None = None
    ref_source: str | None = None
    try:
        ref_path, ref_source = get_reference(
            pdb_id, refs_dir, allow_download=download_missing)
    except Exception as e:
        log.error("[%s/%s] reference unavailable: %s", pdb_id, producer, e)
        for r in records:
            out.append(ScoreRecord(
                pdb_id=pdb_id, producer=producer, model_idx=r.model_idx or 0,
                confidence=r.confidence, raw_scores=r.raw_scores,
                status="no_reference", error=str(e),
            ).as_dict())
        return out

    # Score each model
    for r in records:
        idx = r.model_idx if r.model_idx is not None else 0
        try:
            sanitized = sanitize(
                structure_path=r.structure_path,
                pdb_id=pdb_id,
                ligand_path=r.ligand_path,
            )
        except Exception as e:
            log.exception("[%s/%s#%d] sanitize failed", pdb_id, producer, idx)
            out.append(ScoreRecord(
                pdb_id=pdb_id, producer=producer, model_idx=idx,
                confidence=r.confidence, raw_scores=r.raw_scores,
                status="sanitize_failed", error=f"{type(e).__name__}: {e}",
                reference_source=ref_source,
            ).as_dict())
            continue

        # Optionally persist the canonical ranked model on disk
        if write_canonical:
            try:
                dest = _canonical_path(output_dir, pdb_id, producer, idx)
                sanitized.write(dest)
                r.canonical_path = dest
            except Exception as e:
                log.warning("[%s/%s#%d] write canonical failed: %s",
                             pdb_id, producer, idx, e)

        # Score
        try:
            metrics = _scoring_mod.score_model(sanitized, ref_path)
            status = "ok" if any(v is not None for v in metrics.values()) else "ost_no_metrics"
            out.append(ScoreRecord(
                pdb_id=pdb_id, producer=producer, model_idx=idx,
                confidence=r.confidence,
                bisy_rmsd=metrics.get("bisy_rmsd"),
                lddt_pli=metrics.get("lddt_pli"),
                qs_global=metrics.get("qs_global"),
                status=status,
                reference_source=ref_source,
                raw_scores=r.raw_scores,
            ).as_dict())
        except Exception as e:
            tb = traceback.format_exc(limit=3)
            log.error("[%s/%s#%d] scoring crashed:\n%s", pdb_id, producer, idx, tb)
            out.append(ScoreRecord(
                pdb_id=pdb_id, producer=producer, model_idx=idx,
                confidence=r.confidence, raw_scores=r.raw_scores,
                status="ost_failed", error=f"{type(e).__name__}: {e}",
                reference_source=ref_source,
            ).as_dict())
    return out


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def _collect_tasks(cfg: PipelineConfig) -> list[tuple]:
    tasks: list[tuple] = []
    producers = cfg.producers or list(PARSERS.keys())
    for pdb_id, pdb_dir in iter_pdb_dirs(cfg.input_root):
        records = parse_all(pdb_dir, pdb_id, producers=producers)
        if not records:
            log.info("[%s] no producer outputs found — skipping", pdb_id)
            continue
        records = rank_and_assign(records)
        # Group by producer so the worker amortizes reference load
        by_prod: dict[str, list[ModelRecord]] = {}
        for r in records:
            by_prod.setdefault(r.producer, []).append(r)
        for producer, group in by_prod.items():
            tasks.append((
                pdb_id, producer,
                [r.as_dict() for r in group],
                str(cfg.refs_dir), str(cfg.output_dir),
                cfg.download_missing_refs, cfg.write_canonical_models,
            ))
    return tasks


def _write_output(rows: list[dict], cfg: PipelineConfig) -> Path:
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    df = pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["pdb_id", "producer", "model_idx", "confidence",
                 "bisy_rmsd", "lddt_pli", "qs_global", "status"])
    df = df.sort_values(["pdb_id", "producer", "model_idx"], ignore_index=True)

    if cfg.output_format == "parquet":
        dest = cfg.output_dir / "benchmark.parquet"
        try:
            df.to_parquet(dest, index=False)
        except Exception as e:
            log.warning("parquet write failed (%s); falling back to CSV", e)
            dest = cfg.output_dir / "benchmark.csv"
            df.to_csv(dest, index=False)
    else:
        dest = cfg.output_dir / "benchmark.csv"
        df.to_csv(dest, index=False)
    log.info("wrote %d rows -> %s", len(df), dest)
    return dest


def run(cfg: PipelineConfig) -> Path:
    """Run the full pipeline end-to-end and return the output file path."""
    logging.basicConfig(level=logging.INFO,
                         format="%(asctime)s %(levelname)s %(message)s")
    cfg.input_root = Path(cfg.input_root)
    cfg.refs_dir = Path(cfg.refs_dir)
    cfg.output_dir = Path(cfg.output_dir)

    tasks = _collect_tasks(cfg)
    log.info("collected %d (pdb_id, producer) tasks", len(tasks))
    if not tasks:
        return _write_output([], cfg)

    rows: list[dict] = []

    if cfg.use_ray:
        rows = _run_with_ray(tasks)
    else:
        n = cfg.n_workers or os.cpu_count() or 1
        # ProcessPoolExecutor: true CPU parallelism, respects the try/except
        # contract inside each worker.
        with ProcessPoolExecutor(max_workers=n) as ex:
            futures = {ex.submit(_process_group, t): t for t in tasks}
            for fut in as_completed(futures):
                pdb_id, producer = futures[fut][0], futures[fut][1]
                try:
                    rows.extend(fut.result())
                except Exception as e:
                    log.exception("[%s/%s] worker crashed outright", pdb_id, producer)
                    rows.append(ScoreRecord(
                        pdb_id=pdb_id, producer=producer, model_idx=-1,
                        confidence=None,
                        status="worker_crashed", error=f"{type(e).__name__}: {e}",
                    ).as_dict())

    return _write_output(rows, cfg)


def _run_with_ray(tasks: list[tuple]) -> list[dict]:
    import ray
    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True, log_to_driver=False)
    remote_fn = ray.remote(_process_group)
    futures = [remote_fn.remote(t) for t in tasks]
    rows: list[dict] = []
    for f in futures:
        try:
            rows.extend(ray.get(f))
        except Exception as e:
            log.exception("ray task crashed: %s", e)
    return rows
