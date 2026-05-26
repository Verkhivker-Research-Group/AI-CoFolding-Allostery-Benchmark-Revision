"""Scoring pipeline that reads the normalized tree produced by ``normalize.py``.

Layout expected::

    data/
      <producer>/
        <pdb_id>/
          model_000.cif
          model_001.cif
          manifest.json

This is a drop-in alternative to the raw-crawl ``pipeline.run``. It does
exactly the work of Phase 4 (reference retrieval + OST scoring) because
Phases 1, 2, 3 (parse, rank, sanitize) have already run during normalization.

Key benefit: the scorer doesn't care which producer is missing PDB X; each
``<producer>/<pdb_id>/`` directory is independent and self-describing. You can
add, remove, or re-normalize any single (producer, pdb_id) group and rescore
just that slice.
"""
from __future__ import annotations

import json
import logging
import os
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .schema import ScoreRecord
from .references import get_reference, read_reference_text
from . import scoring as _scoring_mod

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
@dataclass
class ScoreConfig:
    normalized_root: Path              # data/ directory from normalize_all
    refs_dir: Path
    output_dir: Path
    producers: list[str] | None = None      # None → every <producer> subdir present
    pdb_ids: list[str] | None = None        # None → all found
    n_workers: int | None = None
    use_ray: bool = False
    download_missing_refs: bool = True
    output_format: str = "parquet"


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------
def _iter_groups(cfg: ScoreConfig) -> Iterable[tuple[str, str, Path]]:
    """Yield (producer, pdb_id, dir_path) for each (producer, pdb_id) in the tree."""
    root = Path(cfg.normalized_root)
    producer_filter = set(cfg.producers) if cfg.producers else None
    pdb_filter = {p.lower() for p in cfg.pdb_ids} if cfg.pdb_ids else None

    if not root.is_dir():
        raise FileNotFoundError(f"normalized_root not a directory: {root}")

    for prod_dir in sorted(root.iterdir()):
        if not prod_dir.is_dir() or prod_dir.name.startswith("."):
            continue
        producer = prod_dir.name
        if producer_filter and producer not in producer_filter:
            continue
        for pdb_dir in sorted(prod_dir.iterdir()):
            if not pdb_dir.is_dir() or pdb_dir.name.startswith("."):
                continue
            pdb_id = pdb_dir.name.lower()
            if pdb_filter and pdb_id not in pdb_filter:
                continue
            yield producer, pdb_id, pdb_dir


# ---------------------------------------------------------------------------
# Worker (one per (producer, pdb_id) group — amortize reference load)
# ---------------------------------------------------------------------------
def _score_group(args: tuple) -> list[dict]:
    (producer, pdb_id, pdb_dir, refs_dir, download_missing_refs) = args
    pdb_dir = Path(pdb_dir)
    refs_dir = Path(refs_dir)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(process)d] %(levelname)s %(message)s",
    )

    # Load manifest for raw_scores + confidence metadata
    manifest_path = pdb_dir / "manifest.json"
    manifest = {}
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
        except Exception as e:
            log.warning("[%s/%s] bad manifest: %s", producer, pdb_id, e)
    by_idx = {m["model_idx"]: m for m in manifest.get("models", [])}

    rows: list[dict] = []

    # Reference: fetch once for the group
    try:
        ref_path, ref_source = get_reference(
            pdb_id, refs_dir, allow_download=download_missing_refs)
    except Exception as e:
        log.error("[%s/%s] reference unavailable: %s", producer, pdb_id, e)
        for cif in sorted(pdb_dir.glob("model_*.cif")):
            idx = _parse_idx(cif)
            m = by_idx.get(idx, {})
            rows.append(ScoreRecord(
                pdb_id=pdb_id, producer=producer, model_idx=idx,
                confidence=m.get("confidence"),
                raw_scores=m.get("raw_scores", {}) or {},
                status="no_reference", error=str(e),
            ).as_dict())
        return rows

    cif_files = sorted(pdb_dir.glob("model_*.cif"))
    if not cif_files:
        log.warning("[%s/%s] no model_*.cif in %s", producer, pdb_id, pdb_dir)
        return rows

    for cif in cif_files:
        idx = _parse_idx(cif)
        m = by_idx.get(idx, {})
        conf = m.get("confidence")
        raw_scores = m.get("raw_scores", {}) or {}

        try:
            # We don't need to sanitize — the normalized CIFs were sanitized at
            # normalization time. Just wrap the file as a SanitizedCIF for the
            # scorer's signature.
            from .sanitizer import SanitizedCIF
            san = SanitizedCIF(
                text=cif.read_text(),
                pdb_id=pdb_id.upper(),
                n_chains=m.get("n_chains", 0),
                n_ligand_residues=m.get("n_ligand_residues", 0),
                ligand_chain_ids=tuple(m.get("ligand_chain_ids", [])),
            )
            # For producers that store a separate ligand SDF (e.g. DynamicBind),
            # pass it to score_model so it can load the ligand with proper bond
            # connectivity (OST cannot infer bonds for unknown residue names like
            # LIG from the CIF alone).  For CIF-native producers this is None and
            # the existing path is unchanged.
            sdf_path = m.get("source_ligand_path")
            metrics = _scoring_mod.score_model(san, ref_path,
                                               ligand_sdf_path=sdf_path)
            status = "ok" if any(v is not None for v in metrics.values()) \
                     else "ost_no_metrics"
            rows.append(ScoreRecord(
                pdb_id=pdb_id, producer=producer, model_idx=idx,
                confidence=conf,
                bisy_rmsd=metrics.get("bisy_rmsd"),
                lddt_pli=metrics.get("lddt_pli"),
                qs_global=metrics.get("qs_global"),
                status=status, reference_source=ref_source,
                raw_scores=raw_scores,
            ).as_dict())
        except Exception as e:
            tb = traceback.format_exc(limit=3)
            log.error("[%s/%s/model_%03d] scoring crashed:\n%s",
                      producer, pdb_id, idx, tb)
            rows.append(ScoreRecord(
                pdb_id=pdb_id, producer=producer, model_idx=idx,
                confidence=conf, raw_scores=raw_scores,
                status="ost_failed", error=f"{type(e).__name__}: {e}",
                reference_source=ref_source,
            ).as_dict())
    return rows


def _parse_idx(cif_path: Path) -> int:
    # model_000.cif → 0
    try:
        return int(cif_path.stem.split("_")[-1])
    except (ValueError, IndexError):
        return -1


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------
def _write_output(rows: list[dict], cfg: ScoreConfig) -> Path:
    import pandas as pd
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    df = (pd.DataFrame(rows) if rows
          else pd.DataFrame(columns=[
              "pdb_id", "producer", "model_idx", "confidence",
              "bisy_rmsd", "lddt_pli", "qs_global", "status"]))
    df = df.sort_values(["producer", "pdb_id", "model_idx"], ignore_index=True)

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


def run(cfg: ScoreConfig) -> Path:
    logging.basicConfig(level=logging.INFO,
                         format="%(asctime)s %(levelname)s %(message)s")
    cfg.normalized_root = Path(cfg.normalized_root)
    cfg.refs_dir = Path(cfg.refs_dir)
    cfg.output_dir = Path(cfg.output_dir)

    tasks = [
        (prod, pid, str(pdir), str(cfg.refs_dir), cfg.download_missing_refs)
        for prod, pid, pdir in _iter_groups(cfg)
    ]
    log.info("collected %d (producer, pdb_id) tasks from %s",
              len(tasks), cfg.normalized_root)
    if not tasks:
        return _write_output([], cfg)

    rows: list[dict] = []

    if cfg.use_ray:
        rows = _run_with_ray(tasks)
    else:
        n = cfg.n_workers or os.cpu_count() or 1
        with ProcessPoolExecutor(max_workers=n) as ex:
            futures = {ex.submit(_score_group, t): t for t in tasks}
            for fut in as_completed(futures):
                prod, pid = futures[fut][0], futures[fut][1]
                try:
                    rows.extend(fut.result())
                except Exception as e:
                    log.exception("[%s/%s] worker crashed: %s", prod, pid, e)
                    rows.append(ScoreRecord(
                        pdb_id=pid, producer=prod, model_idx=-1,
                        confidence=None,
                        status="worker_crashed",
                        error=f"{type(e).__name__}: {e}",
                    ).as_dict())

    return _write_output(rows, cfg)


def _run_with_ray(tasks: list[tuple]) -> list[dict]:
    import ray
    if not ray.is_initialized():
        ray.init(ignore_reinit_error=True, log_to_driver=False)
    remote_fn = ray.remote(_score_group)
    futures = [remote_fn.remote(t) for t in tasks]
    rows: list[dict] = []
    for f in futures:
        try:
            rows.extend(ray.get(f))
        except Exception as e:
            log.exception("ray task crashed: %s", e)
    return rows
