"""Unified data schema for the pipeline.

Every model producer (AF3, Proteinx, Chai, DynamicBind) emits radically different
confidence metadata. We normalize everything to a single ``ModelRecord`` before
ranking and scoring. Downstream stages only ever see this shape.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any


# Canonical score key used for ranking. Producers populate this from whatever
# their native "best overall" metric is (ranking_score / model_confidence / etc.).
CANONICAL_RANK_KEY = "confidence"


@dataclass
class ModelRecord:
    """One predicted complex (pre-ranking)."""
    pdb_id: str
    producer: str                    # "af3" | "proteinx" | "chai" | "dynamicbind"
    structure_path: Path             # .cif / .pdb / .sdf — whatever the tool emitted
    raw_scores: dict[str, float] = field(default_factory=dict)
    # populated during Phase 1 standardization:
    confidence: float | None = None   # the canonical rank key (higher = better)
    ligand_path: Path | None = None   # only set when ligand is a separate SDF
    # populated during Phase 1 ranking:
    model_idx: int | None = None     # 0 = best, 1 = next, ...
    canonical_path: Path | None = None  # where model_000.cif was written

    def as_dict(self) -> dict[str, Any]:
        d = asdict(self)
        for k in ("structure_path", "ligand_path", "canonical_path"):
            if d.get(k) is not None:
                d[k] = str(d[k])
        return d


@dataclass
class ScoreRecord:
    """One row of the final output table (one per scored model)."""
    pdb_id: str
    producer: str
    model_idx: int
    confidence: float | None
    bisy_rmsd: float | None = None
    lddt_pli: float | None = None
    qs_global: float | None = None
    # bookkeeping
    status: str = "ok"               # "ok" | "sanitize_failed" | "ost_failed" | "no_reference" | ...
    error: str | None = None
    reference_source: str | None = None  # "local" | "rcsb"
    # keep original scores around for later analysis
    raw_scores: dict[str, float] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        d = asdict(self)
        # flatten raw_scores into prefixed columns for a flat CSV/Parquet table
        raw = d.pop("raw_scores", {}) or {}
        for k, v in raw.items():
            d[f"raw_{k}"] = v
        return d
