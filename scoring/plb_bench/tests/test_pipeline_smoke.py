"""End-to-end pipeline smoke test.

We stub out the OST-dependent ``score_model`` so this test can run on any
machine. It verifies: crawl -> parse -> rank -> sanitize -> reference -> score
-> parquet/csv output, plus the robust error-handling path.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Reuse fixtures from the sibling test module
from tests.test_phases_1_to_3 import fake_tree, MINI_CIF, _tiny_sdf  # noqa: F401


@pytest.fixture
def cached_reference(tmp_path: Path) -> Path:
    """Pre-populate the refs dir to avoid any network call."""
    refs = tmp_path / "refs"
    refs.mkdir()
    (refs / "1ABC.cif").write_text(MINI_CIF)
    return refs


def _fake_score(sanitized, ref_path):
    # Deterministic-ish fake metrics keyed off the sanitized CIF length so
    # different models get different numbers.
    h = abs(hash(sanitized.text)) % 1000 / 1000.0
    return {"bisy_rmsd": 1.0 + h, "lddt_pli": 0.5 + h / 2, "qs_global": 0.7 + h / 3}


def test_pipeline_end_to_end(fake_tree, cached_reference, tmp_path):
    from plb_bench.pipeline import run, PipelineConfig

    out_dir = tmp_path / "out"
    cfg = PipelineConfig(
        input_root=fake_tree, refs_dir=cached_reference, output_dir=out_dir,
        n_workers=2, download_missing_refs=False, output_format="csv",
    )

    # Patch inside the worker module path so the subprocess picks it up.
    # Use a single process for the patch to take effect.
    cfg.n_workers = 1
    with patch("plb_bench.scoring.score_model", side_effect=_fake_score):
        # Force in-process execution: override ProcessPoolExecutor with a
        # thread pool for testability.
        import plb_bench.pipeline as pipe
        from concurrent.futures import ThreadPoolExecutor

        with patch.object(pipe, "ProcessPoolExecutor", ThreadPoolExecutor):
            path = run(cfg)

    assert path.exists()
    df = pd.read_csv(path)
    # 4 producers x 2 models = 8 rows
    assert len(df) == 8
    assert set(df["producer"]) == {"af3", "proteinx", "chai", "dynamicbind"}
    # Each producer ranked 0 and 1
    for prod, grp in df.groupby("producer"):
        assert sorted(grp["model_idx"]) == [0, 1]
    # All metrics populated (our fake always returns numbers)
    assert df["bisy_rmsd"].notna().all()
    assert df["lddt_pli"].notna().all()
    assert df["qs_global"].notna().all()
    # Status clean
    assert (df["status"] == "ok").all()


def test_pipeline_handles_missing_reference(fake_tree, tmp_path):
    """If reference is missing and download is disabled, every model should
    emit a row with status=no_reference rather than crashing the worker."""
    from plb_bench.pipeline import run, PipelineConfig
    import plb_bench.pipeline as pipe
    from concurrent.futures import ThreadPoolExecutor

    out_dir = tmp_path / "out"
    refs_empty = tmp_path / "refs_empty"
    refs_empty.mkdir()

    cfg = PipelineConfig(
        input_root=fake_tree, refs_dir=refs_empty, output_dir=out_dir,
        n_workers=1, download_missing_refs=False, output_format="csv",
    )
    with patch.object(pipe, "ProcessPoolExecutor", ThreadPoolExecutor):
        path = run(cfg)
    df = pd.read_csv(path)
    assert len(df) == 8
    assert (df["status"] == "no_reference").all()
    assert df["bisy_rmsd"].isna().all()


def test_pipeline_survives_sanitize_failure(fake_tree, cached_reference, tmp_path):
    """A corrupt model in one producer must not take down the other rows."""
    from plb_bench.pipeline import run, PipelineConfig
    import plb_bench.pipeline as pipe
    from concurrent.futures import ThreadPoolExecutor

    # Nuke the content of one AF3 model to garbage
    bad = next((fake_tree / "1abc" / "af3").rglob("*.cif"))
    bad.write_text("this is not a cif\n")

    cfg = PipelineConfig(
        input_root=fake_tree, refs_dir=cached_reference, output_dir=tmp_path / "out",
        n_workers=1, download_missing_refs=False, output_format="csv",
    )
    with patch("plb_bench.scoring.score_model", side_effect=_fake_score), \
         patch.object(pipe, "ProcessPoolExecutor", ThreadPoolExecutor):
        path = run(cfg)
    df = pd.read_csv(path)
    # The broken AF3 model still produces a row (status=sanitize_failed);
    # the other 7 rows are fine.
    assert len(df) == 8
    bad_rows = df[df["status"] == "sanitize_failed"]
    ok_rows = df[df["status"] == "ok"]
    assert len(bad_rows) == 1
    assert len(ok_rows) == 7
    assert bad_rows.iloc[0]["producer"] == "af3"


def test_canonical_model_files_written(fake_tree, cached_reference, tmp_path):
    from plb_bench.pipeline import run, PipelineConfig
    import plb_bench.pipeline as pipe
    from concurrent.futures import ThreadPoolExecutor

    out_dir = tmp_path / "out"
    cfg = PipelineConfig(
        input_root=fake_tree, refs_dir=cached_reference, output_dir=out_dir,
        n_workers=1, download_missing_refs=False, output_format="csv",
        write_canonical_models=True,
    )
    with patch("plb_bench.scoring.score_model", side_effect=_fake_score), \
         patch.object(pipe, "ProcessPoolExecutor", ThreadPoolExecutor):
        run(cfg)
    models_dir = out_dir / "models" / "1abc"
    for producer in ("af3", "proteinx", "chai", "dynamicbind"):
        assert (models_dir / producer / "model_000.cif").exists()
        assert (models_dir / producer / "model_001.cif").exists()
