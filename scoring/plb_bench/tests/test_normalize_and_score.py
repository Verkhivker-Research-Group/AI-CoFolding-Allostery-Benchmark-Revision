"""Tests for the normalize → score_tree flow."""
import json
import sys
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from plb_bench.normalize import (
    DISCOVERERS, discover_proteinx, normalize_producer, normalize_one,
    normalize_all,
)
from plb_bench.score_tree import ScoreConfig, run as score_run


# Reuse fixtures / fake_tree from the other test module
from tests.test_phases_1_to_3 import (  # noqa: F401
    fake_tree, MINI_CIF, MINI_PDB, _tiny_sdf,
)


# -----------------------------------------------------------------------------
# Discovery against a real-layout Protenix tree
# -----------------------------------------------------------------------------
def _build_proteinx_raw(tmp: Path, pdb_id: str, scores: list[float]) -> Path:
    """Build raw/<pdb_id>/seed_101/predictions/ matching bytedance/Protenix layout."""
    pred = tmp / pdb_id / "seed_101" / "predictions"
    pred.mkdir(parents=True)
    for i, s in enumerate(scores):
        (pred / f"{pdb_id}_101_sample_{i}.cif").write_text(MINI_CIF)
        (pred / f"{pdb_id}_101_summary_confidence_sample_{i}.json").write_text(json.dumps({
            "ranking_score": s, "plddt": 70 + i, "iptm": 0.5 + i * 0.1,
            "ptm": 0.55, "has_clash": False,
        }))
    return pred


def test_discover_proteinx_real_layout(tmp_path: Path):
    raw_root = tmp_path / "raw" / "proteinx"
    raw_root.mkdir(parents=True)
    _build_proteinx_raw(raw_root, "1abc", [0.66, 0.81])
    _build_proteinx_raw(raw_root, "6dql", [0.55, 0.72, 0.90])

    groups = dict(discover_proteinx(raw_root))
    assert set(groups) == {"1abc", "6dql"}
    assert len(groups["1abc"]) == 2
    assert len(groups["6dql"]) == 3
    # Raw scores captured correctly
    for r in groups["1abc"]:
        assert "ranking_score" in r.raw_scores
        assert "plddt" in r.raw_scores


# -----------------------------------------------------------------------------
# Normalization writes model_NNN.cif + manifest
# -----------------------------------------------------------------------------
def test_normalize_producer_proteinx(tmp_path: Path):
    raw_root = tmp_path / "raw" / "proteinx"
    raw_root.mkdir(parents=True)
    _build_proteinx_raw(raw_root, "1abc", [0.66, 0.81])
    _build_proteinx_raw(raw_root, "6dql", [0.55, 0.90])
    out_root = tmp_path / "data"

    entries = normalize_producer("proteinx", raw_root, out_root)
    assert len(entries) == 2

    # Structure: data/proteinx/1abc/model_000.cif, model_001.cif, manifest.json
    for pdb_id in ("1abc", "6dql"):
        d = out_root / "proteinx" / pdb_id
        assert (d / "model_000.cif").exists()
        assert (d / "model_001.cif").exists()
        assert (d / "manifest.json").exists()

    # Ranking correctness: 1abc's 0.81 sample should be model_000
    manifest = json.loads((out_root / "proteinx" / "1abc" / "manifest.json").read_text())
    best = next(m for m in manifest["models"] if m["model_idx"] == 0)
    assert best["confidence"] == pytest.approx(0.81)
    assert best["status"] == "ok"
    # raw_scores preserved
    assert "ranking_score" in best["raw_scores"]


def test_normalize_overwrite_clears_stale(tmp_path: Path):
    """If a rerun produces fewer models, old model_*.cif should be removed."""
    raw_root = tmp_path / "raw" / "proteinx"
    raw_root.mkdir(parents=True)
    _build_proteinx_raw(raw_root, "1abc", [0.1, 0.2, 0.3, 0.4])
    out_root = tmp_path / "data"

    normalize_producer("proteinx", raw_root, out_root)
    first_files = set((out_root / "proteinx" / "1abc").glob("model_*.cif"))
    assert len(first_files) == 4

    # Shrink the raw tree to 2 samples and renormalize
    import shutil
    shutil.rmtree(raw_root)
    raw_root.mkdir(parents=True)
    _build_proteinx_raw(raw_root, "1abc", [0.5, 0.6])
    normalize_producer("proteinx", raw_root, out_root, overwrite=True)

    after = set((out_root / "proteinx" / "1abc").glob("model_*.cif"))
    assert len(after) == 2


def test_normalize_survives_single_sanitize_failure(tmp_path: Path):
    """A broken CIF in a group shouldn't kill the rest of the group."""
    raw_root = tmp_path / "raw" / "proteinx"
    raw_root.mkdir(parents=True)
    pred = _build_proteinx_raw(raw_root, "1abc", [0.3, 0.6, 0.9])

    # Corrupt the middle sample
    (pred / "1abc_101_sample_1.cif").write_text("this is not a cif\n")
    out_root = tmp_path / "data"

    entries = normalize_producer("proteinx", raw_root, out_root)
    assert len(entries) == 1
    entry = entries[0]
    # 2 ok + 1 failed
    assert entry.n_models == 2
    assert entry.failed == 1

    manifest = json.loads((out_root / "proteinx" / "1abc" / "manifest.json").read_text())
    statuses = sorted(m["status"] for m in manifest["models"])
    assert statuses == ["ok", "ok", "sanitize_failed"]
    # Only 2 CIFs on disk
    assert len(list((out_root / "proteinx" / "1abc").glob("model_*.cif"))) == 2


# -----------------------------------------------------------------------------
# Scoring the normalized tree
# -----------------------------------------------------------------------------
def _mock_score(san, ref):
    h = abs(hash(san.text)) % 1000 / 1000.0
    return {"bisy_rmsd": 1.0 + h, "lddt_pli": 0.5 + h / 2, "qs_global": 0.7 + h / 3}


def test_score_normalized_tree_single_producer(tmp_path: Path):
    # Build + normalize a Protenix-only tree
    raw_root = tmp_path / "raw" / "proteinx"
    raw_root.mkdir(parents=True)
    _build_proteinx_raw(raw_root, "1abc", [0.66, 0.81])
    _build_proteinx_raw(raw_root, "6dql", [0.55, 0.90])
    data_root = tmp_path / "data"
    normalize_producer("proteinx", raw_root, data_root)

    # Pre-populate refs
    refs = tmp_path / "refs"
    refs.mkdir()
    (refs / "1ABC.cif").write_text(MINI_CIF)
    (refs / "6DQL.cif").write_text(MINI_CIF)

    out_dir = tmp_path / "out"
    cfg = ScoreConfig(
        normalized_root=data_root, refs_dir=refs, output_dir=out_dir,
        n_workers=1, download_missing_refs=False, output_format="csv",
    )

    # Score with mock. Also patch ProcessPoolExecutor to ThreadPoolExecutor
    # so the in-process mock sticks.
    import plb_bench.score_tree as st
    import plb_bench.scoring as scoring
    from concurrent.futures import ThreadPoolExecutor

    with patch.object(scoring, "score_model", side_effect=_mock_score), \
         patch.object(st, "ProcessPoolExecutor", ThreadPoolExecutor):
        out_path = score_run(cfg)

    df = pd.read_csv(out_path)
    # 2 PDBs × 2 models = 4 rows, all Protenix
    assert len(df) == 4
    assert set(df["producer"]) == {"proteinx"}
    assert (df["status"] == "ok").all()
    # Confidence copied through from manifest
    assert df["confidence"].notna().all()
    # Ranked: model_idx 0 has higher confidence than model_idx 1 for each PDB
    for pid, grp in df.groupby("pdb_id"):
        top = grp[grp["model_idx"] == 0].iloc[0]
        second = grp[grp["model_idx"] == 1].iloc[0]
        assert top["confidence"] > second["confidence"]


def test_score_handles_missing_reference(tmp_path: Path):
    raw_root = tmp_path / "raw" / "proteinx"
    raw_root.mkdir(parents=True)
    _build_proteinx_raw(raw_root, "1abc", [0.5, 0.8])
    data_root = tmp_path / "data"
    normalize_producer("proteinx", raw_root, data_root)

    refs_empty = tmp_path / "refs"
    refs_empty.mkdir()

    cfg = ScoreConfig(
        normalized_root=data_root, refs_dir=refs_empty,
        output_dir=tmp_path / "out",
        n_workers=1, download_missing_refs=False, output_format="csv",
    )
    import plb_bench.score_tree as st
    from concurrent.futures import ThreadPoolExecutor
    with patch.object(st, "ProcessPoolExecutor", ThreadPoolExecutor):
        out_path = score_run(cfg)

    df = pd.read_csv(out_path)
    assert len(df) == 2
    assert (df["status"] == "no_reference").all()


def test_score_pdb_filter(tmp_path: Path):
    raw_root = tmp_path / "raw" / "proteinx"
    raw_root.mkdir(parents=True)
    _build_proteinx_raw(raw_root, "1abc", [0.5, 0.8])
    _build_proteinx_raw(raw_root, "6dql", [0.4, 0.7])
    data_root = tmp_path / "data"
    normalize_producer("proteinx", raw_root, data_root)

    refs = tmp_path / "refs"
    refs.mkdir()
    for pid in ("1ABC", "6DQL"):
        (refs / f"{pid}.cif").write_text(MINI_CIF)

    cfg = ScoreConfig(
        normalized_root=data_root, refs_dir=refs,
        output_dir=tmp_path / "out",
        pdb_ids=["6dql"],   # score just one PDB
        n_workers=1, download_missing_refs=False, output_format="csv",
    )
    import plb_bench.score_tree as st
    import plb_bench.scoring as scoring
    from concurrent.futures import ThreadPoolExecutor
    with patch.object(scoring, "score_model", side_effect=_mock_score), \
         patch.object(st, "ProcessPoolExecutor", ThreadPoolExecutor):
        out_path = score_run(cfg)

    df = pd.read_csv(out_path)
    assert len(df) == 2
    assert set(df["pdb_id"]) == {"6dql"}


# -----------------------------------------------------------------------------
# Different models on different PDBs — the scenario the user cares about
# -----------------------------------------------------------------------------
def test_different_producers_different_pdbs(tmp_path: Path):
    """proteinx covers 1abc+2def; chai covers 2def+3ghi. Make sure scoring
    reports exactly those (producer, pdb_id) pairs."""
    # Build proteinx raw with 1abc + 2def
    px_raw = tmp_path / "raw" / "proteinx"
    px_raw.mkdir(parents=True)
    _build_proteinx_raw(px_raw, "1abc", [0.5, 0.8])
    _build_proteinx_raw(px_raw, "2def", [0.6, 0.9])

    # Build chai raw with 2def + 3ghi
    import numpy as np
    chai_raw = tmp_path / "raw" / "chai"
    chai_raw.mkdir(parents=True)
    for pdb_id in ("2def", "3ghi"):
        d = chai_raw / pdb_id
        d.mkdir()
        for k, ag in enumerate([0.45, 0.60]):
            (d / f"pred.model_idx_{k}.cif").write_text(MINI_CIF)
            np.savez(d / f"scores.model_idx_{k}.npz",
                     aggregate_score=np.array([ag]))

    data_root = tmp_path / "data"
    normalize_producer("proteinx", px_raw, data_root)
    normalize_producer("chai",      chai_raw, data_root)

    # Check the tree shape
    assert (data_root / "proteinx" / "1abc" / "model_000.cif").exists()
    assert (data_root / "proteinx" / "2def" / "model_000.cif").exists()
    assert not (data_root / "proteinx" / "3ghi").exists()
    assert not (data_root / "chai" / "1abc").exists()
    assert (data_root / "chai" / "2def" / "model_000.cif").exists()
    assert (data_root / "chai" / "3ghi" / "model_000.cif").exists()

    # Score everything
    refs = tmp_path / "refs"
    refs.mkdir()
    for pid in ("1ABC", "2DEF", "3GHI"):
        (refs / f"{pid}.cif").write_text(MINI_CIF)

    cfg = ScoreConfig(
        normalized_root=data_root, refs_dir=refs,
        output_dir=tmp_path / "out",
        n_workers=1, download_missing_refs=False, output_format="csv",
    )
    import plb_bench.score_tree as st
    import plb_bench.scoring as scoring
    from concurrent.futures import ThreadPoolExecutor
    with patch.object(scoring, "score_model", side_effect=_mock_score), \
         patch.object(st, "ProcessPoolExecutor", ThreadPoolExecutor):
        out_path = score_run(cfg)

    df = pd.read_csv(out_path)
    # proteinx: 2 pdbs × 2 models = 4; chai: 2 pdbs × 2 models = 4. Total 8.
    assert len(df) == 8
    pairs = set(zip(df["producer"], df["pdb_id"]))
    expected = {("proteinx", "1abc"), ("proteinx", "2def"),
                ("chai", "2def"), ("chai", "3ghi")}
    assert pairs == expected
