"""End-to-end demo: build a synthetic input tree, run the pipeline, show output.

This script is also a good smoke test for the package as installed.
Patches OST scoring since OpenStructure isn't installed in this environment.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Borrow the synthetic tree builder from the test module
from tests.test_phases_1_to_3 import (  # noqa: E402
    MINI_CIF, MINI_PDB, _tiny_sdf,
)


def build_tree(root: Path) -> Path:
    """Build a realistic-looking multi-PDB input tree."""
    inputs = root / "inputs"
    for pdb_id in ("1abc", "2def"):
        base = inputs / pdb_id

        # AF3: three seeds
        for seed, rank in zip((1, 2, 3), (0.91, 0.87, 0.75)):
            d = base / "af3" / f"seed_{seed}"
            d.mkdir(parents=True, exist_ok=True)
            (d / "model.cif").write_text(MINI_CIF)
            (d / "summary_confidences.json").write_text(json.dumps({
                "ranking_score": rank, "iptm": rank - 0.1, "ptm": rank - 0.2,
            }))

        # ProteinX: two samples
        d = base / "proteinx"
        d.mkdir(parents=True, exist_ok=True)
        for k, rc in enumerate([0.66, 0.71]):
            (d / f"sample_{k}.cif").write_text(MINI_CIF)
            (d / f"confidence_sample_{k}.json").write_text(json.dumps({
                "ranking_confidence": rc, "complex_plddt": 78 + k * 2,
            }))

        # Chai: two models
        d = base / "chai"
        d.mkdir(parents=True, exist_ok=True)
        for k, ag in enumerate([0.45, 0.60]):
            (d / f"pred.model_idx_{k}.cif").write_text(MINI_CIF)
            np.savez(d / f"scores.model_idx_{k}.npz",
                     aggregate_score=np.array([ag]),
                     iptm=np.array([ag - 0.05]), ptm=np.array([ag - 0.1]))

        # DynamicBind: two ranks
        d = base / "dynamicbind"
        d.mkdir(parents=True, exist_ok=True)
        (d / "rank1_receptor_1.pdb").write_text(MINI_PDB)
        (d / "rank1_ligand_1.sdf").write_text(_tiny_sdf())
        (d / "rank2_receptor_1.pdb").write_text(MINI_PDB)
        (d / "rank2_ligand_1.sdf").write_text(_tiny_sdf())
        (d / "complex_confidence.csv").write_text(
            "rank,confidence,affinity\n1,0.82,6.5\n2,0.71,6.1\n"
        )
    return inputs


def fake_score(sanitized, ref_path):
    h = abs(hash(sanitized.text)) % 1000 / 1000.0
    return {"bisy_rmsd": 1.0 + h, "lddt_pli": 0.5 + h / 2, "qs_global": 0.7 + h / 3}


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        inputs = build_tree(tmp)
        refs = tmp / "refs"
        refs.mkdir()
        # Pre-populate both references
        for pid in ("1ABC", "2DEF"):
            (refs / f"{pid}.cif").write_text(MINI_CIF)
        out = tmp / "out"

        # Patch OST since it's not installed in this sandbox, and use a thread
        # pool so the patch survives across "workers".
        from plb_bench.pipeline import run, PipelineConfig
        import plb_bench.pipeline as pipe
        from concurrent.futures import ThreadPoolExecutor

        cfg = PipelineConfig(
            input_root=inputs, refs_dir=refs, output_dir=out,
            n_workers=2, download_missing_refs=False, output_format="parquet",
        )
        with patch("plb_bench.scoring.score_model", side_effect=fake_score), \
             patch.object(pipe, "ProcessPoolExecutor", ThreadPoolExecutor):
            out_path = run(cfg)

        print(f"\n== Output file: {out_path}")
        df = pd.read_parquet(out_path) if out_path.suffix == ".parquet" else pd.read_csv(out_path)
        print(f"== Total rows: {len(df)}")
        print(f"== Status counts:")
        print(df["status"].value_counts().to_string())
        print(f"\n== First few rows:")
        with pd.option_context("display.max_columns", 12, "display.width", 180):
            print(df.head(10).to_string(index=False))
        print(f"\n== Canonical ranked models emitted:")
        for p in sorted(out.rglob("model_*.cif")):
            print(f"  {p.relative_to(out)}")


if __name__ == "__main__":
    main()
