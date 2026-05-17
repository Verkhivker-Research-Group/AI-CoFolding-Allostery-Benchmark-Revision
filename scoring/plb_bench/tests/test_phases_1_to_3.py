"""Unit tests for phases 1-3 (schema, parsers, sanitizer, references).

Phase 4 (OST) requires OpenStructure to be installed separately; those tests
are in test_pipeline_smoke.py and are skipped when ost is unavailable.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

# Allow `import plb_bench` without install
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from plb_bench.schema import ModelRecord, ScoreRecord
from plb_bench.registry import (
    PARSERS, parse_af3, parse_chai, parse_dynamicbind, parse_proteinx,
    rank_and_assign, iter_pdb_dirs,
)


# ---------------------------------------------------------------------------
# Fixtures: synthetic producer output directories
# ---------------------------------------------------------------------------
MINI_PDB = """\
HEADER    TEST
ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.104  13.207  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      12.500  14.500  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      13.500  14.700  10.000  1.00 20.00           O
ATOM      5  N   GLY A   2      11.800  15.500  10.000  1.00 20.00           N
ATOM      6  CA  GLY A   2      12.000  16.800  10.000  1.00 20.00           C
ATOM      7  C   GLY A   2      13.500  17.000  10.000  1.00 20.00           C
ATOM      8  O   GLY A   2      14.000  18.100  10.000  1.00 20.00           O
HETATM    9  C1  LIG B   1      15.000  15.000  10.000  1.00 20.00           C
HETATM   10  O1  LIG B   1      16.200  15.100  10.000  1.00 20.00           O
END
"""

MINI_CIF = """\
data_TEST
#
_entry.id TEST
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N ALA A A 1 1 11.104 13.207 10.000 1.00 20.00
ATOM 2 C CA ALA A A 1 1 12.104 13.207 10.000 1.00 20.00
ATOM 3 C C ALA A A 1 1 12.500 14.500 10.000 1.00 20.00
ATOM 4 O O ALA A A 1 1 13.500 14.700 10.000 1.00 20.00
ATOM 5 N N GLY A A 2 2 11.800 15.500 10.000 1.00 20.00
ATOM 6 C CA GLY A A 2 2 12.000 16.800 10.000 1.00 20.00
ATOM 7 C C GLY A A 2 2 13.500 17.000 10.000 1.00 20.00
ATOM 8 O O GLY A A 2 2 14.000 18.100 10.000 1.00 20.00
HETATM 9 C C1 LIG B B . 1 15.000 15.000 10.000 1.00 20.00
HETATM 10 O O1 LIG B B . 1 16.200 15.100 10.000 1.00 20.00
#
"""


@pytest.fixture
def fake_tree(tmp_path: Path) -> Path:
    """Build a fake input tree with outputs from all four producers."""
    root = tmp_path / "inputs"
    pdb_id = "1abc"
    base = root / pdb_id

    # -- AF3 -----------------------------------------------------------------
    af3 = base / "af3" / "seed_1"
    af3.mkdir(parents=True)
    (af3 / "model.cif").write_text(MINI_CIF)
    (af3 / "summary_confidences.json").write_text(json.dumps({
        "ranking_score": 0.87, "iptm": 0.72, "ptm": 0.65,
        "fraction_disordered": 0.02, "has_clash": 0,
    }))
    af3_2 = base / "af3" / "seed_2"
    af3_2.mkdir(parents=True)
    (af3_2 / "model.cif").write_text(MINI_CIF)
    (af3_2 / "summary_confidences.json").write_text(json.dumps({
        "ranking_score": 0.91, "iptm": 0.80, "ptm": 0.70,
    }))

    # -- ProteinX / Protenix (real file naming per bytedance/Protenix docs) --
    # <name>/<seed>/<name>_<seed>_sample_N.cif
    # <name>/<seed>/<name>_<seed>_summary_confidence_sample_N.json
    px = base / "proteinx" / pdb_id / "seed_101"
    px.mkdir(parents=True)
    (px / f"{pdb_id}_101_sample_0.cif").write_text(MINI_CIF)
    (px / f"{pdb_id}_101_summary_confidence_sample_0.json").write_text(json.dumps({
        "ranking_score": 0.66, "plddt": 78.2, "iptm": 0.55, "ptm": 0.60,
        "has_clash": False,
    }))
    (px / f"{pdb_id}_101_sample_1.cif").write_text(MINI_CIF)
    (px / f"{pdb_id}_101_summary_confidence_sample_1.json").write_text(json.dumps({
        "ranking_score": 0.71, "plddt": 80.1, "iptm": 0.60, "ptm": 0.63,
        "has_clash": False,
    }))

    # -- Chai-1 --------------------------------------------------------------
    chai = base / "chai"
    chai.mkdir(parents=True)
    (chai / "pred.model_idx_0.cif").write_text(MINI_CIF)
    np.savez(chai / "scores.model_idx_0.npz",
              aggregate_score=np.array([0.45]),
              iptm=np.array([0.5]), ptm=np.array([0.55]))
    (chai / "pred.model_idx_1.cif").write_text(MINI_CIF)
    np.savez(chai / "scores.model_idx_1.npz",
              aggregate_score=np.array([0.60]),
              iptm=np.array([0.58]), ptm=np.array([0.62]))

    # -- DynamicBind ---------------------------------------------------------
    db = base / "dynamicbind"
    db.mkdir(parents=True)
    (db / "rank1_receptor_1.pdb").write_text(MINI_PDB)
    (db / "rank1_ligand_1.sdf").write_text(_tiny_sdf())
    (db / "rank2_receptor_1.pdb").write_text(MINI_PDB)
    (db / "rank2_ligand_1.sdf").write_text(_tiny_sdf())
    (db / "complex_confidence.csv").write_text(
        "rank,confidence,affinity\n1,0.82,6.5\n2,0.71,6.1\n"
    )
    return root


def _tiny_sdf() -> str:
    # Two-atom SDF: just C1-O1, enough for the sanitizer to exercise RDKit.
    return (
        "LIG\n"
        "     RDKit          3D\n\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "   15.0000   15.0000   10.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   16.2000   15.1000   10.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0\n"
        "M  END\n"
        "$$$$\n"
    )


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------
def test_model_record_serializes_paths():
    r = ModelRecord(pdb_id="1abc", producer="af3",
                    structure_path=Path("/tmp/x.cif"),
                    raw_scores={"iptm": 0.5}, confidence=0.5, model_idx=0)
    d = r.as_dict()
    assert d["structure_path"] == "/tmp/x.cif"
    assert d["raw_scores"] == {"iptm": 0.5}


def test_score_record_flattens_raw_scores():
    sr = ScoreRecord(
        pdb_id="1abc", producer="af3", model_idx=0, confidence=0.9,
        bisy_rmsd=1.2, lddt_pli=0.8, qs_global=0.95,
        raw_scores={"iptm": 0.7, "ptm": 0.65},
    )
    d = sr.as_dict()
    assert "raw_iptm" in d and d["raw_iptm"] == 0.7
    assert "raw_ptm" in d and d["raw_ptm"] == 0.65
    assert "raw_scores" not in d
    assert d["bisy_rmsd"] == 1.2


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------
def test_af3_parser(fake_tree: Path):
    pdb_dir = fake_tree / "1abc"
    recs = parse_af3(pdb_dir, "1abc")
    assert len(recs) == 2
    # Both models should have populated confidence
    assert all(r.confidence is not None for r in recs)
    # Higher ranking_score wins
    ranked = rank_and_assign(recs)
    top = [r for r in ranked if r.model_idx == 0][0]
    assert top.raw_scores["ranking_score"] == pytest.approx(0.91)


def test_proteinx_parser(fake_tree: Path):
    recs = parse_proteinx(fake_tree / "1abc", "1abc")
    assert len(recs) == 2
    rank_and_assign(recs)
    top = [r for r in recs if r.model_idx == 0][0]
    assert top.raw_scores["ranking_score"] == pytest.approx(0.71)


def test_chai_parser(fake_tree: Path):
    recs = parse_chai(fake_tree / "1abc", "1abc")
    assert len(recs) == 2
    rank_and_assign(recs)
    top = [r for r in recs if r.model_idx == 0][0]
    assert top.raw_scores["aggregate_score"] == pytest.approx(0.60)


def test_dynamicbind_parser(fake_tree: Path):
    recs = parse_dynamicbind(fake_tree / "1abc", "1abc")
    assert len(recs) == 2
    assert all(r.ligand_path is not None for r in recs)
    rank_and_assign(recs)
    top = [r for r in recs if r.model_idx == 0][0]
    assert top.raw_scores["confidence"] == pytest.approx(0.82)


def test_missing_producer_silently_skipped(tmp_path: Path):
    # Empty dir — parsers should return []
    (tmp_path / "1abc").mkdir()
    for fn in (parse_af3, parse_proteinx, parse_chai, parse_dynamicbind):
        assert fn(tmp_path / "1abc", "1abc") == []


def test_ranking_handles_none_confidence():
    recs = [
        ModelRecord(pdb_id="x", producer="af3",
                    structure_path=Path("/a"), confidence=None),
        ModelRecord(pdb_id="x", producer="af3",
                    structure_path=Path("/b"), confidence=0.5),
    ]
    rank_and_assign(recs)
    by_idx = {r.model_idx: r for r in recs}
    # The one with a real confidence wins; the None-confidence is last.
    assert by_idx[0].confidence == 0.5
    assert by_idx[1].confidence is None


def test_iter_pdb_dirs_skips_hidden(tmp_path: Path):
    (tmp_path / "1abc").mkdir()
    (tmp_path / "2def").mkdir()
    (tmp_path / ".hidden").mkdir()
    (tmp_path / "afile").write_text("x")
    names = [pid for pid, _ in iter_pdb_dirs(tmp_path)]
    assert names == ["1abc", "2def"]


# ---------------------------------------------------------------------------
# Sanitizer
# ---------------------------------------------------------------------------
gemmi = pytest.importorskip("gemmi")


def test_sanitize_pdb_roundtrip(tmp_path: Path):
    from plb_bench.sanitizer import sanitize
    p = tmp_path / "x.pdb"
    p.write_text(MINI_PDB)
    out = sanitize(p, pdb_id="test")
    assert out.n_chains >= 2           # protein chain A + ligand chain B
    assert "data_TEST" in out.text
    assert "_atom_site" in out.text
    # The sanitizer must produce entity metadata that the raw PDB lacked
    assert "_entity" in out.text


def test_sanitize_cif_roundtrip(tmp_path: Path):
    from plb_bench.sanitizer import sanitize
    p = tmp_path / "x.cif"
    p.write_text(MINI_CIF)
    out = sanitize(p, pdb_id="test")
    assert out.n_chains >= 1
    assert "_atom_site" in out.text


def test_sanitize_appends_sdf_ligand(tmp_path: Path):
    from plb_bench.sanitizer import sanitize
    pdb = tmp_path / "rec.pdb"
    pdb.write_text(MINI_PDB)
    sdf = tmp_path / "lig.sdf"
    sdf.write_text(_tiny_sdf())
    out = sanitize(pdb, pdb_id="db", ligand_path=sdf)
    # Expect an extra ligand chain beyond what was in the PDB
    assert out.n_chains >= 2
    # The merged CIF mentions the SDF-derived atoms
    assert "_atom_site" in out.text


def test_sanitize_in_memory_streaming(tmp_path: Path):
    """The sanitizer must not leak temp files to the output directory."""
    from plb_bench.sanitizer import sanitize
    p = tmp_path / "x.pdb"
    p.write_text(MINI_PDB)
    before = set(tmp_path.iterdir())
    out = sanitize(p, pdb_id="test")
    after = set(tmp_path.iterdir())
    assert before == after
    buf = out.as_buffer()
    assert buf.read().startswith("data_")


# ---------------------------------------------------------------------------
# References
# ---------------------------------------------------------------------------
def test_reference_local_hit(tmp_path: Path):
    from plb_bench.references import get_reference
    (tmp_path / "1ABC.cif").write_text("data_1ABC\n")
    path, source = get_reference("1abc", tmp_path, allow_download=False)
    assert source == "local"
    assert path.name == "1ABC.cif"


def test_reference_download_disabled_raises(tmp_path: Path):
    from plb_bench.references import get_reference
    with pytest.raises(FileNotFoundError):
        get_reference("9zzz", tmp_path, allow_download=False)
