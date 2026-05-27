"""
rescore_pla_dynamicbind.py
==========================
Standalone replacement for the three notebook cells that re-score
DynamicBind PLA and write evalspreadsheets/pla/ CSVs.

Equivalent to running cells:
  86ec0fa1  (Stage 2 — score with OST)
  3f56e1e6  (Stage 3 — pocket RMSD, incremental)
  0bb119b0  (Stage 4 — export CSVs)

Forces the repo's own plb_bench into sys.path so the kernel environment
never matters.

Run:
    conda activate plb
    python -u scoring/rescore_pla_dynamicbind.py 2>&1 | tee /tmp/pla_rescore.log
"""

from __future__ import annotations

import sys
import logging
from pathlib import Path

# ── Force our plb_bench to the front of sys.path ────────────────────────────
_BENCH_ROOT = Path("/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark")
_PKG_ROOT   = _BENCH_ROOT / "scoring" / "plb_bench"

# Remove any stale plb_bench already cached in sys.modules
for _k in list(sys.modules.keys()):
    if "plb_bench" in _k:
        del sys.modules[_k]

# Insert our package root at position 0 so it wins over any installed copy
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))
else:
    sys.path.remove(str(_PKG_ROOT))
    sys.path.insert(0, str(_PKG_ROOT))

# Verify we're using the right one
import plb_bench as _pb
_pb_file = getattr(_pb, "__file__", None)
assert _pb_file and str(_PKG_ROOT) in _pb_file, (
    f"Wrong plb_bench loaded from: {_pb_file}\n"
    f"Expected under: {_PKG_ROOT}"
)

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    stream=sys.stdout,
    force=True,
)
log = logging.getLogger(__name__)
log.info("plb_bench loaded from: %s", _pb_file)

# ── Paths (mirror plb_bench_run_pla.ipynb cell 4d3ee718) ────────────────────
DATA_ROOT  = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive/plb_bench_data")
REFS_DIR   = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive/references_ref_cifs")
OUTPUT_DIR = _BENCH_ROOT / "plb_bench_output" / "pla"
EVAL_DIR   = _BENCH_ROOT / "evalspreadsheets" / "pla"
POCKET_RADIUS    = 5.0
POCKET_CLEAN_THR = 20.0  # not used for PLA CSVs but kept for reference

PRODUCERS_TO_SCORE = ["dynamicbind_pla"]

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
EVAL_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# DynamicBind-correct score_model — monkey-patched into plb_bench.scoring
# ============================================================
# The shared plb_bench/scoring.py uses the legacy LigandScorer.score API which
# gives ~40% bisy_rmsd/lddt_pli coverage and binary QS scores for DynamicBind.
# We replace score_model here with an implementation that mirrors the canonical
# PLB pipeline (4_score/run_all_metrics.py) exactly:
#   • MMCIFPrep for loading  (same-entity ligands, avoids cross-entity bug)
#   • SCRMSDScorer / LDDTPLIScorer via assignment[0] + score_matrix
#   • QSEntity + QSScorer with DTP→LIG reference rename for chain pairing
#   • Bond wiring for DynamicBind LIG (SDF position-match first, CIF fallback)
# This patch is ONLY applied to this process; plb_bench/scoring.py on disk is
# unchanged so the main PLB bench pipeline for af3/chai/protenix/boltz is not
# affected.

import math as _math
import os as _os
import re as _re
import tempfile as _tempfile

import plb_bench.scoring as _scoring_mod


def _db_load_complex(cif_text: str):
    with _tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as f:
        f.write(cif_text)
        tmp = f.name
    try:
        try:
            from ost.mol.alg.scoring_base import MMCIFPrep
            ent, ligs = MMCIFPrep(tmp, extract_nonpoly=True)
            return ent, list(ligs)
        except (ImportError, Exception) as _e:
            log.debug("MMCIFPrep unavailable/failed (%s); falling back", _e)
        from ost import io as _ost_io
        result = _ost_io.LoadMMCIF(tmp, seqres=True, info=True, fault_tolerant=True)
        ent = result[0] if isinstance(result, tuple) else result
        ligs = _scoring_mod._get_ligands(ent)
        return ent, ligs
    finally:
        try:
            _os.unlink(tmp)
        except OSError:
            pass


_DB_IONS = _scoring_mod._IONS_AND_COFACTORS


def _db_filter_ligs(raw_ligs):
    out = []
    for r in raw_ligs:
        name = r.name if hasattr(r, "name") else ""
        if name in _DB_IONS:
            continue
        v = r.Select("ele != H")
        cnt = getattr(v, "atom_count", None)
        if cnt is None:
            try:
                cnt = v.GetAtomCount()
            except Exception:
                cnt = 0
        if cnt > 0:
            out.append(v)
    return out


def _db_safe_score(score_matrix, assignment):
    if not assignment:
        return None
    try:
        trg_i, mdl_i = assignment[0]
        v = score_matrix[trg_i, mdl_i]
        if v is None:
            return None
        f = float(v)
        return None if (_math.isnan(f) or _math.isinf(f)) else f
    except Exception:
        return None


def _db_wire_from_sdf(model_ent, sdf_path):
    try:
        from rdkit import Chem as _Chem
        import ost.mol as _omol
        suppl = _Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
        rdmol = next((m for m in suppl if m is not None and m.GetNumConformers() > 0), None)
        if rdmol is None:
            return 0
        try:
            _Chem.SanitizeMol(rdmol)
        except Exception:
            pass
        conf = rdmol.GetConformer(0)
        lig_residues = [r for ch in model_ent.chains for r in ch.residues if r.name.startswith("LIG")]
        if not lig_residues:
            return 0
        _BO = {
            _Chem.rdchem.BondType.SINGLE: 1, _Chem.rdchem.BondType.DOUBLE: 2,
            _Chem.rdchem.BondType.TRIPLE: 3, _Chem.rdchem.BondType.AROMATIC: 1,
        }
        total = 0
        for lig_res in lig_residues:
            pos_map = {(round(a.pos.x, 2), round(a.pos.y, 2), round(a.pos.z, 2)): a
                       for a in lig_res.atoms}
            idx_map = {}
            for i in range(rdmol.GetNumAtoms()):
                p = conf.GetAtomPosition(i)
                h = pos_map.get((round(p.x, 2), round(p.y, 2), round(p.z, 2)))
                if h is not None:
                    idx_map[i] = h
            if not idx_map:
                continue
            import ost.mol._ost_mol as _ost_mol_impl
            edi = model_ent.EditXCS(_ost_mol_impl.EditMode.UNBUFFERED_EDIT)
            wired = 0
            for bond in rdmol.GetBonds():
                a1 = idx_map.get(bond.GetBeginAtomIdx())
                a2 = idx_map.get(bond.GetEndAtomIdx())
                if a1 and a2:
                    try:
                        edi.Connect(a1, a2, _BO.get(bond.GetBondType(), 1))
                        wired += 1
                    except Exception:
                        pass
            del edi
            total += wired
        return total
    except Exception as _e:
        log.warning("_db_wire_from_sdf failed: %s", _e)
        return 0


def _db_wire_from_cif(model_ent, cif_text: str) -> int:
    """Wire _chem_comp_bond bonds from cif_text into model_ent in-place.

    Replicates _wire_lig_bonds_into_entity from plb_bench.scoring but uses
    EditXCS(0) instead of EditXCS(STANDARD_EDIT) to avoid AttributeError on
    this OST version where ost.mol.STANDARD_EDIT does not exist.
    """
    try:
        bonds = _scoring_mod._parse_cif_lig_bonds(cif_text)
        if not bonds:
            return 0
        lig_residues = [
            res for ch in model_ent.chains
            for res in ch.residues
            if res.name.startswith("LIG")
        ]
        if not lig_residues:
            return 0
        import ost.mol._ost_mol as _ost_mol_impl
        edi = model_ent.EditXCS(_ost_mol_impl.EditMode.UNBUFFERED_EDIT)
        wired = 0
        for lig_res in lig_residues:
            atom_map = {a.name: a for a in lig_res.atoms}
            for a1_name, a2_name, bo in bonds:
                a1 = atom_map.get(a1_name)
                a2 = atom_map.get(a2_name)
                if a1 is None or a2 is None:
                    continue
                try:
                    edi.Connect(a1, a2, bo)
                    wired += 1
                except Exception:
                    pass
        del edi
        log.debug("_db_wire_from_cif: wired %d bonds", wired)
        return wired
    except Exception as _e:
        log.warning("_db_wire_from_cif failed: %s", _e)
        return 0


def _db_scrmsd(mdl_ent, trg_ent, mdl_ligs_ent, trg_ligs, sdf_path=None):
    from ost.mol.alg.ligand_scoring_scrmsd import SCRMSDScorer
    if not trg_ligs:
        return None
    # Use entity handles with bonds wired via _db_wire_from_sdf / _db_wire_from_cif.
    # RDKit Mol is NOT used — this OST version rejects rdkit.Chem.Mol as model_ligands.
    if not mdl_ligs_ent:
        return None
    for sm in (True, False):
        try:
            sc = SCRMSDScorer(mdl_ent, trg_ent, mdl_ligs_ent, trg_ligs,
                              substructure_match=sm)
            v = _db_safe_score(sc.score_matrix, sc.assignment)
            if v is not None:
                log.info("BiSyRMSD=%.3f (sm=%s)", v, sm)
                return v
        except Exception as _e:
            log.info("BiSyRMSD sm=%s: %s", sm, _e)
    return None


def _db_lddt_pli(mdl_ent, trg_ent, mdl_ligs_ent, trg_ligs, sdf_path=None):
    from ost.mol.alg.ligand_scoring_lddtpli import LDDTPLIScorer
    if not trg_ligs:
        return None
    # Use entity handles with bonds wired via _db_wire_from_sdf / _db_wire_from_cif.
    # RDKit Mol is NOT used — this OST version rejects rdkit.Chem.Mol as model_ligands.
    if not mdl_ligs_ent:
        return None
    for sm in (True, False):
        try:
            sc = LDDTPLIScorer(mdl_ent, trg_ent, mdl_ligs_ent, trg_ligs,
                               substructure_match=sm)
            v = _db_safe_score(sc.score_matrix, sc.assignment)
            if v is not None:
                log.info("lDDT-PLI=%.3f (sm=%s)", v, sm)
                return v
        except Exception as _e:
            log.info("lDDT-PLI sm=%s: %s", sm, _e)
    return None


def _db_qs_global(mdl_ent, trg_ent, ref_cif_text=None):
    qs_ref = trg_ent
    model_has_lig = any(r.name == "LIG" for ch in mdl_ent.chains for r in ch.residues)
    if model_has_lig and ref_cif_text is not None:
        ref_lig = _scoring_mod._get_ref_lig_name(trg_ent)
        if ref_lig and ref_lig != "LIG":
            try:
                renamed = _re.sub(r'\b' + _re.escape(ref_lig) + r'\b', "LIG", ref_cif_text)
                from ost import io as _ost_io
                with _tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as f:
                    f.write(renamed)
                    tmp2 = f.name
                try:
                    result = _ost_io.LoadMMCIF(tmp2, seqres=True, info=True, fault_tolerant=True)
                    qs_ref = result[0] if isinstance(result, tuple) else result
                    log.info("QS: renamed ref %s→LIG", ref_lig)
                finally:
                    try:
                        _os.unlink(tmp2)
                    except OSError:
                        pass
            except Exception as _e:
                log.warning("QS ref rename failed (%s); using original ref", _e)
    try:
        from ost.mol.alg.qsscore import QSEntity, QSScorer as _QSScorer
        result = _QSScorer(QSEntity(qs_ref), QSEntity(mdl_ent)).Score()
        v = getattr(result, "qs_global", None)
        if v is not None:
            return float(v)
    except Exception as _e:
        log.info("QS approach 0 failed: %s", _e)
    try:
        from ost.mol.alg.scoring import Scorer
        v = Scorer(mdl_ent, qs_ref).qs_global
        if v is not None:
            return float(v)
    except Exception:
        pass
    try:
        from ost.mol.alg.qsscore import QSScorer
        from ost.mol.alg import chain_mapping as _cm
        mapper = _cm.ChainMapper(qs_ref, n_max_naive=1)
        mapping = mapper.GetMapping(mdl_ent)
        return float(QSScorer(mdl_ent, mapping.alns, qs_ref).global_score)
    except Exception:
        pass
    return None


def _db_score_model(sanitized_model, reference_cif_path, ligand_sdf_path=None):
    result = {"bisy_rmsd": None, "lddt_pli": None, "qs_global": None}
    try:
        from plb_bench.references import read_reference_text
        ref_text = read_reference_text(reference_cif_path)
    except Exception as _e:
        log.error("cannot read reference %s: %s", reference_cif_path, _e)
        return result
    model_text = sanitized_model.text
    try:
        mdl_ent, mdl_ligs_raw = _db_load_complex(model_text)
    except Exception as _e:
        log.error("cannot load model: %s", _e)
        return result
    try:
        trg_ent, trg_ligs_raw = _db_load_complex(ref_text)
    except Exception as _e:
        log.error("cannot load reference: %s", _e)
        return result
    # Wire bonds into entity handles (used as fallback candidate).
    # _db_wire_from_sdf uses position-matching against the SDF file.
    # _db_wire_from_cif parses _chem_comp_bond from the CIF text.
    # Both use EditXCS(0) to avoid the STANDARD_EDIT AttributeError in this OST version.
    bonds = 0
    if ligand_sdf_path and Path(ligand_sdf_path).exists():
        bonds = _db_wire_from_sdf(mdl_ent, ligand_sdf_path)
    if bonds == 0:
        bonds = _db_wire_from_cif(mdl_ent, model_text)
    # Filter entity handle ligands (Candidate B)
    mdl_ligs_ent = _db_filter_ligs(mdl_ligs_raw)
    trg_ligs = _db_filter_ligs(trg_ligs_raw)
    log.info("Ligands: model_ent=%d, target=%d", len(mdl_ligs_ent), len(trg_ligs))
    if trg_ligs:
        result["bisy_rmsd"] = _db_scrmsd(mdl_ent, trg_ent, mdl_ligs_ent, trg_ligs,
                                          sdf_path=ligand_sdf_path)
        result["lddt_pli"]  = _db_lddt_pli(mdl_ent, trg_ent, mdl_ligs_ent, trg_ligs,
                                             sdf_path=ligand_sdf_path)
    result["qs_global"] = _db_qs_global(mdl_ent, trg_ent, ref_cif_text=ref_text)
    return result


# Patch into plb_bench.scoring so score_tree._score_group uses it
_scoring_mod.score_model = _db_score_model
log.info("Patched plb_bench.scoring.score_model with DynamicBind-correct implementation")

# ============================================================
# Stage 2 — Score with OpenStructure
# ============================================================
log.info("=" * 60)
log.info("Stage 2: Score dynamicbind_pla with OST")

import multiprocessing as _multiprocessing          # noqa: E402
import concurrent.futures as _cf                    # noqa: E402

from plb_bench.score_tree import (                  # noqa: E402
    ScoreConfig, _iter_groups, _score_group, _write_output,
)


def _db_patched_score_group(task):
    """Fork-safe worker: re-applies the monkey-patch then calls _score_group.

    With fork the parent's sys.modules (including the patched score_model) is
    inherited, but any lazy `from plb_bench.scoring import score_model` inside
    _score_group would re-bind to whatever is in sys.modules at call time.
    Re-patching here guarantees the right implementation regardless of how
    score_tree resolves the reference.
    """
    import plb_bench.scoring as _sm
    _sm.score_model = _db_score_model
    return _score_group(task)


cfg = ScoreConfig(
    normalized_root       = DATA_ROOT,
    refs_dir              = REFS_DIR,
    output_dir            = OUTPUT_DIR,
    producers             = PRODUCERS_TO_SCORE,
    pdb_ids               = None,
    n_workers             = None,
    use_ray               = False,
    download_missing_refs = True,
    output_format         = "parquet",
)

tasks = [
    (prod, pid, str(pdir), str(cfg.refs_dir), cfg.download_missing_refs)
    for prod, pid, pdir in _iter_groups(cfg)
]
log.info("Collected %d (producer, pdb_id) tasks", len(tasks))

# Parallel scoring via fork — workers inherit the patched score_model.
# _db_patched_score_group re-applies the patch inside each worker to handle
# any lazy imports inside _score_group.
_n_workers = max(1, (_os.cpu_count() or 4) - 1)
log.info("Parallel scoring with %d workers (fork)", _n_workers)

_fork_ctx = _multiprocessing.get_context("fork")
rows: list[dict] = []
with _cf.ProcessPoolExecutor(max_workers=_n_workers, mp_context=_fork_ctx) as _ex:
    _futs = {_ex.submit(_db_patched_score_group, t): i for i, t in enumerate(tasks)}
    _done = 0
    for _fut in _cf.as_completed(_futs):
        try:
            rows.extend(_fut.result())
        except Exception as _e:
            log.warning("task failed: %s", _e)
        _done += 1
        if _done % 20 == 0:
            log.info("  completed %d / %d tasks ...", _done, len(tasks))

new_parquet = _write_output(rows, cfg)
log.info("Scoring complete -> %s", new_parquet)

import pandas as pd  # noqa: E402

new_df = pd.read_parquet(new_parquet)
log.info("Scored %d rows for %s", len(new_df), PRODUCERS_TO_SCORE)
log.info("bisy_rmsd coverage : %.1f%%", 100 * new_df["bisy_rmsd"].notna().mean())
log.info("lddt_pli  coverage : %.1f%%", 100 * new_df["lddt_pli"].notna().mean())
log.info("qs_global coverage : %.1f%%", 100 * new_df["qs_global"].notna().mean())

# ============================================================
# Stage 3 — Pocket RMSD (incremental merge)
# ============================================================
log.info("=" * 60)
log.info("Stage 3: Pocket RMSD (incremental)")

import gzip as _gzip
import os as _os
import tempfile as _tmp
import numpy as _np  # noqa: E402

_SOLVENT_EXT = {
    "HOH","DOD","WAT","EDO","GOL","PEG","PO4","SO4","ACT",
    "MES","TRIS","BME","DTT","MPD","FMT","ACE","NH2",
    "MG","ZN","CA","NA","CL","K","MN","FE","CU","CO",
}


def _ost_load_cif(path: str):
    from ost import io as _ost_io
    p = Path(path)
    text = _gzip.open(p, "rt").read() if p.suffix == ".gz" else p.read_text(errors="replace")
    with _tmp.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as f:
        f.write(text); tmp = f.name
    try:
        result = _ost_io.LoadMMCIF(tmp, fault_tolerant=True)
        return result[0] if isinstance(result, tuple) else result
    finally:
        try: _os.unlink(tmp)
        except OSError: pass


def _poly_types() -> set:
    from ost import mol
    types: set = set()
    for name in ("CHAINTYPE_POLY_PEPTIDE_L","CHAINTYPE_POLY_PEPTIDE_D",
                  "CHAINTYPE_POLY_DN","CHAINTYPE_POLY_RN","CHAINTYPE_POLY",
                  "CHAINTYPE_POLY_SAC_D","CHAINTYPE_POLY_SAC_L","CHAINTYPE_WATER"):
        ct = getattr(mol, name, None)
        if ct is not None: types.add(ct)
    return types


def _collect_ca_by_chain(ent) -> dict:
    poly = _poly_types()
    result: dict = {}
    for chain in ent.chains:
        if chain.chain_type not in poly: continue
        residues = []
        for res in chain.residues:
            ca = res.FindAtom("CA")
            if ca.IsValid():
                residues.append((res.number.num, _np.array([ca.pos.x, ca.pos.y, ca.pos.z])))
        if residues:
            result[chain.name] = residues
    return result


def _get_all_lig_instances(ent, lig_name=None) -> list:
    poly = _poly_types()
    named, others = [], []
    for chain in ent.chains:
        if chain.chain_type in poly: continue
        for res in chain.residues:
            if res.name in _SOLVENT_EXT: continue
            pts = _np.array([[at.pos.x, at.pos.y, at.pos.z] for at in res.atoms])
            if pts.shape[0] == 0: continue
            if lig_name and res.name.upper() == lig_name.upper():
                named.append(pts)
            else:
                others.append(pts)
    if named: return named
    if others: return [max(others, key=lambda x: x.shape[0])]
    return []


def _pick_best_ref_lig(ref_copies: list, anchor: _np.ndarray) -> _np.ndarray:
    if not ref_copies: return _np.empty((0, 3))
    if len(ref_copies) == 1: return ref_copies[0]
    return min(ref_copies, key=lambda pts: _np.linalg.norm(pts.mean(axis=0) - anchor))


def _match_ca(ref_by_chain: dict, mdl_by_chain: dict):
    ref_flat = {(c, r): p for c, res in ref_by_chain.items() for r, p in res}
    mdl_flat = {(c, r): p for c, res in mdl_by_chain.items() for r, p in res}
    common = sorted(set(ref_flat) & set(mdl_flat))
    if len(common) >= 10:
        return _np.array([ref_flat[k] for k in common]), _np.array([mdl_flat[k] for k in common])
    ref_chains = sorted(ref_by_chain.items(), key=lambda x: len(x[1]), reverse=True)
    mdl_chains = sorted(mdl_by_chain.items(), key=lambda x: len(x[1]), reverse=True)
    rp, mp = [], []
    for (rcid, rres), (mcid, mres) in zip(ref_chains, mdl_chains):
        rdict, mdict = dict(rres), dict(mres)
        common_resnums = sorted(set(rdict) & set(mdict))
        if len(common_resnums) >= 5:
            for rnum in common_resnums:
                rp.append(rdict[rnum]); mp.append(mdict[rnum])
        else:
            n = min(len(rres), len(mres))
            for i in range(n):
                rp.append(sorted(rres)[i][1]); mp.append(sorted(mres)[i][1])
    if len(rp) < 3: return None, None
    return _np.array(rp), _np.array(mp)


def _kabsch(ref_pts, mdl_pts):
    rc, mc = ref_pts.mean(0), mdl_pts.mean(0)
    H = (mdl_pts - mc).T @ (ref_pts - rc)
    U, _, Vt = _np.linalg.svd(H)
    d = _np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ _np.diag([1., 1., d]) @ U.T
    return R, rc - R @ mc


def compute_pocket_rmsd(mdl_cif_path: str, ref_cif_path: str,
                         pdb_id: str = "", radius: float = POCKET_RADIUS):
    try:
        mdl = _ost_load_cif(mdl_cif_path)
        ref = _ost_load_cif(ref_cif_path)
        parts = str(pdb_id).split("_")
        lig_name = parts[-1].upper() if len(parts) >= 2 else None
        ref_copies = _get_all_lig_instances(ref, lig_name)
        if not ref_copies: return None
        mdl_copies = _get_all_lig_instances(mdl, lig_name)
        if not mdl_copies: return None
        ref_by_chain = _collect_ca_by_chain(ref)
        mdl_by_chain = _collect_ca_by_chain(mdl)
        ref_pts, mdl_pts = _match_ca(ref_by_chain, mdl_by_chain)
        if ref_pts is None or len(ref_pts) < 3: return None
        R, t = _kabsch(ref_pts, mdl_pts)
        ref_ca_centroid = ref_pts.mean(axis=0)
        ref_lig = _pick_best_ref_lig(ref_copies, ref_ca_centroid)
        if ref_lig.shape[0] == 0: return None
        mask = _np.array([_np.linalg.norm(ref_lig - pos, axis=1).min() <= radius for pos in ref_pts])
        if mask.sum() < 3: return None
        pocket_ref = ref_pts[mask]
        all_mdl_t = (R @ mdl_pts.T).T + t
        pocket_mdl_t = _np.array([
            all_mdl_t[_np.argmin(_np.linalg.norm(all_mdl_t - ref_pos, axis=1))]
            for ref_pos in pocket_ref
        ])
        return float(_np.sqrt(_np.mean(_np.sum((pocket_ref - pocket_mdl_t) ** 2, axis=1))))
    except Exception:
        return None


# Load existing pocket RMSD parquet to cache values for all other producers
prmsd_path = OUTPUT_DIR / "benchmark_with_pocket_rmsd.parquet"
pocket_cache: dict = {}
other_rows = pd.DataFrame()

if prmsd_path.exists():
    existing = pd.read_parquet(prmsd_path)
    # Keep all non-dynamicbind_pla rows intact (fresh from their own scoring parquet)
    other_rows = existing[~existing["producer"].isin(PRODUCERS_TO_SCORE)].copy()
    # Build pocket_rmsd lookup for dynamicbind_pla rows from the cache
    db_cached = existing[existing["producer"].isin(PRODUCERS_TO_SCORE)]
    pocket_cache = (
        db_cached.set_index(["producer", "pdb_id", "model_idx"])["pocket_rmsd"]
        .dropna().to_dict()
    )
    log.info("Loaded cached pocket_rmsd for %d existing dynamicbind_pla rows",
             len(pocket_cache))

# Apply pocket RMSD: use cache where available, compute where missing
pocket_vals = []
n_cached = n_computed = n_skip = n_miss_ref = n_miss_mdl = 0

for i, (_, row) in enumerate(new_df.iterrows()):
    if i % 500 == 0:
        log.info("  pocket RMSD %d / %d ...", i, len(new_df))

    if row["status"] != "ok":
        pocket_vals.append(None); n_skip += 1; continue

    key = (row["producer"], row["pdb_id"], int(row["model_idx"]))
    if key in pocket_cache:
        pocket_vals.append(pocket_cache[key]); n_cached += 1; continue

    pdb_id    = row["pdb_id"]
    model_idx = int(row["model_idx"])

    # Find reference CIF
    ref_cif = None
    ref_src = row.get("reference_source", None)
    if ref_src and Path(str(ref_src)).exists():
        ref_cif = Path(str(ref_src))
    if ref_cif is None:
        for ext in (f"{pdb_id}.cif", f"{pdb_id}.cif.gz"):
            c = REFS_DIR / ext
            if c.exists(): ref_cif = c; break
    if ref_cif is None:
        pocket_vals.append(None); n_miss_ref += 1; continue

    mdl_cif = DATA_ROOT / row["producer"] / pdb_id / f"model_{model_idx:03d}.cif"
    if not mdl_cif.exists():
        pocket_vals.append(None); n_miss_mdl += 1; continue

    pocket_vals.append(compute_pocket_rmsd(str(mdl_cif), str(ref_cif), pdb_id=pdb_id))
    n_computed += 1

new_df["pocket_rmsd"] = pocket_vals
log.info("Pocket RMSD: %d from cache, %d computed, %d skipped, "
         "%d missing ref, %d missing mdl",
         n_cached, n_computed, n_skip, n_miss_ref, n_miss_mdl)
log.info("pocket_rmsd coverage (dynamicbind_pla): %.1f%%",
         100 * new_df["pocket_rmsd"].notna().mean())

# Merge with the other producers and save
merged = pd.concat([other_rows, new_df], ignore_index=True)
merged.to_parquet(prmsd_path, index=False)
log.info("Saved benchmark_with_pocket_rmsd.parquet: %d rows -> %s", len(merged), prmsd_path)

# ============================================================
# Stage 4 — Export CSVs (mirror notebook cell 0bb119b0)
# ============================================================
log.info("=" * 60)
log.info("Stage 4: Export CSVs -> %s", EVAL_DIR)

_PLA_RENAME = {
    "bisy_rmsd":   "pose rmsd",
    "pocket_rmsd": "pocket rmsd",
    "lddt_pli":    "lddt-pli",
    "qs_global":   "qs score",
}
_PLA_METRIC_COLS = ["bisy_rmsd", "pocket_rmsd", "qs_global", "lddt_pli"]
# Column order matches af3_plamainligand_*.csv exactly.
_PLA_OUT_COLS = ["id", "pose rmsd", "pocket rmsd", "qs score", "lddt-pli", "confidence"]


# Reload the full merged parquet — only export dynamicbind_pla CSVs.
# Other producers (af3, chai, protenix, boltz) have correct CSVs from the
# main PLB bench pipeline and must NOT be overwritten here.
full = pd.read_parquet(prmsd_path)
pla_ok = full[(full["status"] == "ok") & (full["producer"] == "dynamicbind_pla")].copy()
# FULL compound key ("uniprot_ligand", e.g. "b6ywb8_dtp"), matching af3/chai —
# do NOT truncate to the UniProt prefix (that collapses multiple ligands per
# protein and breaks cross-producer joins).
pla_ok["id"] = pla_ok["pdb_id"].str.lower()
avail = [c for c in _PLA_METRIC_COLS if c in pla_ok.columns]


def _pla_make_best(grp):
    scored = grp.dropna(subset=["bisy_rmsd"])
    if len(scored) > 0:
        best_idx = scored.groupby("id")["bisy_rmsd"].idxmin().dropna().astype(int)
        best_scored = grp.loc[best_idx].copy()
    else:
        best_scored = grp.iloc[0:0].copy()
    unscored_ids = set(grp["id"].unique()) - set(best_scored["id"])
    if unscored_ids:
        fb = grp[grp["id"].isin(unscored_ids)].groupby("id", as_index=False).first()
        return pd.concat([best_scored, fb], ignore_index=True)
    return best_scored


def _pla_make_avg(grp):
    return grp.groupby("id")[avail + ["confidence"]].mean().reset_index()


def _pla_finalise(frame):
    frame = frame.rename(columns=_PLA_RENAME)
    cols = [c for c in _PLA_OUT_COLS if c in frame.columns]
    return frame[cols].sort_values("id").reset_index(drop=True)


for producer, grp in pla_ok.groupby("producer"):
    best = _pla_finalise(_pla_make_best(grp))
    best.to_csv(EVAL_DIR / f"{producer}mainligand_best.csv", index=False)
    log.info("best        %s: %d rows", producer, len(best))

    avg = _pla_finalise(_pla_make_avg(grp))
    avg.to_csv(EVAL_DIR / f"{producer}mainligand_avg.csv", index=False)
    log.info("avg         %s: %d rows", producer, len(avg))

    # clean variants: models with pocket_rmsd <= POCKET_CLEAN_THR (20 A)
    grp_clean = grp[grp["pocket_rmsd"].notna()
                    & (grp["pocket_rmsd"] <= POCKET_CLEAN_THR)]
    clean_best = (_pla_finalise(_pla_make_best(grp_clean)) if len(grp_clean)
                  else pd.DataFrame(columns=_PLA_OUT_COLS))
    clean_best.to_csv(EVAL_DIR / f"{producer}mainligand_clean_best.csv", index=False)
    log.info("clean_best  %s: %d rows", producer, len(clean_best))

    clean_avg = (_pla_finalise(_pla_make_avg(grp_clean)) if len(grp_clean)
                 else pd.DataFrame(columns=_PLA_OUT_COLS))
    clean_avg.to_csv(EVAL_DIR / f"{producer}mainligand_clean_avg.csv", index=False)
    log.info("clean_avg   %s: %d rows", producer, len(clean_avg))

log.info("All CSVs written to: %s", EVAL_DIR)

# ── Quick summary ─────────────────────────────────────────────────────────────
log.info("=" * 60)
db = pla_ok[pla_ok["producer"] == "dynamicbind_pla"]
log.info("dynamicbind_pla summary: %d compounds, %d models",
         db["id"].nunique(), len(db))
log.info("  bisy_rmsd  coverage: %.1f%%", 100 * db["bisy_rmsd"].notna().mean())
log.info("  lddt_pli   coverage: %.1f%%", 100 * db["lddt_pli"].notna().mean())
log.info("  pocket_rmsd coverage: %.1f%%", 100 * db["pocket_rmsd"].notna().mean())
log.info("Done.")
