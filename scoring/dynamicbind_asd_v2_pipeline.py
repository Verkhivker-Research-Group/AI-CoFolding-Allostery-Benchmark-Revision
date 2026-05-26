"""
dynamicbind_asd_v2_pipeline.py
==============================
Full end-to-end pipeline for DynamicBind ASD v2 runs.

Data source: AI-CoFolding-Archive/asd_dynamicbind_results_v2/asd_dynamicbind_results/
Folder convention: asd_{PDB}_{LIG}/index0_idx_0/rank{N}_*_relaxed.{pdb,sdf}
                   e.g.  asd_3HRF_P47  ->  compound_key = 3hrf_p47

Pipeline stages
---------------
1. Discovery + Normalization
   Parses  asd_{PDB}_{LIG}  folders.
   compound_key = {pdb}_{lig}.lower()  (e.g. 3hrf_p47)
   Sanitizes receptor PDB + ligand SDF into a unified CIF.
   Writes:  plb_bench_data/dynamicbind_asd_v2/{compound_key}/model_NNN.cif
                                                               manifest.json

2. Reference pre-warm
   Extracts 4-char PDB code from compound_key (before first '_').
   Downloads {PDB}-assembly1.cif.gz from RCSB and saves it as
   {compound_key}.cif.gz  in REFS_DIR so get_reference() finds it by
   standard name lookup.  If a local {PDB}*.cif.gz already exists, it is
   copied instead of re-downloaded (avoids duplicate network requests for
   compounds sharing the same PDB).

3. Scoring with OpenStructure
   Uses plb_bench.score_tree.run() with all fixes:
     - SDF-based ligand loading for proper bond connectivity
       (bisy_rmsd / lddt-pli work for LIG residue names)
     - QSEntity + QSScorer Approach 0 (continuous QS including ligand chain)
   Output: plb_bench_output/asd_v2/benchmark.parquet

4. Pocket RMSD (Bug-A + Bug-B fixed version)
   Bug A  — reference ligand copy selected by proximity to ref Cα centroid
            (correct for multi-subunit assemblies).
   Bug B  — geometric nearest-neighbour matching after Kabsch superposition,
            fully independent of residue-number offsets (e.g. Boltz starts at 1).
   Verbatim algorithm from plb_bench_run_pla.ipynb cell 83ecd346.

5. CSV export  (evalspreadsheets/main/)
   dynamicbind_asd_v2alloligand_best.csv       — min pose RMSD per compound
   dynamicbind_asd_v2alloligand_avg.csv        — mean over all models
   dynamicbind_asd_v2alloligand_clean_best.csv — best with pocket_rmsd ≤ 20 Å
   dynamicbind_asd_v2alloligand_clean_avg.csv  — mean with pocket_rmsd ≤ 20 Å

   Column order: id, pose rmsd, pocket rmsd, qs score, lddt-pli, confidence

Run from WSL with the plb conda env active:
    PYTHONPATH=/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark/scoring/plb_bench:$PYTHONPATH \\
    python -u scoring/dynamicbind_asd_v2_pipeline.py 2>&1 | tee /mnt/c/Temp/dynamicbind_asd_v2.log
"""

from __future__ import annotations

import csv as _csv
import gzip as _gzip
import logging
import os as _os
import re
import shutil
import sys
import tempfile as _tmp
from collections import defaultdict
from pathlib import Path

import numpy as _np
import pandas as pd

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    stream=sys.stdout,
    force=True,
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BENCH_ROOT = Path("/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark")
ARCHIVE    = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive")

RAW_ROOT   = ARCHIVE / "asd_dynamicbind_results_v2" / "asd_dynamicbind_results"
DATA_ROOT  = ARCHIVE / "plb_bench_data"
REFS_DIR   = ARCHIVE / "references_ref_cifs"
OUTPUT_DIR = BENCH_ROOT / "plb_bench_output" / "asd_v2"
EVAL_DIR   = BENCH_ROOT / "evalspreadsheets" / "main"

PRODUCER         = "dynamicbind_asd_v2"
POCKET_RADIUS    = 5.0   # Å — matching ASD + PLA pipelines
POCKET_CLEAN_THR = 20.0  # Å — models above this omitted from clean CSVs

# ---------------------------------------------------------------------------
# plb_bench bootstrap
# ---------------------------------------------------------------------------
_pkg_root = BENCH_ROOT / "scoring" / "plb_bench"
if _pkg_root.exists() and str(_pkg_root) not in sys.path:
    sys.path.insert(0, str(_pkg_root))

from plb_bench.normalize import register_discoverer, normalize_producer  # noqa: E402
from plb_bench.schema import ModelRecord                                   # noqa: E402
from plb_bench.score_tree import ScoreConfig, run as score_run             # noqa: E402
from plb_bench.references import _search_local                             # noqa: E402


# ===========================================================================
# Stage 1 — Discovery + Normalization
# ===========================================================================

@register_discoverer(PRODUCER)
def _discover_dynamicbind_asd_v2(raw_root: Path):
    """
    Discover DynamicBind ASD v2 runs.

    Folder convention:  asd_{PDB}_{LIG}/index0_idx_0/rank{N}_*_relaxed.pdb
    compound_key      = {pdb}_{lig}.lower()   e.g. 3hrf_p47

    Uses *_relaxed.pdb exclusively to avoid counting the non-relaxed / relaxed
    pair as two models (matching dynamicbind_asd discoverer behaviour).
    Relaxed SDF is preferred when available; otherwise the non-relaxed SDF for
    the same rank is used as the source_ligand_path for SDF-based scoring.

    Confidence is taken from complete_affinity_prediction.csv column 'lddt'
    (DynamicBind's structural quality score, higher = better).
    Falls back to 'affinity', then to negative rank.
    """
    all_records: dict[str, list[ModelRecord]] = defaultdict(list)

    for folder in sorted(raw_root.iterdir()):
        if not folder.is_dir() or folder.name.startswith("."):
            continue

        # --- parse folder name: asd_{PDB}_{LIG} ----------------------------
        parts = folder.name.split("_")
        if len(parts) < 3 or parts[0].lower() != "asd":
            continue
        pdb_code     = parts[1]               # e.g. 3HRF (4 chars)
        lig_code     = "_".join(parts[2:])    # e.g. P47  (handles multi-part codes)
        compound_key = f"{pdb_code}_{lig_code}".lower()

        # --- per-rank confidence from complete_affinity_prediction.csv ------
        conf_map: dict[int, dict] = {}
        csv_path = folder / "complete_affinity_prediction.csv"
        if csv_path.exists():
            try:
                with csv_path.open() as fh:
                    for row in _csv.DictReader(fh):
                        rank_str = row.get("rank")
                        if rank_str is None:
                            continue
                        try:
                            rank = int(rank_str)
                        except ValueError:
                            continue
                        bucket = conf_map.setdefault(rank, {})
                        for k, v in row.items():
                            if k in ("name", "rank"):
                                continue
                            try:
                                bucket[k] = float(v)
                            except (TypeError, ValueError):
                                pass
            except Exception as exc:
                log.warning("[%s] csv parse error: %s", compound_key, exc)

        # --- collect relaxed receptor PDBs + paired ligand SDFs -------------
        for pdb_file in sorted(folder.rglob("rank*_receptor*_relaxed.pdb")):
            m = re.search(r"rank(\d+)", pdb_file.name)
            if not m:
                continue
            rank = int(m.group(1))

            # Prefer relaxed SDF; fall back to non-relaxed for same rank
            relaxed_sdfs = list(pdb_file.parent.glob(f"rank{rank}_ligand*_relaxed.sdf"))
            any_sdfs     = list(pdb_file.parent.glob(f"rank{rank}_ligand*.sdf"))
            sdf = (relaxed_sdfs[0] if relaxed_sdfs
                   else any_sdfs[0] if any_sdfs else None)

            raw        = conf_map.get(rank, {})
            confidence = raw.get("lddt", raw.get("affinity", -float(rank)))

            all_records[compound_key].append(ModelRecord(
                pdb_id=compound_key, producer=PRODUCER,
                structure_path=pdb_file,
                ligand_path=sdf,
                raw_scores=raw,
                confidence=confidence,
            ))

    for compound_key, records in sorted(all_records.items()):
        if records:
            yield compound_key, records


def run_normalization() -> list:
    log.info("=" * 60)
    log.info("Stage 1: Normalize")
    log.info("  raw_root : %s", RAW_ROOT)
    log.info("  data_root: %s", DATA_ROOT / PRODUCER)

    if not RAW_ROOT.is_dir():
        log.error("RAW_ROOT not found: %s", RAW_ROOT)
        sys.exit(1)

    entries = normalize_producer(PRODUCER, RAW_ROOT, DATA_ROOT, overwrite=True)
    n_ok   = sum(1 for e in entries if not e.failed)
    n_fail = sum(e.failed for e in entries)
    n_mdl  = sum(e.n_models for e in entries)
    log.info("Normalized %d compound groups — %d ok, %d failed — %d models total",
             len(entries), n_ok, n_fail, n_mdl)
    return entries


# ===========================================================================
# Stage 1b — Reference pre-warm
# ===========================================================================

def prewarm_refs() -> None:
    """
    For each normalized compound_key (e.g. '3hrf_p47'), extract the 4-char
    PDB code ('3hrf') and ensure REFS_DIR/{compound_key}.cif.gz exists so that
    score_tree's get_reference(pdb_id=compound_key) finds it by standard search.

    Download order:
      1. If {pdb_code}*.cif.gz already in REFS_DIR: copy → {compound_key}.cif.gz
      2. Else: download from RCSB (assembly1 then plain), save as {compound_key}.cif.gz
    """
    import requests as _req

    log.info("=" * 60)
    log.info("Stage 1b: Pre-warm references -> %s", REFS_DIR)
    REFS_DIR.mkdir(parents=True, exist_ok=True)

    prod_dir = DATA_ROOT / PRODUCER
    if not prod_dir.is_dir():
        log.warning("Normalized tree not found at %s — skipping ref pre-warm", prod_dir)
        return

    compound_dirs = sorted(
        p for p in prod_dir.iterdir() if p.is_dir() and not p.name.startswith(".")
    )

    skipped = ok_local = ok_rcsb = failed = 0

    for pdb_dir in compound_dirs:
        compound_key = pdb_dir.name
        dest = REFS_DIR / f"{compound_key}.cif.gz"

        if dest.exists():
            skipped += 1
            continue

        pdb_code = compound_key.split("_")[0].upper()  # e.g. 3HRF

        # --- 1. look for local copy by PDB code -----------------------------
        local = _search_local(REFS_DIR, pdb_code)
        if local is not None:
            shutil.copy2(str(local), str(dest))
            ok_local += 1
            log.debug("  local copy  %s  ->  %s", local.name, dest.name)
            continue

        # --- 2. download from RCSB ------------------------------------------
        downloaded = False
        for url in [
            f"https://files.rcsb.org/download/{pdb_code}-assembly1.cif.gz",
            f"https://files.rcsb.org/download/{pdb_code}.cif.gz",
        ]:
            try:
                r = _req.get(url, timeout=60)
                if r.status_code == 200:
                    dest.write_bytes(r.content)
                    ok_rcsb += 1
                    log.info("  downloaded  %s  (%s)", compound_key, pdb_code)
                    downloaded = True
                    break
            except Exception as exc:
                log.warning("  download failed  %s: %s", url, exc)

        if not downloaded:
            failed += 1
            log.warning("  FAILED  %s  (pdb=%s)", compound_key, pdb_code)

    log.info(
        "Pre-warm complete: %d from local cache, %d from RCSB, %d skipped, %d failed",
        ok_local, ok_rcsb, skipped, failed,
    )


# ===========================================================================
# Stage 2 — Score with OpenStructure (all fixes applied via scoring.py)
# ===========================================================================

def run_scoring() -> Path:
    log.info("=" * 60)
    log.info("Stage 2: Score with OST  ->  %s", OUTPUT_DIR)

    cfg = ScoreConfig(
        normalized_root       = DATA_ROOT,
        refs_dir              = REFS_DIR,
        output_dir            = OUTPUT_DIR,
        producers             = [PRODUCER],
        pdb_ids               = None,           # all compounds
        n_workers             = None,           # os.cpu_count()
        use_ray               = False,
        download_missing_refs = False,          # pre-warmed in stage 1b
        output_format         = "parquet",
    )
    out = score_run(cfg)
    log.info("Scoring complete  ->  %s", out)
    return out


# ===========================================================================
# Stage 3 — Pocket RMSD
# Bug-A + Bug-B fixed version (plb_bench_run_pla.ipynb cell 83ecd346)
# ===========================================================================

_SOLVENT_EXT = {
    "HOH", "DOD", "WAT", "EDO", "GOL", "PEG", "PO4", "SO4", "ACT",
    "MES", "TRIS", "BME", "DTT", "MPD", "FMT", "ACE", "NH2",
    "MG", "ZN", "CA", "NA", "CL", "K", "MN", "FE", "CU", "CO",
}


def _ost_load_cif(path: str):
    from ost import io as _ost_io
    p = Path(path)
    text = _gzip.open(p, "rt").read() if p.suffix == ".gz" else p.read_text(errors="replace")
    with _tmp.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as f:
        f.write(text)
        tmp = f.name
    try:
        result = _ost_io.LoadMMCIF(tmp, fault_tolerant=True)
        return result[0] if isinstance(result, tuple) else result
    finally:
        try:
            _os.unlink(tmp)
        except OSError:
            pass


def _poly_types() -> set:
    from ost import mol
    types: set = set()
    for name in (
        "CHAINTYPE_POLY_PEPTIDE_L", "CHAINTYPE_POLY_PEPTIDE_D",
        "CHAINTYPE_POLY_DN",        "CHAINTYPE_POLY_RN",
        "CHAINTYPE_POLY",           "CHAINTYPE_POLY_SAC_D",
        "CHAINTYPE_POLY_SAC_L",     "CHAINTYPE_WATER",
    ):
        ct = getattr(mol, name, None)
        if ct is not None:
            types.add(ct)
    return types


def _collect_ca_by_chain(ent) -> dict:
    poly = _poly_types()
    result: dict = {}
    for chain in ent.chains:
        if chain.chain_type not in poly:
            continue
        residues = []
        for res in chain.residues:
            ca = res.FindAtom("CA")
            if ca.IsValid():
                residues.append((res.number.num,
                                 _np.array([ca.pos.x, ca.pos.y, ca.pos.z])))
        if residues:
            result[chain.name] = residues
    return result


def _get_all_lig_instances(ent, lig_name: str | None = None) -> list:
    """Return all non-solvent ligand copies as xyz arrays.
    Named copies (matching lig_name) are preferred; otherwise the single
    largest non-solvent residue is returned (Bug-A fix)."""
    poly = _poly_types()
    named, others = [], []
    for chain in ent.chains:
        if chain.chain_type in poly:
            continue
        for res in chain.residues:
            if res.name in _SOLVENT_EXT:
                continue
            pts = _np.array([[at.pos.x, at.pos.y, at.pos.z] for at in res.atoms])
            if pts.shape[0] == 0:
                continue
            if lig_name and res.name.upper() == lig_name.upper():
                named.append(pts)
            else:
                others.append(pts)
    if named:
        return named
    if others:
        return [max(others, key=lambda x: x.shape[0])]
    return []


def _pick_best_ref_lig(ref_copies: list, anchor: _np.ndarray) -> _np.ndarray:
    """Bug-A fix: select the reference ligand copy closest to a reference-space anchor."""
    if not ref_copies:
        return _np.empty((0, 3))
    if len(ref_copies) == 1:
        return ref_copies[0]
    return min(ref_copies,
               key=lambda pts: _np.linalg.norm(pts.mean(axis=0) - anchor))


def _match_ca(ref_by_chain: dict, mdl_by_chain: dict):
    ref_flat = {(cid, rnum): pos
                for cid, residues in ref_by_chain.items()
                for rnum, pos in residues}
    mdl_flat = {(cid, rnum): pos
                for cid, residues in mdl_by_chain.items()
                for rnum, pos in residues}
    # Strategy A: exact (chain_name, resnum) match
    common = sorted(set(ref_flat) & set(mdl_flat))
    if len(common) >= 10:
        return (_np.array([ref_flat[k] for k in common]),
                _np.array([mdl_flat[k] for k in common]))
    # Strategy B: pair chains by descending size, match by resnum within pairs
    ref_chains = sorted(ref_by_chain.items(), key=lambda x: len(x[1]), reverse=True)
    mdl_chains = sorted(mdl_by_chain.items(), key=lambda x: len(x[1]), reverse=True)
    rp, mp = [], []
    for (rcid, rres), (mcid, mres) in zip(ref_chains, mdl_chains):
        rdict = dict(rres)
        mdict = dict(mres)
        common_resnums = sorted(set(rdict) & set(mdict))
        if len(common_resnums) >= 5:
            for rnum in common_resnums:
                rp.append(rdict[rnum])
                mp.append(mdict[rnum])
        else:
            n = min(len(rres), len(mres))
            for i in range(n):
                rp.append(sorted(rres)[i][1])
                mp.append(sorted(mres)[i][1])
    if len(rp) < 3:
        return None, None
    return _np.array(rp), _np.array(mp)


def _kabsch(ref_pts: _np.ndarray, mdl_pts: _np.ndarray):
    rc, mc = ref_pts.mean(0), mdl_pts.mean(0)
    H = (mdl_pts - mc).T @ (ref_pts - rc)
    U, _, Vt = _np.linalg.svd(H)
    d = _np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ _np.diag([1.0, 1.0, d]) @ U.T
    return R, rc - R @ mc


def compute_pocket_rmsd(mdl_cif_path: str, ref_cif_path: str,
                        pdb_id: str = "", radius: float = POCKET_RADIUS) -> float | None:
    """
    Pocket RMSD — Bug-A + Bug-B fixed version.

    Bug A: reference ligand copy selected by proximity to reference Cα centroid
           (correct for multi-subunit assemblies with multiple ligand copies).
    Bug B: after Kabsch, transform ALL model Cα atoms into reference space and
           find the nearest model Cα for each reference pocket Cα geometrically.
           This is completely independent of residue-number offsets.
    """
    try:
        mdl = _ost_load_cif(mdl_cif_path)
        ref = _ost_load_cif(ref_cif_path)

        # Ligand name = last element of compound_key (e.g. 'P47' from '3hrf_p47')
        parts = str(pdb_id).split("_")
        lig_name = parts[-1].upper() if len(parts) >= 2 else None

        # Step 1: collect all ligand copies
        ref_copies = _get_all_lig_instances(ref, lig_name)
        if not ref_copies:
            return None
        mdl_copies = _get_all_lig_instances(mdl, lig_name)
        if not mdl_copies:
            return None

        # Step 2: Cα matching + Kabsch superposition
        ref_by_chain = _collect_ca_by_chain(ref)
        mdl_by_chain = _collect_ca_by_chain(mdl)
        ref_pts, mdl_pts = _match_ca(ref_by_chain, mdl_by_chain)
        if ref_pts is None or len(ref_pts) < 3:
            return None
        R, t = _kabsch(ref_pts, mdl_pts)

        # Step 3: Bug-A — pick the ref ligand copy closest to ref Cα centroid
        ref_ca_centroid = ref_pts.mean(axis=0)
        ref_lig = _pick_best_ref_lig(ref_copies, ref_ca_centroid)
        if ref_lig.shape[0] == 0:
            return None

        # Step 4: define reference pocket (Cα within radius of ref ligand)
        mask = _np.array([
            _np.linalg.norm(ref_lig - pos, axis=1).min() <= radius
            for pos in ref_pts
        ])
        if mask.sum() < 3:
            return None
        pocket_ref = ref_pts[mask]

        # Step 5: transform ALL model Cα atoms into reference space
        all_mdl_t = (R @ mdl_pts.T).T + t  # shape (N_mdl, 3)

        # Step 6: Bug-B — geometric nearest-neighbour, independent of resnum
        pocket_mdl_t = _np.array([
            all_mdl_t[_np.argmin(_np.linalg.norm(all_mdl_t - ref_pos, axis=1))]
            for ref_pos in pocket_ref
        ])

        # Step 7: pocket RMSD
        return float(_np.sqrt(_np.mean(_np.sum((pocket_ref - pocket_mdl_t) ** 2, axis=1))))

    except Exception:
        return None


def run_pocket_rmsd(bench_parquet: Path) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("Stage 3: Pocket RMSD")

    df = pd.read_parquet(bench_parquet)
    log.info("Loaded %d rows from %s", len(df), bench_parquet.name)

    pocket_vals = []
    total = len(df)
    n_ok = n_skip = n_missing_ref = n_missing_mdl = 0

    for i, (_, row) in enumerate(df.iterrows()):
        if i % 500 == 0:
            log.info("  pocket RMSD %d/%d ...", i, total)

        if row["status"] != "ok":
            pocket_vals.append(None)
            n_skip += 1
            continue

        pdb_id    = row["pdb_id"]
        model_idx = int(row["model_idx"])

        ref_path = _search_local(REFS_DIR, pdb_id)
        if ref_path is None:
            pocket_vals.append(None)
            n_missing_ref += 1
            continue

        mdl_path = DATA_ROOT / PRODUCER / pdb_id / f"model_{model_idx:03d}.cif"
        if not mdl_path.exists():
            pocket_vals.append(None)
            n_missing_mdl += 1
            continue

        pocket_vals.append(
            compute_pocket_rmsd(str(mdl_path), str(ref_path), pdb_id=pdb_id)
        )
        n_ok += 1

    df["pocket_rmsd"] = pocket_vals
    n_valid = sum(v is not None for v in pocket_vals)
    log.info(
        "Pocket RMSD done: %d ok, %d missing ref, %d missing mdl, %d non-ok skipped",
        n_ok, n_missing_ref, n_missing_mdl, n_skip,
    )
    log.info("  valid pocket_rmsd: %d / %d  (%.1f%%)",
             n_valid, total, 100 * n_valid / total if total else 0)
    return df


# ===========================================================================
# Stage 4 — Export CSVs
# Exact format and grouping logic of dynamicbind_asd_pocket_rmsd_and_eval.py
# ===========================================================================

_RENAME = {
    "bisy_rmsd":   "pose rmsd",
    "pocket_rmsd": "pocket rmsd",
    "lddt_pli":    "lddt-pli",
    "qs_global":   "qs score",
}
_METRIC_COLS = ["bisy_rmsd", "pocket_rmsd", "qs_global", "lddt_pli"]
_OUT_COLS    = ["id", "pose rmsd", "pocket rmsd", "qs score", "lddt-pli", "confidence"]


def export_csvs(df: pd.DataFrame) -> None:
    log.info("=" * 60)
    log.info("Stage 4: Export CSVs -> %s", EVAL_DIR)
    EVAL_DIR.mkdir(parents=True, exist_ok=True)

    ok = df[df["status"] == "ok"].copy()
    ok["id"] = ok["pdb_id"].str.lower()
    avail = [c for c in _METRIC_COLS if c in ok.columns]

    def _finalise(frame: pd.DataFrame) -> pd.DataFrame:
        frame = frame.rename(columns=_RENAME)
        cols = [c for c in _OUT_COLS if c in frame.columns]
        return frame[cols].sort_values("id").reset_index(drop=True)

    for producer, grp in ok.groupby("producer"):

        # ── best: min bisy_rmsd per compound; unscored → empty metrics row ──
        scored = grp.dropna(subset=["bisy_rmsd"])
        if len(scored) > 0:
            best_idx = scored.groupby("id")["bisy_rmsd"].idxmin().dropna().astype(int)
            best = grp.loc[best_idx].copy()
        else:
            best = grp.iloc[0:0].copy()
        unscored_ids = set(grp["id"].unique()) - set(best["id"].unique())
        if unscored_ids:
            fallback = (grp[grp["id"].isin(unscored_ids)]
                        .groupby("id", as_index=False).first())
            best = pd.concat([best, fallback], ignore_index=True)
        best = _finalise(best)
        path = EVAL_DIR / f"{producer}alloligand_best.csv"
        best.to_csv(path, index=False)
        log.info("best        %s: %d rows -> %s", producer, len(best), path.name)

        # ── avg: mean over all models per compound ───────────────────────────
        avg = (grp.groupby("id")[avail + ["confidence"]]
                  .mean()
                  .rename(columns=_RENAME)
                  .reset_index())
        avg_cols = [c for c in _OUT_COLS if c in avg.columns]
        avg = avg[avg_cols].sort_values("id").reset_index(drop=True)
        path = EVAL_DIR / f"{producer}alloligand_avg.csv"
        avg.to_csv(path, index=False)
        log.info("avg         %s: %d rows -> %s", producer, len(avg), path.name)

        # ── clean best: min bisy_rmsd among models with pocket_rmsd ≤ thr ───
        # Fallback mirrors `best`: compounds with no bisy_rmsd get their first
        # model in the clean set so that clean_best and clean_avg row counts
        # are always consistent.
        grp_clean = grp[
            grp["pocket_rmsd"].notna() & (grp["pocket_rmsd"] <= POCKET_CLEAN_THR)
        ]
        if len(grp_clean) > 0:
            clean_scored = grp_clean.dropna(subset=["bisy_rmsd"])
            if len(clean_scored) > 0:
                ci = (clean_scored.groupby("id")["bisy_rmsd"]
                                  .idxmin().dropna().astype(int))
                clean_best = grp_clean.loc[ci].copy()
            else:
                clean_best = grp_clean.iloc[0:0].copy()
            # fallback: compounds in clean set but missing bisy_rmsd
            unscored_clean = (set(grp_clean["id"].unique())
                              - set(clean_best["id"].unique()))
            if unscored_clean:
                fallback_clean = (grp_clean[grp_clean["id"].isin(unscored_clean)]
                                  .groupby("id", as_index=False).first())
                clean_best = pd.concat([clean_best, fallback_clean], ignore_index=True)
            clean_best = _finalise(clean_best)
        else:
            clean_best = pd.DataFrame(columns=best.columns)
        path = EVAL_DIR / f"{producer}alloligand_clean_best.csv"
        clean_best.to_csv(path, index=False)
        log.info("clean_best  %s: %d rows -> %s", producer, len(clean_best), path.name)

        # ── clean avg: mean over models with pocket_rmsd ≤ thr ───────────────
        if len(grp_clean) > 0:
            clean_avg = (grp_clean.groupby("id")[avail + ["confidence"]]
                                  .mean()
                                  .rename(columns=_RENAME)
                                  .reset_index())
            ca_cols = [c for c in _OUT_COLS if c in clean_avg.columns]
            clean_avg = clean_avg[ca_cols].sort_values("id").reset_index(drop=True)
        else:
            clean_avg = pd.DataFrame(columns=avg.columns)
        path = EVAL_DIR / f"{producer}alloligand_clean_avg.csv"
        clean_avg.to_csv(path, index=False)
        log.info("clean_avg   %s: %d rows -> %s", producer, len(clean_avg), path.name)

    log.info("All CSVs written to: %s", EVAL_DIR)


# ===========================================================================
# Entry point
# ===========================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="DynamicBind ASD v2 pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Stage numbers:\n"
            "  1  Normalize raw outputs\n"
            "  2  Score with OpenStructure\n"
            "  3  Pocket RMSD\n"
            "  4  Export CSVs\n\n"
            "Example — rescore only (skip normalization):\n"
            "  python dynamicbind_asd_v2_pipeline.py --from-stage 2"
        ),
    )
    parser.add_argument(
        "--from-stage", type=int, default=1, choices=[1, 2, 3, 4],
        metavar="{1,2,3,4}",
        help="Start from this stage (default: 1). "
             "Use 2 to skip normalization when data is already prepared.",
    )
    args = parser.parse_args()
    from_stage = args.from_stage

    # Stage 1 — normalize raw DynamicBind ASD v2 outputs
    if from_stage <= 1:
        entries = run_normalization()
        if not entries:
            log.error("No compounds normalized — check RAW_ROOT: %s", RAW_ROOT)
            sys.exit(1)
        # Stage 1b — pre-warm references so score_tree can find them by compound_key
        prewarm_refs()
    else:
        log.info("Skipping Stage 1 (normalization) — starting from Stage %d", from_stage)

    # Stage 2 — score with OpenStructure (SDF fix + QSEntity fix already in scoring.py)
    if from_stage <= 2:
        bench_parquet = run_scoring()
    else:
        bench_parquet = OUTPUT_DIR / "benchmark.parquet"
        if not bench_parquet.exists():
            log.error("benchmark.parquet not found at %s — run Stage 2 first", bench_parquet)
            sys.exit(1)
        log.info("Skipping Stage 2 — using existing %s", bench_parquet)

    # Stage 3 — pocket RMSD (Bug-A + Bug-B fixed)
    if from_stage <= 3:
        df_enriched = run_pocket_rmsd(bench_parquet)
        enriched_path = OUTPUT_DIR / "benchmark_with_pocket_rmsd.parquet"
        df_enriched.to_parquet(enriched_path, index=False)
        log.info("Saved enriched parquet -> %s", enriched_path)
    else:
        enriched_path = OUTPUT_DIR / "benchmark_with_pocket_rmsd.parquet"
        if not enriched_path.exists():
            log.error("benchmark_with_pocket_rmsd.parquet not found — run Stage 3 first")
            sys.exit(1)
        df_enriched = pd.read_parquet(enriched_path)
        log.info("Skipping Stage 3 — loaded %d rows from %s", len(df_enriched), enriched_path)

    # Stage 4 — export CSVs
    export_csvs(df_enriched)

    # ── Quick summary ──────────────────────────────────────────────────────
    log.info("=" * 60)
    ok_df = df_enriched[df_enriched["status"] == "ok"]
    n_compounds = ok_df["pdb_id"].nunique()
    n_models    = len(ok_df)
    bisy_cov    = ok_df["bisy_rmsd"].notna().mean()
    prmsd_cov   = ok_df["pocket_rmsd"].notna().mean()
    log.info("Summary: %d compounds, %d models", n_compounds, n_models)
    log.info("  bisy_rmsd coverage  : %.1f%%", 100 * bisy_cov)
    log.info("  pocket_rmsd coverage: %.1f%%", 100 * prmsd_cov)
    log.info("Pipeline complete.")
