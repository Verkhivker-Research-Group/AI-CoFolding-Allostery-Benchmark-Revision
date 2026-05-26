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
# Stage 2 — Score with OpenStructure
# ============================================================
log.info("=" * 60)
log.info("Stage 2: Score dynamicbind_pla with OST")

from plb_bench.score_tree import ScoreConfig, run as score_run  # noqa: E402

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
# score_run writes OUTPUT_DIR/benchmark.parquet (dynamicbind_pla rows only)
new_parquet = score_run(cfg)
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


def _uniprot_id(compound_key: str) -> str:
    return compound_key.split("_")[0].lower()


# Reload the full merged parquet so all producers get fresh CSVs
full = pd.read_parquet(prmsd_path)
pla_ok = full[full["status"] == "ok"].copy()
pla_ok["id"] = pla_ok["pdb_id"].apply(_uniprot_id)
avail = [c for c in _PLA_METRIC_COLS if c in pla_ok.columns]

for producer, grp in pla_ok.groupby("producer"):
    # best
    scored = grp.dropna(subset=["bisy_rmsd"])
    best_idx = scored.groupby("id")["bisy_rmsd"].idxmin().dropna().astype(int)
    best = grp.loc[best_idx].copy()
    unscored_ids = set(grp["id"].unique()) - set(best["id"].unique())
    if unscored_ids:
        fallback = grp[grp["id"].isin(unscored_ids)].groupby("id", as_index=False).first()
        best = pd.concat([best, fallback], ignore_index=True)
    best = best.rename(columns=_PLA_RENAME)
    best_cols = ["id"] + [_PLA_RENAME[c] for c in avail if _PLA_RENAME[c] in best.columns] + ["confidence"]
    best = best[[c for c in best_cols if c in best.columns]].sort_values("id").reset_index(drop=True)
    best.to_csv(EVAL_DIR / f"{producer}mainligand_best.csv", index=False)
    log.info("best  %s: %d rows", producer, len(best))

    # avg
    avg = (grp.groupby("id")[avail + ["confidence"]].mean()
              .rename(columns=_PLA_RENAME).reset_index())
    avg_cols = ["id"] + [_PLA_RENAME[c] for c in avail if _PLA_RENAME[c] in avg.columns] + ["confidence"]
    avg = avg[[c for c in avg_cols if c in avg.columns]].sort_values("id").reset_index(drop=True)
    avg.to_csv(EVAL_DIR / f"{producer}mainligand_avg.csv", index=False)
    log.info("avg   %s: %d rows", producer, len(avg))

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
