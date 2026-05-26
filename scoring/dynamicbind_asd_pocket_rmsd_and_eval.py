"""
dynamicbind_asd_pocket_rmsd_and_eval.py
========================================
DYNAMICBIND ASD — STANDALONE POST-PROCESSING SCRIPT
----------------------------------------------------
This script is EXCLUSIVELY for the DynamicBind ASD producer (compound keys of
the form {uniprot}_{lig}, e.g. b6ywb8_dtp).  It does NOT touch any other
producer's data.

What it does (matching plb_bench_run.ipynb cells aedfa6bd + d86bf3a6 exactly):
  1. Loads  plb_bench_output/asd/benchmark.parquet  (already has bisy_rmsd,
     lddt_pli, qs_global from the score_normalized_tree step).
  2. Computes pocket_rmsd for every dynamicbind_asd row using the same
     compute_pocket_rmsd() function defined in plb_bench_run.ipynb cell
     aedfa6bd (Kabsch superposition → 5 Å pocket Cα RMSD).
  3. Saves the enriched table to
       plb_bench_output/asd/benchmark_with_pocket_rmsd.parquet
  4. Writes per-producer best/avg eval CSVs to
       evalspreadsheets/allo/dynamicbind_asdalloligand_best.csv
       evalspreadsheets/allo/dynamicbind_asdalloligand_avg.csv
     (same rename dict and grouping logic as plb_bench_run.ipynb cell d86bf3a6)

Run from WSL with the plb conda env active:
    python -u scoring/dynamicbind_asd_pocket_rmsd_and_eval.py 2>&1 | tee /mnt/c/Temp/dynamicbind_asd_pocket.log
"""

import gzip as _gzip
import os as _os
import tempfile as _tmp

import numpy as _np
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BENCH_ROOT = Path("/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark")
DATA_ROOT  = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive/plb_bench_data")
REFS_DIR   = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive/references_ref_cifs")
OUTPUT_DIR = BENCH_ROOT / "plb_bench_output" / "asd"
EVAL_DIR   = BENCH_ROOT / "evalspreadsheets" / "main"

POCKET_RADIUS = 5.0  # Å — matches plb_bench_run.ipynb

# ---------------------------------------------------------------------------
# Solvent / buffer exclusion list (verbatim from plb_bench_run.ipynb)
# ---------------------------------------------------------------------------
_SOLVENT_EXT = {
    "HOH", "DOD", "WAT", "EDO", "GOL", "PEG", "PO4", "SO4", "ACT",
    "MES", "TRIS", "BME", "DTT", "MPD", "FMT", "ACE", "NH2",
    "MG", "ZN", "CA", "NA", "CL", "K", "MN", "FE", "CU", "CO",
}

# ---------------------------------------------------------------------------
# OST helpers (verbatim from plb_bench_run.ipynb cell aedfa6bd)
# ---------------------------------------------------------------------------

def _ost_load_cif(path):
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


def _poly_types():
    from ost import mol
    types = set()
    for name in (
        "CHAINTYPE_POLY_PEPTIDE_L", "CHAINTYPE_POLY_PEPTIDE_D",
        "CHAINTYPE_POLY_DN", "CHAINTYPE_POLY_RN", "CHAINTYPE_POLY",
        "CHAINTYPE_POLY_SAC_D", "CHAINTYPE_POLY_SAC_L", "CHAINTYPE_WATER",
    ):
        ct = getattr(mol, name, None)
        if ct is not None:
            types.add(ct)
    return types


def _collect_ca_by_chain(ent):
    """Returns {chain_name: [(resnum, xyz), ...]}"""
    poly = _poly_types()
    result = {}
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


def _get_ligand_atoms(ent, lig_name=None):
    poly = _poly_types()
    candidates = []
    for chain in ent.chains:
        if chain.chain_type in poly:
            continue
        for res in chain.residues:
            if res.name in _SOLVENT_EXT:
                continue
            pts = _np.array([[at.pos.x, at.pos.y, at.pos.z] for at in res.atoms])
            if pts.shape[0] == 0:
                continue
            candidates.append((res.name.upper(), pts))
    if not candidates:
        return _np.empty((0, 3))
    if lig_name:
        for name, pts in candidates:
            if name == lig_name.upper():
                return pts  # first matching copy only
    return max(candidates, key=lambda x: x[1].shape[0])[1]


def _match_ca(ref_by_chain, mdl_by_chain):
    ref_flat = {(cid, rnum): pos
                for cid, residues in ref_by_chain.items()
                for rnum, pos in residues}
    mdl_flat = {(cid, rnum): pos
                for cid, residues in mdl_by_chain.items()
                for rnum, pos in residues}
    # Strategy A: exact (chain, resnum) match
    common = sorted(set(ref_flat) & set(mdl_flat))
    if len(common) >= 10:
        return (_np.array([ref_flat[k] for k in common]),
                _np.array([mdl_flat[k] for k in common]))
    # Strategy B: pair chains by descending size
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


def _kabsch(ref_pts, mdl_pts):
    rc, mc = ref_pts.mean(0), mdl_pts.mean(0)
    H = (mdl_pts - mc).T @ (ref_pts - rc)
    U, _, Vt = _np.linalg.svd(H)
    d = _np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ _np.diag([1.0, 1.0, d]) @ U.T
    return R, rc - R @ mc


def compute_pocket_rmsd(mdl_cif_path, ref_cif_path, pdb_id="", radius=POCKET_RADIUS):
    """
    Pocket RMSD — verbatim from plb_bench_run.ipynb cell aedfa6bd.
    Kabsch-superpose all Cα, then RMSD of pocket Cα (within `radius` Å of ref ligand).
    """
    try:
        mdl = _ost_load_cif(mdl_cif_path)
        ref = _ost_load_cif(ref_cif_path)
        parts = str(pdb_id).split("_")
        lig_name = parts[-1].upper() if len(parts) >= 2 else None
        ref_lig = _get_ligand_atoms(ref, lig_name)
        if ref_lig.shape[0] == 0:
            return None
        ref_by_chain = _collect_ca_by_chain(ref)
        mdl_by_chain = _collect_ca_by_chain(mdl)
        ref_pts, mdl_pts = _match_ca(ref_by_chain, mdl_by_chain)
        if ref_pts is None or len(ref_pts) < 3:
            return None
        R, t = _kabsch(ref_pts, mdl_pts)
        mask = _np.array([
            _np.linalg.norm(ref_lig - pos, axis=1).min() <= radius
            for pos in ref_pts
        ])
        if mask.sum() < 3:
            return None
        pocket_ref   = ref_pts[mask]
        pocket_mdl_t = (R @ mdl_pts[mask].T).T + t
        return float(_np.sqrt(_np.mean(_np.sum((pocket_ref - pocket_mdl_t) ** 2, axis=1))))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    from plb_bench.references import _search_local

    # 1. Load parquet produced by score_normalized_tree for dynamicbind_asd
    parquet_in = OUTPUT_DIR / "benchmark.parquet"
    print(f"Loading {parquet_in} ...")
    df = pd.read_parquet(parquet_in)
    print(f"  {len(df)} rows, columns: {list(df.columns)}")

    # 2. Compute pocket_rmsd (matching plb_bench_run.ipynb cell 62059fdc)
    pocket_vals = []
    total = len(df)
    for i, (_, row) in enumerate(df.iterrows()):
        if i % 100 == 0:
            print(f"  pocket RMSD {i}/{total} ...")
        if row["status"] != "ok":
            pocket_vals.append(None)
            continue
        mdl_path = (DATA_ROOT / row["producer"] / row["pdb_id"]
                    / f"model_{int(row['model_idx']):03d}.cif")
        ref_path = _search_local(REFS_DIR, row["pdb_id"])
        if not mdl_path.exists() or ref_path is None:
            pocket_vals.append(None)
            continue
        pocket_vals.append(
            compute_pocket_rmsd(mdl_path, ref_path, pdb_id=row["pdb_id"])
        )

    df["pocket_rmsd"] = pocket_vals
    n_ok = sum(v is not None for v in pocket_vals)
    print(f"\npocket_rmsd: {n_ok} ok, {total - n_ok} missing/failed")

    # 3. Save enriched parquet
    enriched = OUTPUT_DIR / "benchmark_with_pocket_rmsd.parquet"
    df.to_parquet(enriched, index=False)
    print(f"Saved -> {enriched}")

    # 4. Eval CSVs (verbatim rename dict + grouping from plb_bench_run.ipynb cell d86bf3a6)
    _RENAME = {
        "bisy_rmsd":   "pose rmsd",
        "pocket_rmsd": "pocket rmsd",
        "lddt_pli":    "lddt-pli",
        "qs_global":   "qs score",
    }
    _METRIC_COLS = ["bisy_rmsd", "pocket_rmsd", "qs_global", "lddt_pli"]

    # "Clean" versions exclude models with pocket_rmsd above this threshold
    # (threshold inferred from pipeline data: 17.12 Å included, 31.1 Å excluded)
    POCKET_CLEAN_THRESHOLD = 20.0  # Å

    # For ASD pipeline the compound key IS the unique identifier (e.g. "b6ywb8_dtp").
    # Do NOT truncate — multiple ligands share the same UniProt prefix and would collapse.
    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    ok = df[df["status"] == "ok"].copy()
    ok["id"] = ok["pdb_id"].str.lower()
    avail = [c for c in _METRIC_COLS if c in ok.columns]

    for producer, grp in ok.groupby("producer"):
        # ── regular best (min bisy_rmsd across all models) ──────────────────
        # Compounds where all models have NaN bisy_rmsd are included as empty rows.
        scored = grp.dropna(subset=["bisy_rmsd"])
        if len(scored) > 0:
            best_idx = scored.groupby("id")["bisy_rmsd"].idxmin().dropna().astype(int)
            best_scored = grp.loc[best_idx].copy()
        else:
            best_scored = grp.iloc[0:0].copy()
        # Compounds with no scored models: pick first model row, metrics will be NaN
        unscored_ids = set(grp["id"].unique()) - set(best_scored["id"])
        if unscored_ids:
            unscored = (grp[grp["id"].isin(unscored_ids)]
                        .groupby("id", as_index=False).first())
            best = pd.concat([best_scored, unscored], ignore_index=True)
        else:
            best = best_scored
        best = best.rename(columns=_RENAME)
        best_cols = (["id"]
                     + [_RENAME[c] for c in avail if _RENAME[c] in best.columns]
                     + ["confidence"])
        best = (best[[c for c in best_cols if c in best.columns]]
                .sort_values("id").reset_index(drop=True))
        best_csv = EVAL_DIR / f"{producer}alloligand_best.csv"
        best.to_csv(best_csv, index=False)
        print(f"best        {producer}: {len(best)} rows -> {best_csv.name}")

        # ── regular avg (mean across all models) ────────────────────────────
        avg = (grp.groupby("id")[avail + ["confidence"]]
                  .mean().rename(columns=_RENAME).reset_index())
        avg_cols = (["id"]
                    + [_RENAME[c] for c in avail if _RENAME[c] in avg.columns]
                    + ["confidence"])
        avg = (avg[[c for c in avg_cols if c in avg.columns]]
               .sort_values("id").reset_index(drop=True))
        avg_csv = EVAL_DIR / f"{producer}alloligand_avg.csv"
        avg.to_csv(avg_csv, index=False)
        print(f"avg         {producer}: {len(avg)} rows -> {avg_csv.name}")

        # ── clean best (min bisy_rmsd among models with pocket_rmsd ≤ threshold) ──
        # Models with NaN pocket_rmsd are excluded from the clean pool.
        # If a compound has NO model surviving the threshold, it is omitted entirely.
        grp_clean = grp[
            grp["pocket_rmsd"].notna() & (grp["pocket_rmsd"] <= POCKET_CLEAN_THRESHOLD)
        ]
        if len(grp_clean) > 0 and "bisy_rmsd" in grp_clean.columns:
            clean_scored = grp_clean.dropna(subset=["bisy_rmsd"])
            if len(clean_scored) > 0:
                clean_best_idx = (clean_scored.groupby("id")["bisy_rmsd"]
                                              .idxmin().dropna().astype(int))
                clean_best = grp_clean.loc[clean_best_idx].copy()
            else:
                clean_best = grp_clean.iloc[0:0].copy()
            clean_best = clean_best.rename(columns=_RENAME)
            clean_best_cols = (["id"]
                               + [_RENAME[c] for c in avail if _RENAME[c] in clean_best.columns]
                               + ["confidence"])
            clean_best = (clean_best[[c for c in clean_best_cols if c in clean_best.columns]]
                          .sort_values("id").reset_index(drop=True))
        else:
            clean_best = pd.DataFrame(columns=best.columns)
        clean_best_csv = EVAL_DIR / f"{producer}alloligand_clean_best.csv"
        clean_best.to_csv(clean_best_csv, index=False)
        print(f"clean_best  {producer}: {len(clean_best)} rows -> {clean_best_csv.name}")

        # ── clean avg (mean over models with pocket_rmsd ≤ threshold) ───────
        if len(grp_clean) > 0:
            clean_avg = (grp_clean.groupby("id")[avail + ["confidence"]]
                                  .mean().rename(columns=_RENAME).reset_index())
            clean_avg_cols = (["id"]
                              + [_RENAME[c] for c in avail if _RENAME[c] in clean_avg.columns]
                              + ["confidence"])
            clean_avg = (clean_avg[[c for c in clean_avg_cols if c in clean_avg.columns]]
                         .sort_values("id").reset_index(drop=True))
        else:
            clean_avg = pd.DataFrame(columns=avg.columns)
        clean_avg_csv = EVAL_DIR / f"{producer}alloligand_clean_avg.csv"
        clean_avg.to_csv(clean_avg_csv, index=False)
        print(f"clean_avg   {producer}: {len(clean_avg)} rows -> {clean_avg_csv.name}")

    print(f"\nAll CSVs written to: {EVAL_DIR}")
