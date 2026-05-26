"""
dynamicbind_pla_pocket_rmsd_and_eval.py
========================================
DYNAMICBIND PLA — STANDALONE POST-PROCESSING SCRIPT
----------------------------------------------------
This script is EXCLUSIVELY for the DynamicBind PLA producer.
The raw data lives in allosteric_dynamicbind_results/ and was scored
under the producer name "dynamicbind_asd" via the normalize pipeline,
but the compound keys ({uniprot}_{lig}, e.g. b6ywb8_dtp) are confirmed
to be the same PLA benchmark compound keys used by af3_pla, boltz_pla,
chai_pla, and protenix_pla (100% overlap with af3_pla's 168 keys; 320
total including extended compounds).

What it does:
  1. Loads  plb_bench_output/asd/benchmark.parquet  (already has bisy_rmsd,
     lddt_pli, qs_global from the score_normalized_tree step).
  2. Computes pocket_rmsd for every row using Kabsch superposition →
     5 Å pocket Cα RMSD (matching plb_bench_run_pla.ipynb).
  3. Saves the enriched table to
       plb_bench_output/asd/benchmark_with_pocket_rmsd.parquet
  4. Writes per-producer eval CSVs to evalspreadsheets/pla/:
       dynamicbind_plamainligand_best.csv
       dynamicbind_plamainligand_avg.csv
       dynamicbind_plamainligand_clean_best.csv
       dynamicbind_plamainligand_clean_avg.csv
     Column order matches existing PLA eval CSVs:
       id, pose rmsd, pocket rmsd, qs score, lddt-pli, confidence
     Clean versions exclude models with pocket_rmsd > 20.0 Å.

Run from WSL with the plb conda env active:
    python -u scoring/dynamicbind_pla_pocket_rmsd_and_eval.py 2>&1 | tee /mnt/c/Temp/dynamicbind_pla_pocket.log
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
OUTPUT_DIR = BENCH_ROOT / "plb_bench_output" / "asd"   # where the parquet was scored
EVAL_DIR   = BENCH_ROOT / "evalspreadsheets" / "pla"   # matches af3_pla, boltz_pla etc.

POCKET_RADIUS = 5.0   # Å
# Clean versions exclude models with pocket_rmsd above this threshold.
# Confirmed from data: 17.13 Å is included in clean, 31.1 Å is excluded.
POCKET_CLEAN_THRESHOLD = 20.0  # Å

# Output producer name — the parquet uses "dynamicbind_asd" internally but
# the eval CSVs are written as "dynamicbind_pla" to match the PLA pipeline.
OUTPUT_PRODUCER = "dynamicbind_pla"

# ---------------------------------------------------------------------------
# Solvent / buffer exclusion list (verbatim from plb_bench_run_pla.ipynb)
# ---------------------------------------------------------------------------
_SOLVENT_EXT = {
    "HOH", "DOD", "WAT", "EDO", "GOL", "PEG", "PO4", "SO4", "ACT",
    "MES", "TRIS", "BME", "DTT", "MPD", "FMT", "ACE", "NH2",
    "MG", "ZN", "CA", "NA", "CL", "K", "MN", "FE", "CU", "CO",
}

# ---------------------------------------------------------------------------
# OST helpers (matching plb_bench_run_pla.ipynb)
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


def _get_all_lig_instances(ent, lig_name=None):
    """
    Return a list of xyz arrays — one per ligand copy found in the structure.
    Named copies (matching lig_name) are preferred over the fallback (largest
    non-solvent residue).  Returns [] if nothing is found.
    (Bug A fix from plb_bench_run_pla.ipynb cell 83ecd346)
    """
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


def _pick_best_ref_lig(ref_copies, anchor):
    """
    From all reference ligand copies, return the one whose centroid is
    closest to `anchor` (a 3D point in reference coordinate space).
    (Bug A fix from plb_bench_run_pla.ipynb cell 83ecd346)
    """
    if not ref_copies:
        return _np.empty((0, 3))
    if len(ref_copies) == 1:
        return ref_copies[0]
    return min(ref_copies,
               key=lambda pts: _np.linalg.norm(pts.mean(axis=0) - anchor))


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
    Pocket RMSD — verbatim from plb_bench_run_pla.ipynb cell 83ecd346.

    Two bugs fixed vs the earlier ASD notebook version:

    Bug A — wrong ligand copy (multi-subunit assemblies):
      The reference ligand copy used to define the pocket is selected by
      proximity to the reference CA centroid (same coordinate space), not
      the raw model ligand centroid.

    Bug B — index-matched pocket atom selection (Boltz / any model with
      non-PDB residue numbering):
      After Kabsch alignment, transform ALL model CAs into reference space,
      then for each reference pocket CA find the NEAREST transformed model CA.
      Completely immune to residue-number offsets.

    Algorithm:
      1. Collect all reference and model ligand copies.
      2. Collect CA atoms; run chain-aware matching and Kabsch.
      3. Select the reference ligand copy closest to the matched reference CA centroid.
      4. Define pocket: reference CA within `radius` Å of that ligand copy.
      5. Transform ALL model CAs into reference space.
      6. For each reference pocket CA, find the nearest transformed model CA.
      7. Compute RMSD between paired atoms.
    """
    try:
        mdl = _ost_load_cif(mdl_cif_path)
        ref = _ost_load_cif(ref_cif_path)

        parts = str(pdb_id).split("_")
        lig_name = parts[-1].upper() if len(parts) >= 2 else None

        # Step 1: collect all ligand copies
        ref_copies = _get_all_lig_instances(ref, lig_name)
        if not ref_copies:
            return None
        mdl_copies = _get_all_lig_instances(mdl, lig_name)
        if not mdl_copies:
            return None

        # Step 2: CA matching and Kabsch
        ref_by_chain = _collect_ca_by_chain(ref)
        mdl_by_chain = _collect_ca_by_chain(mdl)
        ref_pts, mdl_pts = _match_ca(ref_by_chain, mdl_by_chain)
        if ref_pts is None or len(ref_pts) < 3:
            return None
        R, t = _kabsch(ref_pts, mdl_pts)

        # Step 3: pick ref ligand copy closest to ref CA centroid (Bug A fix)
        ref_ca_centroid = ref_pts.mean(axis=0)
        ref_lig = _pick_best_ref_lig(ref_copies, ref_ca_centroid)
        if ref_lig.shape[0] == 0:
            return None

        # Step 4: define reference pocket
        mask = _np.array([
            _np.linalg.norm(ref_lig - pos, axis=1).min() <= radius
            for pos in ref_pts
        ])
        if mask.sum() < 3:
            return None
        pocket_ref = ref_pts[mask]

        # Step 5: transform ALL model CAs into reference space
        all_mdl_t = (R @ mdl_pts.T).T + t   # shape (N_mdl, 3)

        # Step 6: nearest-neighbour matching (Bug B fix)
        pocket_mdl_t = _np.array([
            all_mdl_t[_np.argmin(_np.linalg.norm(all_mdl_t - ref_pos, axis=1))]
            for ref_pos in pocket_ref
        ])

        # Step 7: pocket RMSD
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

    # 2. Compute pocket_rmsd
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

    # 4. Eval CSVs — matching column order of existing PLA eval CSVs:
    #    id, pose rmsd, pocket rmsd, qs score, lddt-pli, confidence
    _RENAME = {
        "bisy_rmsd":   "pose rmsd",
        "pocket_rmsd": "pocket rmsd",
        "lddt_pli":    "lddt-pli",
        "qs_global":   "qs score",
    }
    # Column order matches af3_plamainligand_best.csv exactly
    _OUT_COLS = ["id", "pose rmsd", "pocket rmsd", "qs score", "lddt-pli", "confidence"]
    _METRIC_COLS = ["bisy_rmsd", "pocket_rmsd", "qs_global", "lddt_pli"]

    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    ok = df[df["status"] == "ok"].copy()
    # id = full compound key (e.g. "b6ywb8_dtp") — do NOT truncate; multiple
    # ligands share the same UniProt prefix and would collapse if truncated.
    ok["id"] = ok["pdb_id"].str.lower()
    avail = [c for c in _METRIC_COLS if c in ok.columns]

    def _make_best(grp):
        """Select min bisy_rmsd model per id; include unscored ids as empty rows."""
        scored = grp.dropna(subset=["bisy_rmsd"])
        if len(scored) > 0:
            best_idx = scored.groupby("id")["bisy_rmsd"].idxmin().dropna().astype(int)
            best_scored = grp.loc[best_idx].copy()
        else:
            best_scored = grp.iloc[0:0].copy()
        unscored_ids = set(grp["id"].unique()) - set(best_scored["id"])
        if unscored_ids:
            unscored = (grp[grp["id"].isin(unscored_ids)]
                        .groupby("id", as_index=False).first())
            return pd.concat([best_scored, unscored], ignore_index=True)
        return best_scored

    def _make_avg(grp):
        """Mean over all models per id."""
        return (grp.groupby("id")[avail + ["confidence"]]
                   .mean().reset_index())

    def _finalise(frame, rename=True):
        if rename:
            frame = frame.rename(columns=_RENAME)
        cols = [c for c in _OUT_COLS if c in frame.columns]
        return frame[cols].sort_values("id").reset_index(drop=True)

    # Only one producer in this parquet (dynamicbind_asd), output as dynamicbind_pla
    grp = ok  # if multiple producers ever appear, loop here

    # ── regular best ────────────────────────────────────────────────────────
    best = _finalise(_make_best(grp))
    best_csv = EVAL_DIR / f"{OUTPUT_PRODUCER}mainligand_best.csv"
    best.to_csv(best_csv, index=False)
    print(f"best        {OUTPUT_PRODUCER}: {len(best)} rows -> {best_csv.name}")

    # ── regular avg ─────────────────────────────────────────────────────────
    avg = _finalise(_make_avg(grp))
    avg_csv = EVAL_DIR / f"{OUTPUT_PRODUCER}mainligand_avg.csv"
    avg.to_csv(avg_csv, index=False)
    print(f"avg         {OUTPUT_PRODUCER}: {len(avg)} rows -> {avg_csv.name}")

    # ── clean best (pocket_rmsd ≤ threshold only) ────────────────────────────
    grp_clean = grp[
        grp["pocket_rmsd"].notna() & (grp["pocket_rmsd"] <= POCKET_CLEAN_THRESHOLD)
    ]
    if len(grp_clean) > 0:
        clean_best = _finalise(_make_best(grp_clean))
    else:
        clean_best = pd.DataFrame(columns=_OUT_COLS)
    clean_best_csv = EVAL_DIR / f"{OUTPUT_PRODUCER}mainligand_clean_best.csv"
    clean_best.to_csv(clean_best_csv, index=False)
    print(f"clean_best  {OUTPUT_PRODUCER}: {len(clean_best)} rows -> {clean_best_csv.name}")

    # ── clean avg (pocket_rmsd ≤ threshold only) ─────────────────────────────
    if len(grp_clean) > 0:
        clean_avg = _finalise(_make_avg(grp_clean))
    else:
        clean_avg = pd.DataFrame(columns=_OUT_COLS)
    clean_avg_csv = EVAL_DIR / f"{OUTPUT_PRODUCER}mainligand_clean_avg.csv"
    clean_avg.to_csv(clean_avg_csv, index=False)
    print(f"clean_avg   {OUTPUT_PRODUCER}: {len(clean_avg)} rows -> {clean_avg_csv.name}")

    print(f"\nAll CSVs written to: {EVAL_DIR}")
