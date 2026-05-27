"""
fix_qs_dynamicbind.py
=====================
DynamicBind-ONLY QS-score resolution.

OST 2.11's QSScorer/QSEntity is polymer-only: it strips the ligand, so a
single-chain DynamicBind model gets a degenerate binary QS (1.0 if the one
protein chain maps, else 0.0). This script replaces the DynamicBind qs_global
with a continuous *ligand-inclusive* interface QS that reproduces the intent of
the main PLB-bench pipeline's QS (protein-ligand interface contact conservation).

Method (per model vs the CORRECTED reference, same protein):
  1. Superpose model protein onto reference via matched CA (Kabsch).
  2. Map each model residue -> nearest reference residue (post-superposition).
  3. ref contact vector  : ref residue in contact with REF ligand (heavy-atom
                           min-dist <= contact_d).
  4. model contact vector: mapped model residue in contact with MODEL ligand.
  5. QS = sum_shared w / (sum_shared w + n_nonshared),
         w = max(0, 1 - |d_ref - d_mdl| / contact_d).

Only touches producer == 'dynamicbind_pla'. Writes qs_global back into
benchmark.parquet and benchmark_with_pocket_rmsd.parquet, then regenerates the 4
DynamicBind eval CSVs (full compound key, matching af3).

Run:
    conda activate plb
    python -u scoring/fix_qs_dynamicbind.py
"""
from __future__ import annotations
import sys, os, tempfile, logging, multiprocessing as _mp
import concurrent.futures as _cf
from pathlib import Path
import numpy as np

_PKG = Path("/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark/scoring/plb_bench")
if str(_PKG) not in sys.path:
    sys.path.insert(0, str(_PKG))

import pandas as pd
from plb_bench.references import get_reference, read_reference_text
from plb_bench.scoring import _get_ref_lig_name

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    stream=sys.stdout, force=True)
log = logging.getLogger(__name__)

BENCH      = Path("/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark")
DATA_ROOT  = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive/plb_bench_data")
REFS_DIR   = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive/references_ref_cifs")
OUT_DIR    = BENCH / "plb_bench_output" / "pla"
EVAL_DIR   = BENCH / "evalspreadsheets" / "pla"
PRODUCER   = "dynamicbind_pla"
CONTACT_D  = 12.0
POCKET_CLEAN_THR = 20.0

# ── OST helpers ──────────────────────────────────────────────────────────────
def _load(text):
    from ost import io as ost_io
    with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False) as f:
        f.write(text); tmp = f.name
    try:
        r = ost_io.LoadMMCIF(tmp, seqres=True, info=True, fault_tolerant=True)
        return r[0] if isinstance(r, tuple) else r
    finally:
        try: os.unlink(tmp)
        except OSError: pass


def _poly_cas(ent):
    out = []
    for ch in ent.chains:
        for r in ch.residues:
            ca = r.FindAtom("CA")
            if not ca.IsValid():
                continue
            heavy = np.array([[a.pos.x, a.pos.y, a.pos.z]
                              for a in r.atoms if a.element != "H"])
            out.append(((ch.name, r.number.num),
                        np.array([ca.pos.x, ca.pos.y, ca.pos.z]), heavy))
    return out


def _lig_heavy(ent, names):
    pts = []
    for ch in ent.chains:
        for r in ch.residues:
            if r.name in names:
                for a in r.atoms:
                    if a.element != "H":
                        pts.append([a.pos.x, a.pos.y, a.pos.z])
    return np.array(pts) if pts else None


def _kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Qc - R @ Pc


def _min_dist(res_heavy, lig):
    if res_heavy is None or res_heavy.size == 0 or lig is None:
        return np.inf
    from scipy.spatial import distance
    return distance.cdist(res_heavy, lig).min()


def qs_ligand_inclusive(mdl, ref, ref_lig, contact_d=CONTACT_D):
    rcas = _poly_cas(ref); mcas = _poly_cas(mdl)
    if len(rcas) < 3 or len(mcas) < 3:
        return None
    rl = _lig_heavy(ref, {ref_lig}); ml = _lig_heavy(mdl, {"LIG"})
    if rl is None or ml is None:
        return None
    rmap = {k: (ca, hv) for k, ca, hv in rcas}
    mmap = {k: (ca, hv) for k, ca, hv in mcas}
    common = [k for k in rmap if k in mmap]
    if len(common) >= 3:
        Rca = np.array([rmap[k][0] for k in common])
        Mca = np.array([mmap[k][0] for k in common])
    else:
        n = min(len(rcas), len(mcas))
        Rca = np.array([rcas[i][1] for i in range(n)])
        Mca = np.array([mcas[i][1] for i in range(n)])
    if len(Rca) < 3:
        return None
    R, t = _kabsch(Mca, Rca)                       # model -> ref frame
    ml_t = (R @ ml.T).T + t
    ref_d = {k: _min_dist(hv, rl) for k, _, hv in rcas}
    refkeys = [k for k, _, _ in rcas]
    refca = np.array([ca for _, ca, _ in rcas])
    mdl_d_byref = {}
    for k, ca, hv in mcas:
        ca_t = (R @ ca) + t
        j = int(np.argmin(np.linalg.norm(refca - ca_t, axis=1)))
        d = _min_dist((R @ hv.T).T + t, ml_t)
        rk = refkeys[j]
        if rk not in mdl_d_byref or d < mdl_d_byref[rk]:
            mdl_d_byref[rk] = d
    keys = (set(k for k, d in ref_d.items() if d <= contact_d) |
            set(k for k, d in mdl_d_byref.items() if d <= contact_d))
    if not keys:
        return None
    shared_w = 0.0; nonshared = 0
    for k in keys:
        dr = ref_d.get(k, np.inf); dm = mdl_d_byref.get(k, np.inf)
        if dr <= contact_d and dm <= contact_d:
            shared_w += max(0.0, 1.0 - abs(dr - dm) / contact_d)
        else:
            nonshared += 1
    denom = shared_w + nonshared
    return float(shared_w / denom) if denom > 0 else None


# ── per-target worker (loads reference once, scores all its models) ──────────
def _score_target(args):
    pdb_id, idxs = args
    try:
        ref_path, _ = get_reference(pdb_id, REFS_DIR, allow_download=False)
        ref = _load(read_reference_text(ref_path))
        ref_lig = _get_ref_lig_name(ref)
    except Exception as e:
        log.warning("[%s] reference load failed: %s", pdb_id, e)
        return {}
    out = {}
    for idx in idxs:
        mp = DATA_ROOT / PRODUCER / pdb_id / f"model_{idx:03d}.cif"
        if not mp.exists():
            continue
        try:
            mdl = _load(mp.read_text())
            out[idx] = qs_ligand_inclusive(mdl, ref, ref_lig)
        except Exception as e:
            log.debug("[%s/%d] qs failed: %s", pdb_id, idx, e)
            out[idx] = None
    return {(pdb_id, k): v for k, v in out.items()}


# ── CSV export (full compound key, matching af3) ─────────────────────────────
_RENAME = {"bisy_rmsd": "pose rmsd", "pocket_rmsd": "pocket rmsd",
           "lddt_pli": "lddt-pli", "qs_global": "qs score"}
_OUT_COLS = ["id", "pose rmsd", "pocket rmsd", "qs score", "lddt-pli", "confidence"]
_METRIC_COLS = ["bisy_rmsd", "pocket_rmsd", "qs_global", "lddt_pli"]


def _export_csvs(df):
    ok = df[(df["status"] == "ok") & (df["producer"] == PRODUCER)].copy()
    ok["id"] = ok["pdb_id"].str.lower()
    avail = [c for c in _METRIC_COLS if c in ok.columns]

    def make_best(grp):
        scored = grp.dropna(subset=["bisy_rmsd"])
        if len(scored):
            bi = scored.groupby("id")["bisy_rmsd"].idxmin().dropna().astype(int)
            best = grp.loc[bi].copy()
        else:
            best = grp.iloc[0:0].copy()
        miss = set(grp["id"].unique()) - set(best["id"])
        if miss:
            fb = grp[grp["id"].isin(miss)].groupby("id", as_index=False).first()
            best = pd.concat([best, fb], ignore_index=True)
        return best

    def make_avg(grp):
        return grp.groupby("id")[avail + ["confidence"]].mean().reset_index()

    def fin(frame):
        frame = frame.rename(columns=_RENAME)
        return frame[[c for c in _OUT_COLS if c in frame.columns]].sort_values("id").reset_index(drop=True)

    fin(make_best(ok)).to_csv(EVAL_DIR / f"{PRODUCER}mainligand_best.csv", index=False)
    fin(make_avg(ok)).to_csv(EVAL_DIR / f"{PRODUCER}mainligand_avg.csv", index=False)
    clean = ok[ok["pocket_rmsd"].notna() & (ok["pocket_rmsd"] <= POCKET_CLEAN_THR)]
    cb = fin(make_best(clean)) if len(clean) else pd.DataFrame(columns=_OUT_COLS)
    ca = fin(make_avg(clean)) if len(clean) else pd.DataFrame(columns=_OUT_COLS)
    cb.to_csv(EVAL_DIR / f"{PRODUCER}mainligand_clean_best.csv", index=False)
    ca.to_csv(EVAL_DIR / f"{PRODUCER}mainligand_clean_avg.csv", index=False)
    log.info("CSVs regenerated: best/avg + clean_best/clean_avg")


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    prmsd = OUT_DIR / "benchmark_with_pocket_rmsd.parquet"
    bench = OUT_DIR / "benchmark.parquet"
    df = pd.read_parquet(prmsd)
    db = df[df["producer"] == PRODUCER]
    log.info("DynamicBind rows: %d", len(db))

    # group model indices per target
    targets = {}
    for _, r in db.iterrows():
        if r["status"] != "ok":
            continue
        targets.setdefault(r["pdb_id"], []).append(int(r["model_idx"]))
    tasks = [(pid, idxs) for pid, idxs in targets.items()]
    log.info("targets: %d", len(tasks))

    qs_map = {}
    nworkers = max(1, (os.cpu_count() or 4) - 1)
    log.info("computing ligand-inclusive QS with %d workers ...", nworkers)
    ctx = _mp.get_context("fork")
    with _cf.ProcessPoolExecutor(max_workers=nworkers, mp_context=ctx) as ex:
        futs = [ex.submit(_score_target, t) for t in tasks]
        done = 0
        for fut in _cf.as_completed(futs):
            try:
                qs_map.update(fut.result())
            except Exception as e:
                log.warning("task failed: %s", e)
            done += 1
            if done % 25 == 0:
                log.info("  %d / %d targets", done, len(tasks))

    # write qs_global back for dynamicbind rows
    def newqs(row):
        if row["producer"] != PRODUCER:
            return row["qs_global"]
        return qs_map.get((row["pdb_id"], int(row["model_idx"])), row["qs_global"])

    for path in (prmsd, bench):
        if not path.exists():
            continue
        d = pd.read_parquet(path)
        d["qs_global"] = d.apply(newqs, axis=1)
        d.to_parquet(path, index=False)
        log.info("updated qs_global -> %s", path.name)

    # CSV (from the pocket parquet which has pocket_rmsd)
    _export_csvs(pd.read_parquet(prmsd))

    upd = pd.read_parquet(prmsd)
    q = upd[upd["producer"] == PRODUCER]["qs_global"].dropna()
    log.info("NEW QS: n=%d mean=%.3f frac(0or1)=%.1f%% coverage=%.1f%%",
             len(q), q.mean(), 100 * ((q == 1) | (q == 0)).mean(),
             100 * upd[upd["producer"] == PRODUCER]["qs_global"].notna().mean())
    log.info("Done.")


if __name__ == "__main__":
    main()
