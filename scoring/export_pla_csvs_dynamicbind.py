"""Regenerate the 4 DynamicBind PLA eval CSVs from the corrected parquet.

Reads plb_bench_output/pla/benchmark_with_pocket_rmsd.parquet (dynamicbind_pla,
already has bisy_rmsd / lddt_pli / qs_global / pocket_rmsd against the CORRECTED
references) and writes, keyed by FULL compound id ("uniprot_ligand", matching
af3/chai/boltz):

    dynamicbind_plamainligand_best.csv
    dynamicbind_plamainligand_avg.csv
    dynamicbind_plamainligand_clean_best.csv   (pocket_rmsd <= 20 A)
    dynamicbind_plamainligand_clean_avg.csv    (pocket_rmsd <= 20 A)

Logic mirrors scoring/dynamicbind_pla_pocket_rmsd_and_eval.py exactly. Does NOT
recompute any metric and does NOT touch other producers' CSVs.
"""
from pathlib import Path
import pandas as pd

BENCH = Path("/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark")
PARQUET = BENCH / "plb_bench_output" / "pla" / "benchmark_with_pocket_rmsd.parquet"
EVAL_DIR = BENCH / "evalspreadsheets" / "pla"
OUTPUT_PRODUCER = "dynamicbind_pla"
POCKET_CLEAN_THRESHOLD = 20.0

_RENAME = {
    "bisy_rmsd":   "pose rmsd",
    "pocket_rmsd": "pocket rmsd",
    "lddt_pli":    "lddt-pli",
    "qs_global":   "qs score",
}
_OUT_COLS = ["id", "pose rmsd", "pocket rmsd", "qs score", "lddt-pli", "confidence"]
_METRIC_COLS = ["bisy_rmsd", "pocket_rmsd", "qs_global", "lddt_pli"]

df = pd.read_parquet(PARQUET)
df = df[df["producer"] == OUTPUT_PRODUCER].copy()
ok = df[df["status"] == "ok"].copy()
# FULL compound key (e.g. "b6ywb8_dtp"), matching af3 — do NOT truncate to UniProt
ok["id"] = ok["pdb_id"].str.lower()
avail = [c for c in _METRIC_COLS if c in ok.columns]
print(f"rows: {len(df)}  ok: {len(ok)}  unique compounds: {ok['id'].nunique()}")


def _make_best(grp):
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
    return grp.groupby("id")[avail + ["confidence"]].mean().reset_index()


def _finalise(frame):
    frame = frame.rename(columns=_RENAME)
    cols = [c for c in _OUT_COLS if c in frame.columns]
    return frame[cols].sort_values("id").reset_index(drop=True)


grp = ok

best = _finalise(_make_best(grp))
best.to_csv(EVAL_DIR / f"{OUTPUT_PRODUCER}mainligand_best.csv", index=False)
print(f"best        : {len(best)} rows")

avg = _finalise(_make_avg(grp))
avg.to_csv(EVAL_DIR / f"{OUTPUT_PRODUCER}mainligand_avg.csv", index=False)
print(f"avg         : {len(avg)} rows")

grp_clean = grp[grp["pocket_rmsd"].notna() & (grp["pocket_rmsd"] <= POCKET_CLEAN_THRESHOLD)]
clean_best = _finalise(_make_best(grp_clean)) if len(grp_clean) else pd.DataFrame(columns=_OUT_COLS)
clean_best.to_csv(EVAL_DIR / f"{OUTPUT_PRODUCER}mainligand_clean_best.csv", index=False)
print(f"clean_best  : {len(clean_best)} rows")

clean_avg = _finalise(_make_avg(grp_clean)) if len(grp_clean) else pd.DataFrame(columns=_OUT_COLS)
clean_avg.to_csv(EVAL_DIR / f"{OUTPUT_PRODUCER}mainligand_clean_avg.csv", index=False)
print(f"clean_avg   : {len(clean_avg)} rows")

print("done ->", EVAL_DIR)
