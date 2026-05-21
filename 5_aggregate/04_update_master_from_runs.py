import sys

# Allow running as a standalone script without installing as a package.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

"""Update a master spreadsheet with model-level metrics and/or ligand RMSD measurements.

This consolidates many notebook cells (8-28, 31-37) into a single configurable tool.

Supported update modes:
1) upsert_model_metrics: updates/creates rows for a (source, dataset) from a model-run CSV,
   keeping best (lowest) rmsd per group.
2) update_ligand_rmsd: updates existing rows' RMSD/ligand_rmsd from a ligand RMSD CSV,
   keeping best (lowest) per group.

You control how IDs are extracted via --id_mode.
"""

from __future__ import annotations

import argparse
from typing import Literal, Optional

import pandas as pd

from scripts.utils.csv_ops import best_row_per_group, write_csv


IdMode = Literal[
    "pdb_id_first4",
    "pdb_id_first5",
    "cif_file_first5",
    "complex_seed_first5",
    "pdb_plus_chain",
    "pdb_no_underscore_upper",
]


def _make_id(series: pd.Series, mode: IdMode, chain: Optional[pd.Series] = None) -> pd.Series:
    s = series.astype(str)
    if mode == "pdb_id_first4":
        return s.str[:4].str.upper()
    if mode == "pdb_id_first5":
        return s.str[:5].str.upper()
    if mode == "cif_file_first5":
        return s.str[:5].str.upper()
    if mode == "complex_seed_first5":
        return s.str[:5].str.upper()
    if mode == "pdb_plus_chain":
        if chain is None:
            raise ValueError("pdb_plus_chain requires --chain_col")
        return (series.astype(str).str.strip() + chain.astype(str).str.strip()).str.upper()
    if mode == "pdb_no_underscore_upper":
        return series.astype(str).str.replace("_", "", regex=False).str.upper()
    raise ValueError(f"Unknown id_mode: {mode}")


def upsert_model_metrics(
    master: pd.DataFrame,
    run_df: pd.DataFrame,
    *,
    source: str,
    dataset: str,
    id_col: str,
    id_mode: IdMode,
    rmsd_col: str,
    tm_col: Optional[str],
    plddt_col: Optional[str],
    group_keep: Literal["min", "max"] = "min",
) -> pd.DataFrame:
    # Build PDB_ID
    run_df = run_df.copy()
    run_df["PDB_ID"] = _make_id(run_df[id_col], id_mode)

    # Choose best per PDB_ID
    best = best_row_per_group(run_df, "PDB_ID", rmsd_col, keep=group_keep)

    # Standardize output cols
    best_out = pd.DataFrame({
        "Source": source,
        "dataset": dataset,
        "PDB_ID": best["PDB_ID"],
        "model_rmsd": best[rmsd_col],
    })
    if tm_col and tm_col in best.columns:
        best_out["model_tm_score"] = best[tm_col]
    if plddt_col and plddt_col in best.columns:
        best_out["model_plddt"] = best[plddt_col]

    # Ensure required columns exist in master
    for col in ["model_rmsd", "model_tm_score", "model_plddt", "RMSD", "ligand_rmsd"]:
        if col not in master.columns:
            master[col] = pd.NA

    # Drop existing rows for these PDB_IDs to avoid duplicates, then append
    mask = (
        (master["Source"].astype(str).str.lower() == source.lower())
        & (master["dataset"].astype(str).str.lower() == dataset.lower())
        & (master["PDB_ID"].astype(str).str.upper().isin(best_out["PDB_ID"].astype(str).str.upper()))
    )
    master2 = master.loc[~mask].copy()

    # Align columns
    for col in master2.columns:
        if col not in best_out.columns:
            best_out[col] = pd.NA
    best_out = best_out[master2.columns]

    return pd.concat([master2, best_out], ignore_index=True)


def update_ligand_rmsd(
    master: pd.DataFrame,
    ligand_df: pd.DataFrame,
    *,
    source: str,
    dataset: str,
    id_col: str,
    id_mode: IdMode,
    rmsd_col: str,
    target_col: str = "RMSD",
) -> pd.DataFrame:
    ligand_df = ligand_df.copy()
    ligand_df["PDB_ID"] = _make_id(ligand_df[id_col], id_mode)
    best = best_row_per_group(ligand_df, "PDB_ID", rmsd_col, keep="min")
    update_map = dict(zip(best["PDB_ID"].astype(str).str.upper(), best[rmsd_col]))

    if target_col not in master.columns:
        master[target_col] = pd.NA

    mask = (
        (master["Source"].astype(str).str.lower() == source.lower())
        & (master["dataset"].astype(str).str.lower() == dataset.lower())
    )

    def _upd(row):
        key = str(row["PDB_ID"]).upper()
        return update_map.get(key, row[target_col])

    master2 = master.copy()
    master2.loc[mask, target_col] = master2.loc[mask].apply(_upd, axis=1)
    return master2


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("upsert_model_metrics")
    p1.add_argument("--master_csv", required=True)
    p1.add_argument("--run_csv", required=True)
    p1.add_argument("--out_csv", required=True)
    p1.add_argument("--source", required=True)
    p1.add_argument("--dataset", required=True)
    p1.add_argument("--id_col", required=True)
    p1.add_argument("--id_mode", required=True)
    p1.add_argument("--rmsd_col", required=True)
    p1.add_argument("--tm_col", default=None)
    p1.add_argument("--plddt_col", default=None)

    p2 = sub.add_parser("update_ligand_rmsd")
    p2.add_argument("--master_csv", required=True)
    p2.add_argument("--ligand_csv", required=True)
    p2.add_argument("--out_csv", required=True)
    p2.add_argument("--source", required=True)
    p2.add_argument("--dataset", required=True)
    p2.add_argument("--id_col", required=True)
    p2.add_argument("--id_mode", required=True)
    p2.add_argument("--rmsd_col", required=True)
    p2.add_argument("--target_col", default="RMSD")

    args = ap.parse_args()

    master = pd.read_csv(args.master_csv)

    if args.cmd == "upsert_model_metrics":
        run_df = pd.read_csv(args.run_csv)
        updated = upsert_model_metrics(
            master,
            run_df,
            source=args.source,
            dataset=args.dataset,
            id_col=args.id_col,
            id_mode=args.id_mode,
            rmsd_col=args.rmsd_col,
            tm_col=args.tm_col,
            plddt_col=args.plddt_col,
        )
        write_csv(updated, args.out_csv)
        print(f"[ok] wrote {args.out_csv}")
    else:
        ligand_df = pd.read_csv(args.ligand_csv)
        updated = update_ligand_rmsd(
            master,
            ligand_df,
            source=args.source,
            dataset=args.dataset,
            id_col=args.id_col,
            id_mode=args.id_mode,
            rmsd_col=args.rmsd_col,
            target_col=args.target_col,
        )
        write_csv(updated, args.out_csv)
        print(f"[ok] wrote {args.out_csv}")


if __name__ == "__main__":
    main()
