# build_dynamicbind_inputs_from_chains_and_csv.py
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chains-dir", required=True)
    p.add_argument("--csv", required=True)
    p.add_argument("--out-dir", default="data/dynamicbind_inputs")
    p.add_argument("--ligand-type", default="Allosteric")
    p.add_argument("--include-empty", action="store_true")
    p.add_argument("--use-smiles-col", default="")
    return p.parse_args()


def infer_pdb_id_from_filename(name: str) -> str:
    stem = Path(name).stem
    if "_" in stem:
        return stem.split("_")[0].upper()
    return stem[:4].upper()


def extract_entries(row: pd.Series, ligand_type: str, smiles_col: str) -> list[str]:
    ligands = str(row["Ligand"]).split()
    ligand_types = str(row["LigandType"]).split()
    smiles_entries = str(row[smiles_col]).split() if smiles_col else []

    out = []
    for idx, (lig, ltyp) in enumerate(zip(ligands, ligand_types)):
        if ltyp != ligand_type:
            continue
        if smiles_col:
            out.append(smiles_entries[idx])
        else:
            out.append(lig.split(":")[0])
    return out


def main() -> int:
    args = parse_args()
    chains_dir = Path(args.chains_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.csv)
    df["PDB"] = df["PDB"].astype(str).str.upper().str.strip()

    wrote, skipped = 0, 0
    for pdb_file in sorted(chains_dir.glob("*.pdb")):
        pdb_id = infer_pdb_id_from_filename(pdb_file.name)
        match = df[df["PDB"] == pdb_id]

        lig_entries: list[str] = []
        if not match.empty:
            lig_entries = extract_entries(match.iloc[0], args.ligand_type, args.use_smiles_col)

        if (not lig_entries) and (not args.include_empty):
            skipped += 1
            continue

        tdir = out_dir / pdb_id
        tdir.mkdir(parents=True, exist_ok=True)

        shutil.copy2(pdb_file, tdir / "protein.pdb")

        if args.use_smiles_col:
            (tdir / "ligand.smi").write_text("\n".join(lig_entries) + "\n")
        else:
            (tdir / "ligand.comp_id.txt").write_text("\n".join(lig_entries) + "\n")

        wrote += 1

    print(f"Done. wrote={wrote} skipped={skipped} out={out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
	