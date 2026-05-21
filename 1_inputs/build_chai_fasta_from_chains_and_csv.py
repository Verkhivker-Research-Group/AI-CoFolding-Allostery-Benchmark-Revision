# build_chai_fasta_from_chains_and_csv.py
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chains-dir", required=True)
    p.add_argument("--csv", required=True)
    p.add_argument("--out-dir", default="data/chai_requests")
    p.add_argument("--ligand-type", default="Allosteric")
    p.add_argument("--include-empty", action="store_true")
    return p.parse_args()


def pdb_to_sequence(pdb_path: Path) -> str:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    ppb = PPBuilder()
    seq = ""
    for model in structure:
        for chain in model:
            for pep in ppb.build_peptides(chain):
                seq += str(pep.get_sequence())
    return seq


def infer_pdb_id_from_filename(name: str) -> str:
    stem = Path(name).stem
    if "_" in stem:
        return stem.split("_")[0].upper()
    return stem[:4].upper()


def extract_ligands(row: pd.Series, ligand_type: str) -> list[str]:
    ligands = str(row["Ligand"]).split()
    ligand_types = str(row["LigandType"]).split()
    out = []
    for lig, ltyp in zip(ligands, ligand_types):
        if ltyp == ligand_type:
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
        seq = pdb_to_sequence(pdb_file)

        match = df[df["PDB"] == pdb_id]
        ligs: list[str] = []
        if not match.empty:
            ligs = extract_ligands(match.iloc[0], args.ligand_type)

        if (not ligs) and (not args.include_empty):
            skipped += 1
            continue

        lines = [f">{pdb_id}|protein|A", seq]
        for i, comp_id in enumerate(ligs, start=1):
            lines.append(f">{pdb_id}|ligand|{i}|comp_id={comp_id}")
            lines.append(comp_id)

        out_path = out_dir / f"{pdb_id}.fasta"
        out_path.write_text("\n".join(lines) + "\n")
        wrote += 1

    print(f"Done. wrote={wrote} skipped={skipped} out={out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
	