# Input-generation scripts

This folder contains scripts that generate model input files from reference structures and dataset metadata.

Notebook converted:
- `cif.ipynb` → `build_af3_json_from_chains_and_csv.py`

## Build AF3 request JSONs (protein chain + allosteric ligand IDs)

Requirements:
- A folder of single-chain PDB files, e.g. `content/Chains/1CKJ_A.pdb`
- A CSV with columns: `PDB`, `Ligand`, `LigandType`
  - `Ligand` contains whitespace-separated entries like `Y3I:...
  - `LigandType` is whitespace-separated entries aligned with `Ligand` entries.

Command:
```bash
python scripts/1_inputs/build_af3_json_from_chains_and_csv.py \
  --chains-dir content/Chains \
  --csv output.csv \
  --out-dir data/af3_requests \
  --only-allosteric
```

Output:
- `data/af3_requests/<PDB>.json` for each PDB ID found.
