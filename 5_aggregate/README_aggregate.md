# Aggregation scripts (replacing `conc.ipynb`)

This folder contains the modular replacements for the original ad-hoc notebook.

Typical flow:

1. (Optional) Select best rows from raw model outputs

```bash
python scripts/5_aggregate/01_select_best_rows.py af3_best \
  --in_csv path/to/af3/alloresults.csv \
  --out_csv spread/af3allo_best.csv \
  --seed_col Complex_Seed --rmsd_col RMSD
```

2. Combine a folder of per-source CSVs into a single file

```bash
python scripts/5_aggregate/02_folder_tools.py combine \
  --folder spread --pattern "*.csv" \
  --out_csv combined_spreadsheet.csv
```

3. Normalize to canonical master schema

```bash
python scripts/5_aggregate/03_normalize_master_schema.py \
  --in_csv combined_spreadsheet.csv \
  --out_csv master.csv
```

4. Upsert model-level metrics into the master (repeat per model/dataset)

```bash
python scripts/5_aggregate/04_update_master_from_runs.py upsert_model_metrics \
  --master_csv master.csv \
  --run_csv path/to/chai/chaimain.csv \
  --out_csv master_plus_chai.csv \
  --source chai --dataset main \
  --id_col pdb_id --id_mode pdb_id_first4 \
  --rmsd_col rmsd --tm_col tm_score --plddt_col mean_plddt
```

5. Update ligand RMSD measurements into the master

```bash
python scripts/5_aggregate/04_update_master_from_runs.py update_ligand_rmsd \
  --master_csv master_plus_chai.csv \
  --ligand_csv path/to/af3/alloresults.csv \
  --out_csv master_plus_chai_plus_ligrmsd.csv \
  --source af3 --dataset allo \
  --id_col Complex_Seed --id_mode complex_seed_first5 \
  --rmsd_col RMSD --target_col RMSD
```

6. Finalize schema for paper/figures

```bash
python scripts/5_aggregate/05_finalize_master.py \
  --in_csv master_plus_chai_plus_ligrmsd.csv \
  --out_csv cleanmaster.csv \
  --rename_to_final --dropna
```
