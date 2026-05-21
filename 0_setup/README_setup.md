# Setup scripts

These scripts are scriptified versions of the exploratory notebooks:
- `pdb.ipynb` → `download_pdbs_from_excel.py`
- `jsonprot.ipynb` → `download_cifs_from_excel.py`

## Download PDBs
```bash
python scripts/0_setup/download_pdbs_from_excel.py --excel data.xlsx --out-dir data/pdb_files
```

## Download CIFs
```bash
python scripts/0_setup/download_cifs_from_excel.py --excel data.xlsx --out-dir data/cif_files
```

## Notes
- Default column name is `entryName` and the PDB ID is extracted as the substring before the first underscore.
- Use `--limit` for quick smoke tests.
