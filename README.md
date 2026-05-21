# Benchmarking AI-Based Co-Folding and Docking Models for Predicting Structures of Orthosteric and Allosteric Ligand–Protein Complexes  
## Decoding the Allosteric Blind Spot Using a Landscape-Guided Interpretable AI Framework

## Pipeline Overview

![Benchmarking pipeline diagram](/flowchart.png)


## Models Evaluated

### Allosteric Dataset (ASD)
- AlphaFold 3 (AF3)
- Chai-1
- Boltz-2
- Protenix

### Protein–Ligand Allosteric Dataset (PLA)
- AlphaFold 3 (AF3)
- Chai-1
- Protenix


## Datasets

### 1. Allosteric Dataset (ASD)
- Designed to interrogate the allosteric blind spot in protein–ligand structure prediction
- Ligands bind to conformationally dynamic pockets spatially remote from canonical active sites
- **3,680 unique compounds** across four models (af3, boltz, chai, proteinx)
- Compound keys follow `{pdb4}_{ligand}` format (e.g., `3hrf_p47`)

**Components:**
- **Kinase-Centric Subset (KinCoRe):**
  - Obtained from Dunbrack Lab's Kinase Conformation Resource (KinCoRe)
  - 217 complexes of Type III (back-pocket) and Type IV (distal myristoyl-site) kinase inhibitors

- **Curated Allosteric Kinase Set:**
  - 136 experimentally validated allosteric kinase inhibitors
  - Confirmed via X-ray crystallography (Hu et al., 2021)

- **Diversified Allosteric Subset (Non-Kinase):**
  - 1,613 structures spanning non-kinase protein families
  - 1,277 unique ligands
  - Derived from a large-scale ligand binding pocket dataset for drug discovery (Moine-Franel et al., 2024)
  - Provides a statistically robust probe into the energetic heterogeneity of allostery

### 2. Protein–Ligand Allosteric Dataset (PLA)
- Focused evaluation of allosteric ligand binding across a curated set of protein–ligand complexes
- **652 unique compounds** across three models (af3_pla, chai_pla, protenix_pla)
- Compound keys follow `{uniprot}_{ligand}` format (e.g., `o14936_v5w`) for cross-model comparison
- Multiple PDB structures per UniProt+ligand combination are normal; best model selected by lowest pose RMSD


## Evaluation Metrics

### Ligand-Level Accuracy
- **Ligand Pose RMSD (BiSyRMSD)**
  - Symmetry-corrected RMSD between predicted and reference ligand heavy atoms
  - Computed via global protein superposition followed by element-group permutation matching
  - Handles multiple ligand copies in crystal structures; selects the best-matching (ref, model) pair
  - Primary implementation uses OpenStructure's SCRMSDScorer; fallback to custom superposition-based
    scorer for models that use non-standard ligand residue names (see `PIPELINE_FIXES.txt`)

### Pocket and Interaction Fidelity
- **Pocket RMSD**
  - Binding-site Cα RMSD after chain-aware Kabsch superposition on pocket residues
  - Pocket defined as reference Cα atoms within 5 Å of the correct reference ligand instance
  - Strategy A: exact (chain, residue number) pairing — used when numbering is consistent
  - Strategy B: ligand-anchored nearest-neighbour matching in ligand-centred frame — fallback
    for chain-merged models, multimers, or non-standard chain naming
- **QS-score** — Quaternary structure reproduction score (OpenStructure)
- **lDDT-PLI** — Local Distance Difference Test for protein–ligand interface accuracy (OpenStructure)


## Output Files

Results are in `evalspreadsheets/main/` (ASD) and `evalspreadsheets/pla/` (PLA).  
Two file variants are provided per producer:

| File pattern | Description | Use for |
|---|---|---|
| `{producer}mainligand_best.csv` | Best model per compound (lowest pose RMSD), all metrics | Analysis, tables |
| `{producer}mainligand_avg.csv` | Mean metrics across all 5 models per compound | Analysis, tables |
| `{producer}mainligand_clean_best.csv` | Best, with `pocket_rmsd > 20 Å` or `bisy_rmsd > 20 Å` rows removed | **Figures** |
| `{producer}mainligand_clean_avg.csv` | Avg, same outlier filter | **Figures** |

Columns: `id` (compound key), `pose rmsd`, `pocket rmsd`, `qs score`, `lddt-pli`, `confidence`

> **Use the `_clean_` files for all figures and statistical comparisons.**  
> The clean files remove a small number of residual computation artifacts (< 0.5% of rows per producer)
> that produce non-physical distributions. The regular files retain full data including NaN rows
> (genuine cases where the model did not place the ligand).


## Data Quality Summary

After all pipeline fixes, the final NaN rates for successfully scored compounds are:

| Producer | Pose RMSD NaN | Pocket RMSD NaN | Artifact rows excluded from clean |
|---|---|---|---|
| af3 | 0.4% | 1.5% | 24 (0.4%) |
| boltz | 0.0% | 0.9% | 3 (0.3%) |
| chai | 0.0% | 0.0% | 0 |
| proteinx | 3.3% | 0.5% | 0 |
| af3_pla | 3.9% | 3.9% | 4 (0.4%) |
| chai_pla | 1.1% | 1.1% | 8 (0.2%) |
| protenix_pla | 3.9% | 3.9% | 0 |

All remaining NaN values represent models where the ligand was genuinely absent from the prediction.

> See `PIPELINE_FIXES.txt` for a full description of the 10 bugs identified and fixed in the
> plb_bench scoring pipeline after the initial run.


## Repository Structure

The repository is organized as a modular, stage-based analysis pipeline. Each folder corresponds to a distinct step in the benchmarking workflow, from dataset preparation to evaluation and aggregation.

> **Note:** Raw model outputs, normalized CIF trees, and reference structures (ground truths) are stored in a separate local archive (`AI-CoFolding-Archive/`) and are not included in this repository. The pipeline notebooks reference those paths directly. Reference CIFs are also downloaded automatically during the pre-warm step if not present locally.

### Folder Descriptions

- **`0_setup/`**  
  Scripts for downloading structures, constructing orthosteric and allosteric datasets, and validating dataset integrity.

- **`1_inputs/`**  
  Construction of model-specific inputs, including FASTA files, ligand specifications, and AlphaFold 3 JSON requests.

- **`2_run_models/`**  
  Collection and standardization of prediction outputs from AF3, Chai-1, Boltz-2, and Protenix into a unified format.

- **`3_postprocess_predictions/`**  
  Post-processing of CIF structures, including ligand extraction, pocket and binding-site residue definition, and structure normalization.

- **`4_score/`**  
  Implementation of all evaluation metrics: ligand pose RMSD (BiSyRMSD), pocket RMSD, QS-score, and lDDT-PLI via OpenStructure.

- **`5_aggregate/`**  
  Aggregation of raw metric outputs into unified master tables and summary statistics used for analysis.

- **`scoring/`**  
  End-to-end pipeline notebooks and the `plb_bench` scoring package:
  - `plb_bench_run.ipynb` — full ASD pipeline (normalize → score → pocket RMSD → export)
  - `plb_bench_run_pla.ipynb` — full PLA pipeline (normalize → score → pocket RMSD → export)
  - `plb_bench/` — internal scoring library (discoverers, normalizer, OST scorer, schema)

- **`evalspreadsheets/`**  
  Final benchmark CSVs produced by the pipeline:
  - `main/` — ASD results per producer, regular and `_clean_` variants
  - `pla/` — PLA results per producer, regular and `_clean_` variants
  - Columns: `id`, `pose rmsd`, `pocket rmsd`, `qs score`, `lddt-pli`, `confidence`

- **`references/`**  
  Small reference files included in the repo:
  - `pla_uniprot_map.csv` — UniProt → PDB mapping for PLA compounds
  - Reference CIFs (ground truth structures) are stored in `AI-CoFolding-Archive/references_ref_cifs/` (not in this repo)

- **`Figures/`**  
  Publication-ready figures generated from aggregated results.

- **`utils/`**  
  Shared helper functions and utilities used across multiple pipeline stages.

- **`PIPELINE_FIXES.txt`**  
  Full documentation of 10 bugs identified and corrected in the plb_bench scoring pipeline,
  including pocket RMSD computation errors, id collision bugs, BiSyRMSD NaN recovery, and
  reference CIF lookup failures.

- **`PLA_PIPELINE_PROCEDURE.txt`**  
  Step-by-step procedure for running the PLA pipeline from raw model outputs through to
  final CSVs, including data copying, UniProt→PDB mapping, normalization, scoring, and export.


## Running the Pipeline

### Requirements
- **WSL (Ubuntu)** with `conda activate plb`
- **OpenStructure** installed in the `plb` conda environment
- Raw model outputs and reference CIFs in `AI-CoFolding-Archive/`

### ASD Pipeline
```bash
conda activate plb
cd /mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark/scoring
jupyter nbconvert --to notebook --execute --inplace plb_bench_run.ipynb \
  --ExecutePreprocessor.kernel_name=plb --ExecutePreprocessor.timeout=21600
```

### PLA Pipeline
```bash
conda activate plb
cd /mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark/scoring
jupyter nbconvert --to notebook --execute --inplace plb_bench_run_pla.ipynb \
  --ExecutePreprocessor.kernel_name=plb --ExecutePreprocessor.timeout=21600
```

Output CSVs are written to `evalspreadsheets/main/` (ASD) and `evalspreadsheets/pla/` (PLA).

### Post-Processing Fixes (apply after initial pipeline run)
```bash
conda activate plb
# Fix id collisions, BiSyRMSD NaN, and re-export all CSVs with clean variants
python /mnt/c/Temp/fix_all_issues.py 2>&1 | tee /mnt/c/Temp/fix_all_issues.log

# Fix proteinx ASD pocket RMSD (assembly CIF lookup)
python /mnt/c/Temp/fix_proteinx_pocket_rmsd.py 2>&1 | tee /mnt/c/Temp/fix_proteinx_pocket_rmsd.log
```


## Dependencies

### Python
- **Python 3.10+**

### Core Python Packages
- `pandas`
- `numpy`
- `scipy`
- `pyyaml`
- `biopython`
- `pyarrow`

### OpenStructure (Required for Scoring)
This project relies on **OpenStructure (OST)** and its Python bindings for all structure scoring:
- `ost.mol.alg.ligand_scoring_scrmsd` (BiSyRMSD / SCRMSDScorer)
- `ost.mol.alg.ligand_scoring_lddtpli` (lDDT-PLI)
- `ost.mol.alg.qsscore` (QS-score)

> OpenStructure is Linux-only and must be installed in WSL. All scoring steps run inside WSL via `jupyter nbconvert`.

### Standard Library
- `argparse`, `pathlib`, `json`, `csv`, `re`, `shutil`, `dataclasses`, `typing`


## References

1. **Hu H, Laufkötter O, Miljković F, Bajorath J.**  
   *Data set of competitive and allosteric protein kinase inhibitors confirmed by X-ray crystallography.*  
   Data Brief, 2021, **35**:106816.  
   https://doi.org/10.1016/j.dib.2021.106816

2. **Moine-Franel A, Mareuil F, Nilges M, Ciambur CB, Sperandio O.**  
   *A comprehensive dataset of protein–protein interactions and ligand binding pockets for advancing drug discovery.*  
   Scientific Data, 2024, **11**(1):402.  
   https://doi.org/10.1038/s41597-024-03233-z
