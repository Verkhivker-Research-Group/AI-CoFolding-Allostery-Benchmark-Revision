# plb_bench — Protein-Ligand Benchmarking Pipeline

High-throughput evaluation of protein-ligand pose predictions from
**AlphaFold 3**, **Protenix**, **Chai-1**, and **DynamicBind**, scored with
**OpenStructure** (OST).

---

## Architecture

```
           input root/                                   ┌─ Phase 1: crawl + parse + rank
           ├── 1abc/                                     │  (registry.py)
           │   ├── af3/seed_*/model.cif + *.json ────────┤
           │   ├── proteinx/sample_*.cif + *.json ───────┤
           │   ├── chai/pred.model_idx_*.cif + *.npz ────┤
           │   └── dynamicbind/rank*_receptor*.pdb + sdf ┘
           └── 2def/...
                       │
                       ▼
          ┌───────────────────────────────┐
          │   ModelRecord[] per producer  │  unified schema
          │   sorted by confidence,       │  (schema.py)
          │   model_idx 0..N-1            │
          └───────────────┬───────────────┘
                          │
                          ▼
      ┌──────────────── Phase 2 ─────────────────┐
      │  Robust CIF Sanitizer (in-memory)         │
      │  gemmi → setup_entities() ───────┐        │
      │     ↓ if empty structure         │        │
      │  manual _atom_site rebuild ──────┼──→ clean mmCIF text
      │     ↓ if that fails              │        │
      │  biopython MMCIFIO ──────────────┘        │
      │                                            │
      │  + merges SDF ligand as HETATM chain (RDKit)
      └──────────────────────┬────────────────────┘
                             │
                             ▼
           ┌── Phase 3: reference ──┐
           │  local cache (atomic)  │  with per-PDB lock
           │    ↓ miss              │
           │  RCSB assembly1.cif.gz │
           │    ↓ miss              │
           │  RCSB <pdb>.cif.gz     │
           └───────────┬────────────┘
                       │
                       ▼
         ┌── Phase 4: OpenStructure ──┐
         │  SCRMSDScorer  → BiSyRMSD  │
         │  LDDTPLIScorer → lDDT-PLI  │
         │  qsscore.QSScorer → QS_global
         └───────────┬────────────────┘
                     │
                     ▼
         flat Parquet / CSV  (benchmark.parquet)
         one row per (pdb_id, producer, model_idx)
```

### Scheduling

Work is grouped by **(pdb_id, producer)** and distributed across a
`ProcessPoolExecutor` (or Ray via `--ray`). One group processes 2–5 models
together, amortizing the reference-loading cost.

### Robustness

Every failure mode emits a row rather than crashing the worker:

| `status` column        | Meaning                                      |
|------------------------|----------------------------------------------|
| `ok`                   | All three metrics computed                   |
| `ost_no_metrics`       | OST loaded but all three scorers returned None|
| `sanitize_failed`      | The CIF/PDB could not be rebuilt             |
| `no_reference`         | Reference unavailable (local + RCSB both)    |
| `ost_failed`           | Scorer raised unexpectedly                   |
| `worker_crashed`       | Worker process died (pool-level catch)       |

---

## Install

```bash
pip install -e .
# OpenStructure must be installed separately:
#   conda install -c bioconda openstructure
#   (or use the scicore docker image)
```

## Usage

```bash
plb-bench run \
    --input /data/predictions \
    --refs  /data/refs_cache \
    --out   /data/benchmark_out \
    --workers 32
```

Programmatic:

```python
from plb_bench.pipeline import run, PipelineConfig

run(PipelineConfig(
    input_root="/data/predictions",
    refs_dir="/data/refs_cache",
    output_dir="/data/benchmark_out",
    producers=["af3", "chai"],      # subset
    n_workers=32,
    output_format="parquet",
))
```

## Adding a new producer

```python
from plb_bench.registry import register_parser
from plb_bench.schema import ModelRecord

@register_parser("myprod")
def parse_myprod(pdb_dir, pdb_id):
    # return list[ModelRecord] with .confidence populated
    ...
```

No other code changes needed — the orchestrator discovers it automatically.

## Output schema

One row per scored model, optimised for pandas / polars / duckdb:

```
pdb_id, producer, model_idx, confidence,
bisy_rmsd, lddt_pli, qs_global,
status, error, reference_source,
raw_ranking_score, raw_iptm, raw_ptm,           # AF3
raw_ranking_confidence, raw_complex_plddt,      # Protenix
raw_aggregate_score,                            # Chai
raw_confidence, raw_affinity,                   # DynamicBind
...
```
