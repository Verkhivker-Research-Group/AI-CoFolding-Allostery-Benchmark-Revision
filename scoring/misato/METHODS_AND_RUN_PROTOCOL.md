# MISATO addendum: cohort definition and reproducible scoring

This is a Methods draft for the *delivered* DiffDock and EquiBind runs, not a
claim that the original QM input list or Lucas's preparation logs were
recovered. `build_methods_accounting.py` regenerates the exact counts below
from the local manifests and crystal audit, retaining reported values in a
separate evidence category. Output under `misato_output/` is local and ignored
by Git; use a new versioned output path for each rerun.
The completed local scoring pass and publication caveats are recorded in
[`METHODS_RESULTS.md`](METHODS_RESULTS.md), which is tracked with the pipeline.

## Which part of MISATO did Lucas actually compute?

**Successful outputs are demonstrably incomplete.** The delivered folders
contain 8,419 DiffDock primary (`rank1.sdf`) poses and 8,565 EquiBind primary
(`lig_equibind_corrected.sdf`) poses across 8,606 distinct IDs. There are
8,379 IDs with both primary poses; 40 have DiffDock only and 186 have
EquiBind only. Four of the 8,569 EquiBind target folders lack the expected
primary SDF. Other DiffDock rank/confidence files are retained in the raw-file
manifest but are not pooled with rank 1. These counts are directly observed.

The public `QM.hdf5` passed Zenodo's MD5 checksum
`1199f0af1eac684da3a6c2ddd0f321df` and contains **19,413** root ID groups,
confirming Lucas's reported QM count. Every one of the 8,606 delivered IDs
occurs in that QM file; 10,807 QM IDs have no delivered target directory.
All 8,419 DiffDock and 8,565 EquiBind primary SDFs also match their QM ID's
heavy-element composition. All of these poses also preserve the QM
heavy-atom **order and untyped bond graph**. This does not establish bond
orders, stereochemistry, protonation, starting coordinates, or which QM
conformer Lucas supplied to the models.
Lucas also reported 19,392 written ligand SDFs and 21 sanitization failures.
The [MISATO paper](https://www.nature.com/articles/s43588-024-00627-2)
states that its **starting PDBbind set** had 19,443 protein–ligand structures,
30 more than Lucas's reported QM count; those are different processing stages,
not necessarily a discrepancy. The delivered primary outputs are 43.4% (DiffDock) and 44.1%
(EquiBind) of the **verified** QM count. These fractions are not verified
per-input success rates: we still have no prepared input manifest, per-ID
status log, or copy of Lucas's protein/SDF preparation code. Thus we can say
exactly which successful poses were delivered, and that they do not cover the
full QM set; we cannot yet say exactly which
of the 19,413 entries were submitted, failed, or never attempted.

Lucas additionally reported 14,057 RCSB proteins saved and 5,356 missing;
14,037 DiffDock input pairs with 8,419 posed and 5,610 RDKit-unreadable
ligand skips (the latter two total 14,029, leaving **eight** input pairs
unexplained); and 8,565 EquiBind posed plus 10,827 failed, which totals the
19,392 reported ligand SDFs. Of the EquiBind failures, he attributed 5,354
to missing proteins, leaving 5,473 without a finer verified status. These
are collaborator-reported counts, not reconstructed from failure logs.

The official MISATO MD split lists contain 16,972 IDs. Of the 8,606 delivered
IDs, 7,482 overlap the split lists and 1,124 do not. The MD splits are **not**
the docking denominator: Lucas reports using static RCSB crystal receptors
and one QM-derived ligand conformer, not MD frames or the MD split selection.
The local `C:\Users\Ryan\Downloads\MD.hdf5` currently has the published byte
length but also has an `.aria2` sidecar, and opening its HDF5 root fails with
`wrong B-tree signature`. The length alone does **not** establish a complete
download; this file must not be used for MD-derived analysis until the
transfer completes and its published checksum/HDF5 structure are verified.
MD data are not required for the static-docking analysis here.

## Primary reference and scoring rule

The scoring input is each delivered primary ligand SDF plus the corresponding
RCSB asymmetric-unit crystal CIF. The native ligand is identified by the
heavy-element signature of the predicted ligand and then checked by an
untyped connectivity/coordinate audit. To avoid a prediction-dependent native
copy choice, the primary cohort requires **one** native crystal ligand
candidate and a successful graph/coordinate check. This yields 3,982
DiffDock poses and 3,963 EquiBind poses individually eligible, including
3,953 paired IDs. These are eligibility counts, **not** final metric counts.
Within those 3,953 pairs, 925 have an exact graph/stereochemistry match
between the two predicted SDFs, 2,979 differ in stereochemical representation,
and 49 require other representation review. The 925 exact-between-methods
pairs define the conservative primary paired comparison; all 3,953 form a
separately labeled topology-matched sensitivity analysis. Neither check by
itself establishes the full stereochemistry of the native crystal ligand.
For all other outputs, the status and reason are preserved. In particular,
multi-copy references selected by DiffDock proximity, composite native
ligands, unmatched native ligands, and graph failures are *not* scored in the
primary comparison. They may be added later as explicitly labeled sensitivity
analyses after independent reference adjudication.

The original `plb_bench` scorer runs in WSL Ubuntu's `plb` conda environment,
parallelized by target. The MISATO adapter removes bound nonpolymer ligands
and waters from the static crystal receptor, adds each predicted ligand SDF,
and reuses the same experimental CIF as reference for both methods. Unlike
the generic benchmark scorer, it passes only the independently selected
native residue to OpenStructure and requires a full-ligand (not substructure)
match. A temporary heavy-atom-only SDF is built for model-complex loading,
because the delivered QM-derived SDFs often contain explicit hydrogens while
the crystal reference ligands do not. The original SDFs are never modified;
heavy-atom coordinates are preserved in the temporary representation.
The adapter records BiSyRMSD and lDDT-PLI separately, with an explicit failure row
if either metric cannot be calculated. No Docker run or oracle-best selection
is used. Pocket Cα RMSD and whole-complex QS are omitted from the primary
comparison because the model receptor is copied from the crystal structure;
they would not independently evaluate docking quality. Pairwise method
comparisons use only IDs where **both** methods have that specific metric.
Formal pose RMSD is also compared with the independent preliminary direct
coordinate RMSD; disagreements greater than 0.5 Å are flagged for review,
not silently removed. A separately labeled direct-coordinate RMSD sensitivity
analysis retains strict-reference pairs whose formal metric fails, so the
potential direction of complete-case bias is visible.

OpenStructure's BiSyRMSD performs a local binding-site superposition and may
map chemically equivalent protein chains. Consequently, it can be small even
when the ligand is far from the selected crystal ligand in the **original
static coordinate frame**. Formal BiSyRMSD/lDDT-PLI and fixed-frame direct
RMSD therefore answer different questions and must be reported separately.
Every >0.5 Å formal-versus-direct RMSD difference is listed for review. No
claim of final fixed-frame docking accuracy should be made from the formal
metric alone. `5EQQ` is explicitly excluded from the WSL run after a
>20-minute scorer stall; its two poses remain in the delivered denominator.
The Python-level runtime alarm is best effort and may not interrupt a long
C++ call.

This is an independent crystal-reference evaluation, not exact reproduction
of Lucas's prepared receptors or ligand conformers. His original receptor
files/scripts and the complete QM/input/failure logs remain the outstanding
requirements for exact input-level provenance.

## Commands

From PowerShell in this repository, rebuild the observed/reported Methods
accounting into a **new** output directory:

```powershell
python scoring/misato/extract_unpaired_signatures.py --out-csv misato_output/inventory/unpaired_pose_signatures_v2.csv
C:\Users\Ryan\AI-CoFolding-Allostery-Benchmark-Revision-W-Misato\misato_output\.venv\Scripts\python.exe scoring/misato/audit_qm_coverage.py --qm C:\Users\Ryan\Downloads\QM.hdf5 --unpaired-signatures misato_output/inventory/unpaired_pose_signatures_v2.csv --out-dir misato_output/qm_coverage_v4
C:\Users\Ryan\AI-CoFolding-Allostery-Benchmark-Revision-W-Misato\misato_output\.venv\Scripts\python.exe scoring/misato/export_qm_heavy_graphs.py --qm C:\Users\Ryan\Downloads\QM.hdf5 --out-jsonl misato_output/inventory/qm_heavy_graphs_v2.jsonl
python scoring/misato/audit_pose_qm_graphs.py --qm-graphs misato_output/inventory/qm_heavy_graphs_v2.jsonl --out-dir misato_output/pose_qm_graph_audit_v2
python scoring/misato/build_methods_accounting.py --qm-summary misato_output/qm_coverage_v4/qm_coverage_summary.json --qm-graph-summary misato_output/pose_qm_graph_audit_v2/pose_qm_graph_summary.json --out-dir misato_output/methods_accounting_v7
```

Run the formal score in the existing local WSL environment (the path below is
this checkout; update `cd` if the checkout moves):

```powershell
wsl -d Ubuntu -- bash -lic 'conda activate plb && cd /mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark-Revision-W-Misato && python scoring/misato/score_static_docking_wsl.py --out-csv misato_output/misato_score_all_heavy_selected_v4.csv --workers 8 --skip-target-id 5EQQ'
```

The score CSV is flushed after every completed target. If interrupted, run
the same command with `--resume` appended; existing target/method rows are
left intact. For a new scoring version, choose a new CSV name. When the full
run is finished, generate the metric-specific denominator and paired tables:

```powershell
python scoring/misato/summarize_static_docking.py --scores misato_output/misato_score_all_heavy_selected_v3.csv --out-dir misato_output/score_summary_heavy_selected_v3
```

The summarizer refuses a partial score file by default. Use
`--allow-partial` only for smoke tests. The generated
`METHODS_ACCOUNTING.md`, `METHODS_RESULTS.md`, JSON summaries, full score CSV,
per-ID paired CSV, and exhaustive RMSD disagreement CSV are the quantitative sources to transfer into the
manuscript after review.

To regenerate the repository's portable delivered-ID/status index from a
completed score file, choose a **new** output filename and review it before
replacing the currently tracked snapshot:

```powershell
python scoring/misato/export_delivered_cohort.py --scores misato_output/misato_score_all_heavy_selected_v3.csv --out-csv misato_output/delivered_cohort_status_regenerated.csv
```

The tracked [`delivered_cohort_status.csv`](delivered_cohort_status.csv)
contains no local paths or ligand coordinates. It identifies exactly which
QM IDs received each primary pose and whether the corresponding static score
was obtained; it does not imply that every other QM ID was attempted.
