# MISATO delivered-pose static metrics: steps 3–5

This stage works only from Lucas's 8,419 DiffDock and 8,565 EquiBind
**delivered** ligand poses and the RCSB crystal references. It never redocks,
optimizes a pose, or uses an MD frame. The score base is the versioned
`misato_output/misato_score_recovered_v1.csv`; the original v3 score strings
are not modified. All output files under `misato_output/` are local and
Git-ignored.

## Gnina and confidence

`rescore_gnina_wsl.py` supplies each delivered SDF to official gnina v1.3.3
with `--score_only` and reports **CNNscore** (range 0–1). The receptor is a
temporary PDB made from the crystal asymmetric unit's polymer, without
hydrogens, waters, non-polymer ligands, or graph-verified short-polymer native
ligand copies. It is the same receptor for the two delivered methods for a
target. No source SDF/CIF is edited. A two-pose batch was checked against
separate invocations and returned the same CNNscores; batching does not mix
target receptors. The run is resumable and records per-pose failures rather
than interpreting them as zero.

The installed executable is the official
[v1.3.3 `gnina.cuda12.8.static` release asset](https://github.com/gnina/gnina/releases/tag/v1.3.3),
stored locally as `misato_output/tools/gnina.1.3.3.cuda12.8.static`.
Its local SHA-256 is
`3340c1f49cd3c7c84d8699182a1c6af13c7fa2a22448d1204640446106f72172`.
CUDA 12.8 runtime and cuDNN 9 shared libraries are local under
`misato_output/tools/gnina_cuda_libs/`; the script adds only those paths to
the child gnina process. An older v1.1 CPU smoke was **not** pooled: it uses
a different CNN ensemble and yielded different scores. For reproducibility,
the final table records the gnina version and every row's rescore status.
For throughput, the exact same binary and libraries can be copied to
`/home/ryan/.cache/misato_gnina_runtime/` inside WSL and supplied with
`--gnina` and `--cuda-libs`. A two-pose check yielded identical CNNscores
with the native-filesystem copy. The cache is not a second model version.

The `confidence` column remains **method-native**: DiffDock's raw rank-1
confidence logit is recovered by the existing exact SHA-256 match between
the primary `rank1.sdf` and a scored `rank1_confidence-*.sdf`; EquiBind has
no native confidence and stays blank. The distinct `rescore_confidence`
column holds gnina CNNscore for **both** methods. Do not treat a DiffDock
logit as a probability or compare it numerically with CNNscore.

Gnina's training data include PDBbind-derived complexes, while MISATO was
built from PDBbind 2022. The resulting overlap/similarity can make CNNscore
optimistic on MISATO. It is suitable as a same-tool ranking/calibration
diagnostic, not an independent absolute accuracy claim. Sources:
[gnina repository](https://github.com/gnina/gnina),
[gnina 1.0 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8191141/),
[MISATO study](https://www.nature.com/articles/s43588-024-00627-2).

## Static-receptor QS

`compute_static_qs_wsl.py` adapts the ligand-inclusive contact QS in
`scoring/fix_qs_dynamicbind.py`. For each receptor residue, it calculates
the nearest heavy-atom distance to the assigned native crystal ligand and
to the delivered pose. Contacts use a 12 Å cutoff. Shared contacts receive
`max(0, 1 - |d_native - d_pose|/12)` weight; each contact present in only
one placement contributes one nonshared count. The score is the sum of
shared weights divided by that sum plus nonshared count. Because the
receptor is identical, no superposition or residue remapping is used.
For graph-verified peptide/glycan/composite references, all copies are
available to OST, but QS uses the copy assigned in the formal score.

The completed `misato_output/static_qs_all_v1.csv` has a score for every
formally scored pose: **15,679/15,679**, without failures. Per-row controls
gave exactly 1.0 for native-versus-native and 0.0 after moving the native
ligand far outside the pocket. Spearman correlation with lDDT-PLI is 0.653
overall (0.496 DiffDock; 0.864 EquiBind). A score of 1.0 means the set of
receptor contacts is identical; this QS formula is not a coordinate-placement
metric and may equal 1.0 even if shared contact distances differ. Label it
**static-receptor QS** everywhere; never pool it with co-folding QS, which
also includes receptor prediction error.

## Pocket RMSD

`pocket_rmsd_angstrom` remains blank, and `receptor_mode` is
`static_crystal`. This is **not applicable** to rigid-receptor docking into
the reference crystal structure. A populated 0.0 would suggest a measured
method accuracy and could contaminate cross-method plots. The placement
of the ligand relative to the static pocket is represented separately by
pose RMSD and static-receptor QS.

## Versioned local commands

Run in PowerShell from this checkout. The WSL environment is Ubuntu's
`plb` conda environment. The first command is complete; the gnina command
is the resumable full run. Do not run either on an output file that already
exists without the script's `--resume` flag.

```powershell
wsl -d Ubuntu -- bash -lic 'conda activate plb && cd /mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark-Revision-W-Misato && python scoring/misato/compute_static_qs_wsl.py --workers 8 --out-csv misato_output/static_qs_all_v1.csv'
wsl -d Ubuntu -- bash -lic 'conda activate plb && cd /mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark-Revision-W-Misato && python scoring/misato/rescore_gnina_wsl.py --gnina /home/ryan/.cache/misato_gnina_runtime/gnina --cuda-libs /home/ryan/.cache/misato_gnina_runtime/libs --workers 10 --batch --resume --out-csv misato_output/gnina_rescore_all_v1.csv'
python scoring/misato/merge_static_metrics.py --gnina misato_output/gnina_rescore_all_v1.csv --out-csv misato_output/misato_static_metrics_v1.csv
```

The final merge enforces complete, unique delivered-pose and formally scored
QS key sets. It cannot silently create a partial output while gnina is still
running. The later step-6 manuscript export/figures remain separate work.
