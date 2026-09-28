# MISATO input-provenance findings

Generated from the delivered output folders and the public MISATO source code.
This is an evidence log, not a statement of final benchmark results.

## Lucas's reported input preparation (pending file-level verification)

Lucas reported that he downloaded the complete `MD.hdf5` and `QM.hdf5`, but
used **neither MD trajectories nor MD-derived receptor frames for docking**.
Both DiffDock and EquiBind received static RCSB crystal protein structures and
one ligand conformer built from MISATO QM coordinates. He did not use the
official MD train/validation/test splits, MD restart/topology files, or
electronic densities. These statements are collaborator-reported provenance;
his preparation scripts and exact input files have not yet been inspected.

Lucas's counts reconcile as follows, with one unresolved discrepancy:

- QM: 19,413 ligand IDs in his HDF5, versus 19,443 stated in the paper.
- Ligand SDFs: 19,392 written + 21 sanitization failures = 19,413.
- RCSB proteins: 14,057 saved + 5,356 missing = 19,413. Lucas says most
  missing proteins reflect failed downloads, not a scientific cutoff.
- DiffDock: 14,037 input pairs; 8,419 posed + 5,610 RDKit-unreadable-ligand
  skips = 14,029, leaving **8 input pairs unaccounted for** in that report.
- EquiBind: 8,565 posed + 10,827 failed = 19,392 ligand SDFs. Lucas attributes
  5,354 of the failures to missing proteins; the remainder need finer labels.

The 8,419 DiffDock and 8,565 EquiBind poses match the delivered primary-file
counts below. No per-ID failure log or input manifest has yet been supplied to
verify the other totals or the reason for each exclusion.

## Delivered-output audit

- EquiBind: 8,569 target directories; 8,565 contain the expected
  `lig_equibind_corrected.sdf`; four are empty (`1A37`, `2Z4W`, `3MBZ`,
  `4GD6`).
- DiffDock: 8,419 target directories; each has a canonical `rank1.sdf`.
  A further 84,936 rank/confidence SDFs are retained as unselected
  candidate/repeat-run evidence.
- Shared primary cohort: 8,379 PDB IDs after the three shared EquiBind
  absences are excluded with an explicit reason.
- All 8,379 paired primary ligands have the same heavy-element composition.
  8,261 have matching connectivity; 1,712 match exactly and 6,549 differ in
  stereochemical representation. The 118 remaining pairs parse as coordinate
  SDFs but do not yield an RDKit hydrogen-normalized graph descriptor.

The generated local manifests in `misato_output/inventory/` contain the full
per-ID disposition and SHA-256 record for every delivered SDF.

## Comparison with the published MISATO MD splits

The checksum-verified `train_MD.txt`, `val_MD.txt`, and `test_MD.txt` lists
contain 16,972 distinct IDs (13,765 train; 1,595 validation; 1,612 test).
The delivered output folders contain 8,606 distinct IDs. Of those, 7,482
appear in the published MD splits and 1,124 do not. Among the split-listed
IDs, 7,368 have both primary ligand candidates; 9,490 split-listed IDs were
not delivered by either method. Of the 1,124 out-of-split IDs, 1,011 have
both primary ligand candidates.

These are ID-level **MD-split overlap** observations, not a reconstruction of
Lucas's filtering denominator: Lucas says he used QM-derived ligands and did
not apply these splits. An out-of-split ID is not proof that it is absent from
QM.hdf5 or every MISATO source file. See
`scoring/misato/compare_published_splits.py` and
`misato_output/split_comparison/` for the reproducible comparison.

## Reference audit: why PDB alone is insufficient

The first 100 target IDs were resolved against RCSB biological assemblies:

- 68 had one single-component non-polymer candidate matching the predicted
  ligand's heavy-element signature;
- 2 had one composite candidate;
- 6 were ambiguous; and
- 24 had no matching component combination of up to three residue types.

For example, `11GS` has a 39-heavy-atom predicted ligand. Its experimental
structure contains separate 20-heavy-atom `GSH` and 19-heavy-atom `EAA`
components; their combined signature matches the prediction. This shows that
the results may represent a ligand assembly rather than one crystallographic
residue.

## What the original MISATO data provides

The public MISATO processing code represents every PDB ID as an HDF5 group
with `trajectory_coordinates`, `atoms_type`, `atoms_number`, `atoms_residue`,
`atoms_element`, and `molecules_begin_atom_index`. Its preprocessing selects
ligand atoms where `atoms_residue == 0`, falling back to the final molecule for
peptide ligands. This selection can naturally produce a multi-component ligand
assembly. See [preprocessing_db.py](https://github.com/t7morgen/misato-dataset/blob/master/src/data/processing/preprocessing_db.py).

This was verified against the repository's supplied `tiny_md.hdf5` sample.
For `10GS`, `11GS`, and `16PK`, the `atoms_residue == 0` atom selection has
the exact same heavy-element signature as the delivered EquiBind SDF. The
`11GS` selection contains 66 total atoms / 39 heavy atoms, consistent with the
combined `GSH` + `EAA` representation above. This chemical-identity match is
not evidence that an MD frame was used; Lucas reports that none was.

## Path to crystal-reference scoring

The repository already contains RCSB download helpers, experimental-reference
retrieval, and ligand-scoring machinery. These can be adapted to the delivered
Misato output IDs without waiting for a per-ID failure table from Lucas. The
reference ligand can be inferred by atom/graph matching where unique;
ambiguous, composite, or unmatched cases must be flagged or reviewed rather
than silently assigned. The old `select_ligand.py` default of choosing the
largest non-water residue is not sufficient for cases such as `11GS`.

For each scored target, first verify that the predicted SDF coordinates and
the chosen RCSB crystal reference share a valid coordinate frame, and that
the selected native ligand corresponds to the predicted chemical entity.
Then crystal-reference ligand pose RMSD and lDDT-PLI can be reported with an
explicit scored/excluded denominator. QS and protein pocket RMSD from an
unchanged crystal receptor would be trivial, not an independent docking result.

Lucas's exact prepared receptor/SDF files or preparation commands are useful
for verifying exact reproduction, especially if coordinate frames do not
match. They are **not a prerequisite** to starting an independent,
crystal-reference evaluation of the delivered poses. A complete per-ID failure
table and an explanation of eight DiffDock input pairs are not required for
that narrower evaluation. The full MD archive and its topology/restart files
are also not required.
