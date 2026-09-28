# Misato docking analysis — scope and provenance

## Current deliverable

This repository is an exploratory fork. The current deliverable is a reproducible analysis of the supplied MISATO **EquiBind** and **DiffDock** outputs, separate from the published ASD/PLA benchmark.

The work is limited to an immutable archive inventory, reconstruction of the evaluated subset and all exclusions, format-specific normalization, validated scoring only after the receptor/reference coordinate frame is established, method coverage and provisional comparison figures, and a concise methods/provenance report.

It does not include new docking runs, TankBind, full MISATO MD reprocessing, edits to legacy outputs, or a manuscript addendum. Those require a later PI decision.

## Provenance constraint

The supplied folders contain ligand-only SDF predictions. They do not contain receptor structures, input preparation scripts, logs, or configurations. Lucas subsequently reported that both docking methods used static RCSB crystal structures and one ligand conformer built from MISATO QM coordinates; the MD trajectories were **not** used for docking. This is a collaborator report, not yet independently verified from his scripts or input files. No protein-aligned docking score may be reported until the actual receptor/reference coordinate frame and ligand mapping are validated. Reproducing these docking inputs does not require the full MD trajectory archive.

The published MD train/validation/test lists are an overlap audit, **not** the denominator for reconstructing Lucas's QM-based filtering. The repository's RCSB retrieval and reference-scoring code can be adapted for an independent crystal-reference evaluation of the delivered poses; this does not require Lucas's complete input archive or per-ID failure table. Crystal-reference ligand pose RMSD and lDDT-PLI are potentially informative after coordinate-frame and ligand-identity validation. Ambiguous references must be flagged or excluded. Pocket/binding-site RMSD and QS score from an unchanged crystal receptor would be trivial or non-discriminative, not comparable to protein cofolding performance claims.

## Decision gate

- **Approved addendum:** freeze source release, manifest, filters, and software configuration, then prepare a reviewed PR with Misato-specific methods and results.
- **Deferred:** retain this fork as the reproducible standalone analysis and do not merge or combine its figures with published claims.
