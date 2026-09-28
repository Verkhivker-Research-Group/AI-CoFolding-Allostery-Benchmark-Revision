# MISATO delivered-run analysis: Methods and results record

This record applies to **Lucas's delivered DiffDock rank-1 and EquiBind
corrected ligand poses**, not to every attempted run or to the complete
MISATO MD data. It was checked against
`misato_output/methods_accounting_v7/` and
`misato_output/score_summary_heavy_selected_v2/` on 25 September 2026.
Those local, Git-ignored directories contain the machine-readable per-ID
accounting, full scores, paired table, and RMSD quality-control list.
`METHODS_AND_RUN_PROTOCOL.md` defines the reference and scoring rules.

## Dataset and run completion

The checksum-verified public `QM.hdf5` contains **19,413 ID groups**. Every
delivered target ID occurs in it. The delivered folders contain 8,419
DiffDock primary poses and 8,565 EquiBind primary poses, spanning 8,606
distinct IDs: 8,379 with both poses, 40 DiffDock-only, and 186
EquiBind-only. Thus only 43.4% and 44.1% of the QM IDs, respectively, have
a delivered primary pose. The exact successful IDs and each method's score
status are in the tracked, path-free
[`delivered_cohort_status.csv`](delivered_cohort_status.csv); the richer
local inventory remains under `misato_output/inventory/`. Neither method's
successful outputs cover the entire QM set. All 16,984 primary poses match their corresponding QM
ligand's heavy-element composition, atom order, and untyped connectivity.
This does not verify their bond orders, stereochemistry, protonation, or
Lucas's starting conformer.

Lucas reported preparing 19,392 ligand SDFs from those 19,413 QM IDs (21
sanitization failures) and saving 14,057 RCSB proteins (5,356 missing).
He reported 14,037 DiffDock input pairs, 8,419 poses and 5,610
RDKit-unreadable SDF skips; these last two counts leave **eight input pairs
unexplained**. He reported 8,565 EquiBind poses and 10,827 failures, of
which 5,354 were attributed to missing proteins and 5,473 have no finer
verified status. These preparation, submission, and failure counts are
**collaborator-reported**, not reconstructed from input files or logs.
Consequently, the delivered successful subset is known exactly, but the
complete attempted cohort and reasons for all absent outputs are not.
The MISATO paper's 19,443 starting PDBbind structures are a different
processing stage from the 19,413 verified QM groups. The 16,972 official
MD-split IDs are not a docking denominator; 7,482 delivered IDs overlap
them, while 1,124 do not. Lucas reports using static RCSB receptors and one
QM-derived ligand conformer rather than MD frames.

## Crystal reference, scoring, and denominators

Only an unambiguous, single native ligand with a successful independent
heavy-atom graph and coordinate audit is eligible for the conservative
comparison. This gives 3,982 individually eligible DiffDock poses, 3,963
EquiBind poses, and **3,953 paired IDs** before OpenStructure scoring.
For 925 of those pairs the two predicted SDFs have the same exact ligand
graph/stereochemical representation; 2,979 pairs differ in stereochemical
representation and 49 require other representation review. The 925 form the
conservative paired subset, not a proof of matching *native-crystal*
stereochemistry. The 3,953 form a broader sensitivity cohort. Multi-copy
native ligands selected using a prediction, composite references, unmatched
ligands, and failed graph checks are excluded from these comparisons.

The reference-selection audit partitions each method's delivered poses as
follows. Selection categories are mutually exclusive; the later graph check
reduces the single-candidate counts to the strict eligible counts above.

| Reference-selection category | DiffDock | EquiBind |
| --- | ---: | ---: |
| Single candidate before graph check | 4,007 | 3,990 |
| Copy chosen by DiffDock proximity (excluded) | 2,753 | 2,743 |
| Copy chosen by EquiBind proximity (excluded) | 0 | 61 |
| Ambiguous crystal copies | 110 | 201 |
| Composite reference requiring review | 131 | 131 |
| No native heavy-element match | 1,418 | 1,439 |

The 25 DiffDock and 27 EquiBind single-candidate poses that did not enter the
strict cohort failed the independent graph/coordinate check. The full audit
retains each target's selection and failure reason rather than treating all
excluded poses as equivalent.

The existing `plb_bench` code was run locally under WSL Ubuntu with Python
3.11.6, OpenStructure 2.11.1, RDKit 2026.03.1, and Gemmi 0.7.5. The
protein coordinates came from the RCSB crystal reference; nonpolymer
ligands and waters were removed before adding each predicted SDF. The
reference passed to OpenStructure was restricted to the independently
selected native residue. Explicit hydrogens were removed **only from a
temporary scoring copy** of the SDF; original poses and heavy-atom
coordinates were preserved. Full-ligand matching, rather than
substructure-only matching or best-of-ranks selection, was required.
BiSyRMSD and lDDT-PLI were calculated separately. Missing values were
retained as missing; the comparison denominator is the number of IDs on
which **both** methods have that metric. Pocket Cα RMSD and whole-complex
QS were not compared because the model receptor was copied from the
crystal. The pathological target `5EQQ` was predeclared as a runtime
exclusion after a prior >20-minute OpenStructure stall; both poses remain
in the delivered and eligible denominators.

| Cohort / outcome | DiffDock | EquiBind |
| --- | ---: | ---: |
| Delivered primary poses | 8,419 | 8,565 |
| Strict-reference eligible | 3,982 | 3,963 |
| Reference/graph excluded | 4,437 | 4,602 |
| Runtime excluded (`5EQQ`) | 1 | 1 |
| Scoring exception | 1 | 1 |
| No BiSyRMSD (or only partial metrics) | 391 | 917 |
| BiSyRMSD obtained | 3,589 | 3,044 |
| lDDT-PLI obtained | 3,589 | 3,047 |

The two scoring exceptions are the same target, `4FBX`, where the selected
native residue was not recovered by the OpenStructure residue/element
check. Of the 1,308 partial-metric rows, 1,305 have neither metric and
three EquiBind rows have lDDT-PLI only. No failed score was filled with
zero or silently counted as success.

| Paired analysis | Eligible IDs | Both BiSyRMSD | Both lDDT-PLI | DiffDock median BiSyRMSD (Å) | EquiBind median BiSyRMSD (Å) | DiffDock median lDDT-PLI | EquiBind median lDDT-PLI |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Exact prediction graph / strict reference | 925 | 749 | 749 | 0.573 | 3.042 | 0.959 | 0.532 |
| Broader topology-matched / strict reference | 3,953 | 2,926 | 2,928 | 0.887 | 3.417 | 0.906 | 0.488 |

In the exact-graph cohort, DiffDock has lower BiSyRMSD in 696 of the 749
scored pairs, and higher lDDT-PLI in 710 of 749 (three lDDT ties).
The corresponding 2 Å BiSyRMSD counts are 662/749 for DiffDock and
178/749 for EquiBind. These are comparisons **among scoreable, paired,
delivered poses**, not dataset-wide success rates.

## Coordinate-frame quality control and interpretation

OpenStructure BiSyRMSD locally superposes the binding site and can map
equivalent protein chains. It therefore is **not** identical to the
fixed-frame pose displacement in the original crystal coordinates. Against
the independently computed provisional, untyped-graph direct RMSD, 171
scored method rows differ by >0.5 Å; 152 differ by >2 Å and 94 by >10 Å.
The full list is
`misato_output/score_summary_heavy_selected_v2/rmsd_qc_disagreements.csv`.
The differences are retained and flagged, not silently removed. For
example, target `4WSK` DiffDock has BiSyRMSD 1.31 Å versus fixed-frame
direct RMSD 113.23 Å; it cannot be described as accurate placement in
the original crystal frame solely from BiSyRMSD.

For all 925 exact-graph/strict-reference pairs, the provisional fixed-frame
direct RMSD median is 0.639 Å for DiffDock and 3.181 Å for EquiBind; these
are **exploratory sensitivity results**, not the original OpenStructure
metric. A post hoc subset of 717 exact-graph pairs has both formal scores
and ≤0.5 Å formal-versus-direct disagreement for both methods. It is
reported only as a frame-consistency sensitivity check, not substituted
for the prespecified 925-pair cohort. Native-ligand chemical identity and
large alignment discrepancies require adjudication before interpreting
the formal numbers as manuscript-ready docking accuracy.

The outstanding provenance request to Lucas is the per-ID prepared input
and status manifest, the exact RCSB receptor files/download and cleaning
procedure, the QM conformer/SDF-generation rule, and run logs for both
models. The full MD HDF5 is **not** needed to score these static poses, but
is needed for any additional MD-derived half of the earlier benchmark.
The current local `MD.hdf5` has a transfer sidecar and fails HDF5 root
parsing despite its apparent full size, so it was not used.

## Updated delivered-pose extension (27 September 2026)

The preceding strict-reference tables are retained as a historical baseline.
The versioned all-copies and multi-residue extension below supersedes their
*coverage* counts; it does not overwrite the old scores or imply that a
manuscript-ready comparison has already been exported. Gnina rescoring is
still in progress locally, and final figures have not been regenerated.

### Delivered denominator and selection limitation

The denominator for every docking coverage percentage is the delivered
primary-pose cohort: 8,419 DiffDock poses and 8,565 EquiBind poses (16,984
method rows). Roughly eleven thousand of the 19,413 QM IDs were never posed
by either method. Lucas attributed missing inputs chiefly to failed RCSB
protein downloads and RDKit ligand-read failures; his input-pipeline counts
are collaborator-reported and do not form a measured model-failure rate.
Neither the official MD train/validation/test split nor trajectory frames
were used for these rigid-receptor docking runs.

### Reference assignment across crystal copies

The extension presents every chemically matching crystal ligand copy to
OpenStructure, allowing its ligand scorer to assign the appropriate copy
as in the co-folding benchmark. This replaces the strict baseline's
prediction-proximity exclusion for multi-copy targets. The assigned copy
is recorded per pose in the versioned score table. Scores from the original
v3 strict run are preserved byte-for-byte when they already exist; an
all-copies result fills only a previously unscored pose.

The reference entity is reduced to heavy atoms before ligand scoring.
This strips deposited hydrogens from native ligands for parity with the
heavy-atom delivered SDFs; no source CIF or SDF is modified. The receptor
is the static crystal polymer, with non-polymer species and waters removed.
These preprocessing choices apply consistently to both docking methods.

### Multi-residue ligand recovery

Short polymer chains (at most 30 residues), branched glycans, and nearby
multi-component native ligands were audited for complete heavy-element and
bond-graph identity with each delivered pose. Graph-valid native assemblies
were collapsed into a single generic pseudo-residue per crystal copy, and
the pose into the same graph, because OpenStructure's ligand scorer requires
a single-residue ligand. Original repeated atom names were replaced by
unique names keyed to verified pose-atom indices; the coordinates and bonds
were retained. All valid native copies were offered for assignment.

The expanded audit examined 3,119 previously excluded method rows over
1,573 IDs. Full graph matches occurred in 1,933 rows: 1,444 short-polymer,
386 branched-glycan, and 103 nearby composite rows. Another 1,159 lacked
an exact elemental candidate, 13 had a bond-count mismatch, and 14 exceeded
the graph-mapping cap. Thus the initial estimate of about 1,590 recoveries
per method was not supported by the chemistry checks. Validation included
native-versus-native controls and short-peptide, glycan, and composite
spot checks; the detailed record is in MULTIRESIDUE_RECOVERY_PROTOCOL.md.

After a higher-limit retry, 15,679 of 16,984 delivered method rows have
both formal OpenStructure BiSyRMSD and lDDT-PLI (92.3%). This comprises
6,633 unchanged v3 rows, 7,118 additional all-copies rows, 1,926
multi-residue rows, and two recovered all-copies timeouts. By method,
formal paired-metric coverage is 7,780/8,419 DiffDock and 7,899/8,565
EquiBind. These are coverage figures among delivered poses, not success
rates over attempted docking inputs or all MISATO complexes.

The remaining 1,305 delivered rows are not silently treated as failures
of a docking model: 1,264 lack an acceptable crystal reference or graph
match, 26 have only partial/no formal metrics, seven have only a separately
labeled fixed-frame direct RMSD, four had scorer exceptions, two failed
reference mapping, and two are the predeclared 5EQQ runtime exclusions.
The multi-residue hard cases 3M3R timed out again and 4YEE exceeded a
runtime guard; they remain among the direct-only diagnostics, not formal
scores. Each non-scored row has an explicit status and reason.

Fixed-frame direct graph RMSD is not substituted into the formal pose-RMSD
column. For multi-residue recovered targets, only 1,631/1,924 initial
formal scores agreed with the direct diagnostic within 0.5 Å; 263 differed
by more than 2 Å. This is expected to be possible because BiSyRMSD
superposes a binding site and can select equivalent crystal copies, whereas
the direct diagnostic remains in the deposited coordinate frame. The seven
direct-only rows retain blank formal RMSD and lDDT-PLI values.

### Confidence and rescoring

The method-native confidence column is DiffDock's rank-1 logit, recovered
for all 8,419 delivered DiffDock poses by exact SHA-256 identity between
the primary rank1.sdf and a confidence-named SDF in the same delivery.
This logit is not a probability. EquiBind did not deliver a native
confidence score, so its method-native confidence stays blank. The
separate rescore_confidence column is reserved for a common gnina
CNNscore on both methods; the two scales must not be conflated.

Gnina v1.3.3 receives the delivered SDF and a temporary receptor PDB
containing the crystal polymer, with waters, non-polymer ligands, and
graph-verified short-polymer native ligands removed. It runs in score-only
mode and reports CNNscore between 0 and 1; the ligand is not redocked or
optimized. The same gnina version, receptor-preparation rule, and score
parser apply to both docking methods. The per-pose status is retained if
gnina cannot score a delivered ligand. The full run is resumable and has
not yet completed at the time of this addendum.

Gnina's CNN models were developed with PDBbind-derived complexes, and
MISATO was itself built from PDBbind 2022. Exact or related complex
overlap may therefore make CNNscore optimistic as an accuracy estimator.
This score is used for common-scale ranking and calibration diagnostics,
not as an independent absolute success claim. Peptide, glycan, and
composite assemblies may also be out of the model's small-molecule
training domain, so their rescoring should be stratified in analysis.
See the gnina software and MISATO data publications linked in
STATIC_METRICS_PROTOCOL.md.

### Static-receptor interface QS

The ligand-inclusive contact score adapts fix_qs_dynamicbind.py to an
unchanged crystal receptor. For each polymer residue, the minimum
heavy-atom distance is measured to the assigned native crystal ligand
and to the delivered pose. Contacts use a 12 Å cutoff. Shared contacts
are weighted by max(0, 1 - absolute distance difference / 12 Å); each
contact present only in one placement adds one nonshared count. QS is
the sum of shared weights divided by that sum plus nonshared count.
No receptor superposition or residue remapping is performed. The
output column is labeled static-receptor QS, never pooled with the
co-folding QS that also reflects receptor prediction error.

Static-receptor QS was obtained for all 15,679 formally scored poses.
Every native-ligand self-control was exactly 1.0; shifting that ligand
well outside the receptor gave 0.0. Across formally scored rows, its
Spearman correlation with lDDT-PLI is 0.653 overall (0.496 DiffDock,
0.864 EquiBind). QS measures interface-contact conservation rather than
atomwise placement: if the two contact sets are identical, this formula
can equal 1.0 even when their shared-contact distances differ.

### Pocket RMSD and receptor representation

Pocket Cα RMSD is not applicable here: both methods placed ligands into
the reference crystal receptor and did not output a redesigned receptor
for comparison. The pocket-RMSD field remains blank, not zero, and
receptor_mode is static_crystal. Comparative figures must filter on this
mode and exclude blank pocket RMSD rather than treating missing values
as perfect receptor prediction. Static-receptor QS and pose RMSD describe
different aspects of ligand placement; neither is renamed pocket RMSD.

### Pending final export and figures

The step-6 evalspreadsheets/misato_v3 export, sidecar provenance CSVs,
and comparative figures must wait for the full 16,984-row gnina status
table. They are not claimed by this addendum. The final plots should
report each method's own formally scored cohort and the paired subset
where both methods scored, with blanks excluded from denominators.
The previously defined exact stereochemistry-match subset remains a
sensitivity analysis. The final Methods/Results numbers should be
regenerated from the completed joined table, not copied from a partial
gnina run or the earlier strict-reference tables.
