# MISATO multi-residue reference recovery: validation record

This is the step 1–2 validation for **delivered poses only**. It does not
regenerate docking inputs or rerun DiffDock/EquiBind. Original v3 scores and
the all-copies run remain unchanged. The recovery results are versioned under
the Git-ignored `misato_output/` directory until reconciled.

## Eligibility and chemistry

The existing crystal audit examined individual HETATM residues. The new
`audit_multiresidue_recovery.py` additionally groups RCSB mmCIF polymer and
branched-oligosaccharide subchains of 2–30 residues and tests nearby pairs of
non-polymer components. A component pair must have atoms within 2.3 Å;
element counts alone must not combine distant copies. A candidate is
eligible only if the complete heavy-element composition **and** full untyped
bond graph match the delivered pose. Equal atom counts or a one-way
substructure match are insufficient: a disconnected two-component candidate
can otherwise match a connected pose by ignoring the pose's extra bond.

The full local audit considered 3,119 previously excluded method rows over
1,573 target IDs. It found 1,933 full-graph matches (1,444 short-polymer,
386 branched-glycan, 103 nearby non-polymer-pair method rows), 1,159 without
an exact elemental candidate, 13 with a bond-count mismatch, and 14
exceeding the graph-mapping cap. These are **eligibility** counts, not
completed OST scores; they do not substantiate the preliminary
~1,590-per-method recovery estimate. The per-pose audit is
`misato_output/multiresidue_recovery_all_v2.csv`.

For graph-valid assemblies, `score_multiresidue_recovery_wsl.py` constructs
one generic OST pseudo-residue per native crystal copy, plus one for the
delivered pose. Since different residues reuse names such as `CA` and `C`,
the pseudo-residue uses unique names keyed to the verified pose atom index;
it preserves a mapping to the original atom coordinates and bonds rather
than reusing ambiguous names. All graph-valid crystal copies are passed to
OST for assignment. The static receptor contains the original crystal
protein/nucleic-acid polymer coordinates, with the candidate ligand
subchains and non-polymer species removed. No source CIF or SDF is changed.
The reference and model pseudo-residues use the same complete heavy-atom
graph; substructure-only matching is disabled.

Method-validation probes on short peptides (`1A07`, `1A30`, `1APV`), a
branched glycan (`1AGM`), and a two-component ligand (`11GS`) returned both
BiSyRMSD and lDDT-PLI. Native-against-itself controls returned lDDT-PLI
1.0 and RMSD within numerical roundoff of zero. The `11GS` check exposed false matches
between crystal components ~18 Å apart; the 2.3 Å component-pair rule and
full bond-graph check remove those combinations. These probes do not replace
the complete per-row failure audit or prove that every candidate is chemically
identical in bond order/stereochemistry to the deposited ligand.

The fixed-frame direct graph RMSD is retained only as a separately labeled
fallback/sensitivity metric if OST cannot score a verified multi-residue
candidate. It must not be pooled with OST's binding-site-superposed BiSyRMSD,
and lDDT-PLI stays missing for a direct-only fallback. In the completed
multi-residue score, only 1,631/1,924 (84.8%) initial formal scores were
within 0.5 Å of their direct-RMSD diagnostic; 263 differed by >2 Å and 174
by >10 Å. The previously observed 97% single-chain agreement does **not**
validate substitution for these assemblies.

The first full multi-residue pass formally scored 1,924/1,933 graph-valid
method rows. A 300-second per-stage retry recovered two more (`3LRH`
EquiBind and `3ODI` EquiBind); their DiffDock partners had already scored.
Both `3M3R` methods still timed out. Three EquiBind rows remain partial,
and both `4YEE` methods are explicitly excluded by a runtime guard after
an OST operation exceeded eight minutes. The retry and first-pass CSVs
record those statuses rather than assigning zero-valued metrics.

The versioned merge `misato_output/misato_score_recovered_v1.csv` contains
all 16,984 delivered pose rows. It has 15,679 formal pose-RMSD/lDDT-PLI
pairs (92.3% of delivered), seven separately labeled direct-only fallbacks,
26 partial rows, 1,264 reference/graph-check exclusions, four scoring
exceptions, two mapping failures, and two runtime guards. Each non-scored
row has an explicit reason. The 6,633 original v3 scored rows retain their
metric strings exactly; 7,118 more formal rows come from the original
all-copies pass, 1,926 from multi-residue scoring/retry, and two from the
old timeout retry. This is a step 1–2 intermediate table, **not** the final
v3 spreadsheet export with QS and confidence columns.

## Residual all-copies failures

The 90-second all-copies run had 34 non-`5EQQ` residual method rows:
26 partial/no-metric rows, four invalid-handle exceptions, two
`no_ost_native_copy` mapping failures, and two runtime timeouts. A 300-second
retry scored both timeouts (`2VXA` and `4QIJ` EquiBind) without redocking;
`5EQQ` remains predeclared as excluded. OST's diagnostic state API explains
the 26 partial rows as 10 disconnected ligands, 10 absent model binding
sites (where lDDT-PLI can still be available), and six identity mismatches.
The per-row diagnostic is
`misato_output/allcopies_residual_diagnosis_v1.csv`. Unresolved rows are
retained as failures with their actual reason, not imputed as zero.

## Reproduction (WSL Ubuntu, `plb` environment)

From PowerShell in this checkout:

```powershell
python scoring/misato/audit_multiresidue_recovery.py --out-csv misato_output/multiresidue_recovery_all_v3.csv
wsl -d Ubuntu -- bash -lic 'conda activate plb && cd /mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark-Revision-W-Misato && python scoring/misato/score_multiresidue_recovery_wsl.py --audit misato_output/multiresidue_recovery_all_v3.csv --workers 6 --timeout-seconds 120 --skip-target-id 4YEE --out-csv misato_output/multiresidue_score_all_v2.csv'
```

Existing outputs are never overwritten. Use `--resume` only with a partially
completed scoring CSV made with the same code and audit version. Merge into a
new, status-preserving score table with `merge_recovery_scores.py` **after**
the complete recovery run; its assertion requires all original v3 scored
values to remain byte-for-byte unchanged.

The actual run used the `v2` audit, `multiresidue_score_all_v1.csv`, and
`multiresidue_retry_v1.csv`. Reproduce that merge without overwriting it:

```powershell
python scoring/misato/merge_recovery_scores.py --recovery misato_output/multiresidue_score_all_v1.csv --recovery-retry misato_output/multiresidue_retry_v1.csv --out-csv misato_output/misato_score_recovered_v2.csv
```
