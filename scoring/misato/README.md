# MISATO result intake

For the delivered-pose extension, see
[METHODS_RESULTS.md](METHODS_RESULTS.md) for the versioned Methods/Results
addendum, [MULTIRESIDUE_RECOVERY_PROTOCOL.md](MULTIRESIDUE_RECOVERY_PROTOCOL.md)
for steps 1–2, and [STATIC_METRICS_PROTOCOL.md](STATIC_METRICS_PROTOCOL.md)
for steps 3–5. The latter's gnina run is resumable and its final join must
wait for every delivered pose; the earlier strict-reference CSVs remain
untouched. [RUN_STEPS_3_TO_5.md](RUN_STEPS_3_TO_5.md) contains the visible
local commands and progress check.

`create_manifest.py` inventories the supplied EquiBind and DiffDock trees without changing them. It writes a SHA-256 raw-file manifest, a per-target cohort manifest, and an inventory summary.

Only `lig_equibind_corrected.sdf` (EquiBind) and `rank1.sdf` (DiffDock) are primary candidates. Other DiffDock `rank*_confidence-*.sdf` files are retained as repeat-run evidence, never silently pooled into the primary result.

```powershell
python scoring/misato/create_manifest.py `
  --equibind-root 'C:\Users\Ryan\Downloads\misato_equibind_results\nfshome\turano\Datasets\misato\equibind_output' `
  --diffdock-root 'C:\Users\Ryan\Downloads\MisatoDiffdockResults\MisatoDiffdockResults' `
  --out-dir misato_output\inventory
```

The delivered files are ligand-only. This tool deliberately does not calculate RMSD: protein-aligned docking metrics require validation of the receptor/reference coordinate frame and ligand identity.

## New crystal-reference path (separate from the old benchmark)

`fetch_crystal_cifs.py` uses the delivered target manifest, strips any suffix
only for the four-character RCSB lookup, and saves **asymmetric-unit** entry
CIFs in `misato_output/rcsb_asymmetric_unit/`. It never changes the existing
ASD/PLA reference cache or overwrites an existing CIF or report. A bounded
test run is:

```powershell
python scoring/misato/fetch_crystal_cifs.py `
  --limit 20 `
  --report misato_output/rcsb_asymmetric_unit/fetch_first20.csv
```

Omit `--limit` for all delivered targets and choose a new `--report` path.
Existing valid CIFs are reused. The script requires `requests`.

`audit_crystal_poses.py` reads the new cache plus the original ligand-only
SDFs, finds native crystal ligands by heavy-element signature, and attempts
element/connectivity-aware direct-coordinate RMSD for single-component
matches. It uses **one selected crystal copy for both methods**. Composite,
ambiguous, unmatched, and unreadable cases are flagged instead of forced into
a score. A bounded run is:

```powershell
python scoring/misato/audit_crystal_poses.py `
  --limit 20 `
  --out-csv misato_output/crystal_pose_first20.csv
```

The output column `pose_rmsd_angstrom_provisional` is exploratory: chemical
bonds in the native ligand are inferred from coordinates, and when multiple
crystal copies exist the nearest copy is selected using DiffDock's position
(or EquiBind's if DiffDock is unavailable). This can favor the anchor method.
The result is **not** the original pipeline's OpenStructure BiSyRMSD or
lDDT-PLI. Review graph/coordinate-frame status and reference selection before
aggregating any values. The script needs `gemmi`, `numpy`, and `rdkit`; both
scripts refuse to overwrite an existing output report.

## Compare target IDs with the published MISATO MD splits

`compare_published_splits.py` compares the target manifest against Zenodo's
`train_MD.txt`, `val_MD.txt`, and `test_MD.txt`. It downloads only these three
small text files, verifies their published MD5 checksums, and saves them with
the outputs so the audit can be repeated offline. It never reads `MD.hdf5`.
Lucas reports that these splits were **not** used to select docking inputs:
his ligands came from QM coordinates and the proteins were static RCSB
structures. This script therefore reports MD-split overlap only, not the
actual filtering rate or denominator for his runs.

```powershell
python scoring/misato/compare_published_splits.py `
  --target-manifest misato_output/inventory/target_manifest.csv `
  --out-dir misato_output/split_comparison
```

To rerun without network access, add
`--split-dir misato_output/split_comparison/published_lists`.
The script writes `target_comparison.csv` (one row per ID in either source),
`split_summary.csv`, and `summary.json`. `paired_primary_candidate` means both
expected ligand files are present; it does **not** certify molecular identity,
MD availability, or docking accuracy. Targets outside the published split lists
are retained and labeled, not silently discarded. An absent target tells us
which published MD IDs were not delivered, but not why Lucas excluded them.

`audit_references.py` performs that next, staged check against RCSB experimental
structures. Start with a bounded sample; it caches references locally and
identifies non-polymer residue candidates whose heavy-element signature matches
the paired predictions.

```powershell
python scoring/misato/audit_references.py `
  --identity-csv misato_output\inventory\primary_ligand_identity.csv `
  --cache-dir misato_output\references `
  --out-csv misato_output\inventory\reference_audit_sample.csv `
  --limit 100
```

`audit_hdf5_entries.py` is an optional audit of the original MISATO MD HDF5,
not part of reproducing Lucas's static-receptor docking inputs. It reads only
metadata for requested target groups (not all trajectory coordinates) and
reports the ligand-selection rule used by the public MISATO preprocessing
code. Run it in an environment with `h5py` and `numpy`:

```powershell
python scoring/misato/audit_hdf5_entries.py `
  --h5 'D:\MISATO\MD.hdf5' `
  --target-manifest misato_output\inventory\target_manifest.csv `
  --out-csv misato_output\inventory\misato_hdf5_audit.csv
```
