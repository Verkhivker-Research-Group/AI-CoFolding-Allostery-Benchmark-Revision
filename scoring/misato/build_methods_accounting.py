"""Rebuild MISATO cohort accounting without treating reported attempts as observed files.

The output is deliberately descriptive, not a claim that Lucas's complete
input/preparation logs were recovered. Counts are grouped by evidence source.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


REPORTED = {
    "qm_hdf5_entries": 19413,
    "qm_ligand_sdfs_written": 19392,
    "qm_ligand_sanitization_failed": 21,
    "rcsb_proteins_saved": 14057,
    "rcsb_proteins_missing": 5356,
    "diffdock_input_pairs": 14037,
    "diffdock_posed": 8419,
    "diffdock_rdkit_skipped": 5610,
    "equibind_posed": 8565,
    "equibind_failed": 10827,
    "equibind_failed_missing_protein": 5354,
}
PUBLISHED = {
    "pdbbind_starting_structures": 19443,
    "md_complexes": 16972,
    "source_url": "https://www.nature.com/articles/s43588-024-00627-2",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def build_accounting(manifest: list[dict[str, str]],
                     identity: list[dict[str, str]],
                     crystal: list[dict[str, str]],
                     split_summary: dict | None,
                     qm_summary: dict | None = None,
                     qm_graph_summary: dict | None = None) -> dict:
    ids = [row["target_id"].upper() for row in manifest]
    if len(ids) != len(set(ids)):
        raise ValueError("Target manifest contains duplicate IDs")
    delivered = {
        method: {row["target_id"].upper() for row in manifest
                 if row[f"{method}_status"] == "candidate"}
        for method in ("diffdock", "equibind")
    }
    raw_dirs = {
        method: sum(row[f"{method}_status"] != "not_delivered" for row in manifest)
        for method in delivered
    }
    paired = delivered["diffdock"] & delivered["equibind"]
    identity_by_id = {row["target_id"].upper(): row for row in identity}
    if len(identity_by_id) != len(identity):
        raise ValueError("Identity audit contains duplicate IDs")
    crystal_by_key = {(row["target_id"].upper(), row["method"]): row for row in crystal}
    if len(crystal_by_key) != len(crystal):
        raise ValueError("Crystal audit contains duplicate target/method pairs")
    crystal_status = {method: Counter() for method in delivered}
    crystal_selection = {method: Counter() for method in delivered}
    strict = {method: set() for method in delivered}
    for method, target_ids in delivered.items():
        for target_id in target_ids:
            row = crystal_by_key.get((target_id, method))
            if row is None:
                crystal_status[method]["not_audited"] += 1
                continue
            crystal_status[method][row["reference_status"]] += 1
            crystal_selection[method][row["reference_selection"]] += 1
            if (row["reference_selection"] == "single_candidate"
                    and row["reference_status"] == "provisional_pose_rmsd"
                    and row["reference_residues"]):
                strict[method].add(target_id)
    strict_paired = {target_id for target_id in strict["diffdock"] & strict["equibind"]
                     if crystal_by_key[(target_id, "diffdock")]["reference_residues"]
                     == crystal_by_key[(target_id, "equibind")]["reference_residues"]}
    identity_status = Counter(identity_by_id[target_id]["disposition"]
                              if target_id in identity_by_id else "not_audited"
                              for target_id in paired)
    reported = dict(REPORTED)
    published = dict(PUBLISHED)
    published["starting_set_minus_lucas_reported_qm_entries"] = (
        PUBLISHED["pdbbind_starting_structures"] - REPORTED["qm_hdf5_entries"])
    reported["diffdock_reported_input_unaccounted"] = (
        REPORTED["diffdock_input_pairs"] - REPORTED["diffdock_posed"]
        - REPORTED["diffdock_rdkit_skipped"])
    reported["equibind_failures_other_than_missing_protein"] = (
        REPORTED["equibind_failed"] - REPORTED["equibind_failed_missing_protein"])
    reported["equibind_reported_outcomes_total"] = (
        REPORTED["equibind_posed"] + REPORTED["equibind_failed"])
    if reported["equibind_reported_outcomes_total"] != REPORTED["qm_ligand_sdfs_written"]:
        raise ValueError("EquiBind reported outcome counts do not reconcile")
    observed = {
        "unique_delivered_target_ids": len(ids),
        "equibind_target_directories": raw_dirs["equibind"],
        "diffdock_target_directories": raw_dirs["diffdock"],
        "equibind_primary_pose_files": len(delivered["equibind"]),
        "diffdock_primary_pose_files": len(delivered["diffdock"]),
        "paired_primary_pose_ids": len(paired),
        "equibind_only_primary_pose_ids": len(delivered["equibind"] - delivered["diffdock"]),
        "diffdock_only_primary_pose_ids": len(delivered["diffdock"] - delivered["equibind"]),
        "paired_identity_dispositions": dict(sorted(identity_status.items())),
        "crystal_status_by_method": {k: dict(sorted(v.items())) for k, v in crystal_status.items()},
        "crystal_selection_by_method": {k: dict(sorted(v.items())) for k, v in crystal_selection.items()},
        "strict_single_reference_diffdock_ids": len(strict["diffdock"]),
        "strict_single_reference_equibind_ids": len(strict["equibind"]),
        "strict_single_reference_paired_ids": len(strict_paired),
        "strict_paired_by_prediction_identity": dict(sorted(Counter(
            identity_by_id[target_id]["disposition"] if target_id in identity_by_id
            else "not_audited" for target_id in strict_paired).items())),
    }
    if observed["equibind_primary_pose_files"] != REPORTED["equibind_posed"]:
        raise ValueError("Delivered EquiBind count differs from Lucas's report")
    if observed["diffdock_primary_pose_files"] != REPORTED["diffdock_posed"]:
        raise ValueError("Delivered DiffDock count differs from Lucas's report")
    if split_summary:
        observed["published_md_split_ids"] = split_summary["published_id_count"]
        if observed["published_md_split_ids"] != PUBLISHED["md_complexes"]:
            raise ValueError("Local split-list count differs from the MISATO paper")
        observed["delivered_ids_outside_md_splits"] = split_summary["outside_published_md_splits"]["delivered_ids"]
        observed["delivered_ids_within_md_splits"] = len(ids) - observed["delivered_ids_outside_md_splits"]
    if qm_summary:
        if qm_summary["delivered_unique_ids"] != len(ids):
            raise ValueError("QM coverage audit and target manifest disagree")
        observed["verified_qm_hdf5_group_ids"] = qm_summary["qm_group_ids"]
        for key in ("delivered_ids_in_qm", "delivered_ids_outside_qm",
                    "qm_ids_without_delivered_directory", "diffdock_primary_in_qm",
                    "equibind_primary_in_qm", "paired_primary_in_qm",
                    "diffdock_primary_outside_qm", "equibind_primary_outside_qm",
                    "paired_primary_outside_qm"):
            observed[key] = qm_summary[key]
        observed["pose_heavy_element_match_qm_by_method"] = qm_summary.get(
            "pose_heavy_element_match_qm_by_method", {})
        observed["qm_hdf5_matches_lucas_reported_count"] = (
            qm_summary["qm_group_ids"] == REPORTED["qm_hdf5_entries"])
    if qm_graph_summary:
        observed["pose_qm_heavy_graph_by_method"] = qm_graph_summary["by_method"]
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "evidence": {
            "published_primary_source": published,
            "collaborator_reported_not_independently_verified": reported,
            "directly_observed_delivered_files_and_reference_audit": observed,
        },
        "interpretation": [
            "Successful poses were not delivered for the complete reported 19,413-entry QM set.",
            ("The public QM entry count, ID membership, and pose heavy-element compositions were checksum-verified; "
             "prepared ligands, protein downloads, attempted inputs, and failure causes "
             "still require Lucas's input/status logs." if qm_summary else
             "The reported QM entry count, prepared ligands, protein downloads, attempted inputs, and failure causes cannot be independently verified without QM.hdf5 and Lucas's input/status logs."),
            "The eight DiffDock input pairs absent from the reported posed/skipped breakdown remain unresolved.",
            "Official MD train/val/test lists are overlap context, not a docking-input denominator: Lucas says he did not use these splits or MD frames.",
            "Strict single-reference counts are eligibility estimates, not completed OpenStructure scores; they exclude multi-copy, composite, unmatched, and graph-failed cases.",
        ],
    }


def methods_text(result: dict) -> str:
    r = result["evidence"]["collaborator_reported_not_independently_verified"]
    p = result["evidence"]["published_primary_source"]
    o = result["evidence"]["directly_observed_delivered_files_and_reference_audit"]
    qm_verified = o.get("verified_qm_hdf5_group_ids")
    qm_sentence = (f"The checksum-verified public `QM.hdf5` contains {qm_verified:,} "
                   f"root groups and {o['delivered_ids_in_qm']:,} of the delivered IDs; "
                   f"{o['delivered_ids_outside_qm']:,} delivered IDs are absent from QM."
                   if qm_verified is not None else
                   "The public `QM.hdf5` has not been locally verified.")
    chemistry = o.get("pose_heavy_element_match_qm_by_method", {})
    qm_chem_sentence = (f"The delivered pose/QM heavy-element matches are "
                        f"{chemistry['diffdock'].get('match', 0):,} for DiffDock and "
                        f"{chemistry['equibind'].get('match', 0):,} for EquiBind. "
                        "This verifies elemental composition, not atom mapping, "
                        "stereochemistry, or Lucas's exact starting conformer."
                        if "diffdock" in chemistry and "equibind" in chemistry else
                        "Pose/QM elemental composition has not been audited.")
    graphs = o.get("pose_qm_heavy_graph_by_method", {})
    qm_graph_sentence = (f"All {graphs['diffdock'].get('same_indexed_graph', 0):,} "
                         f"DiffDock and {graphs['equibind'].get('same_indexed_graph', 0):,} "
                         "EquiBind primary poses retain the QM heavy-atom ordering "
                         "and untyped connectivity. This does not validate bond order, "
                         "stereochemistry, protonation, or starting coordinates."
                         if "diffdock" in graphs and "equibind" in graphs else
                         "Pose/QM heavy-atom connectivity has not been audited.")
    qm_denominator = qm_verified if qm_verified is not None else r["qm_hdf5_entries"]
    return f"""# MISATO docking Methods accounting (draft; generated)

This file separates collaborator-reported preparation/attempt counts from quantities
independently observed in the delivered output folders. It is **not** evidence
that every QM entry was submitted, scored, or delivered.

## Data and input provenance

Lucas reported {r['qm_hdf5_entries']:,} ligand entries in his QM HDF5. The
MISATO paper states that its *starting PDBbind set* had
{p['pdbbind_starting_structures']:,} protein–ligand structures, which is
{p['starting_set_minus_lucas_reported_qm_entries']:,} more; these are different
processing stages, not necessarily an inconsistency.
He reported {r['qm_ligand_sdfs_written']:,} SDFs written and
{r['qm_ligand_sanitization_failed']:,} sanitization failures, and
{r['rcsb_proteins_saved']:,} RCSB protein files saved versus
{r['rcsb_proteins_missing']:,} missing. His stated docking inputs used static
RCSB crystal proteins and one QM-derived ligand conformer, **not** MD frames
or the published MD train/validation/test split. The preparation files and
per-target status logs have not been supplied, so these are reported rather
than independently verified preparation counts.

{qm_sentence}
{qm_chem_sentence}
{qm_graph_sentence}

## Delivered poses and completion

The delivered trees contain {o['diffdock_primary_pose_files']:,} DiffDock
`rank1.sdf` poses in {o['diffdock_target_directories']:,} target directories
and {o['equibind_primary_pose_files']:,} EquiBind
`lig_equibind_corrected.sdf` poses in {o['equibind_target_directories']:,}
target directories. There are {o['unique_delivered_target_ids']:,} distinct
target IDs and {o['paired_primary_pose_ids']:,} IDs with both primary poses.
The DiffDock folder contains other rank/confidence SDFs, which are retained
as provenance but are not pooled into the primary rank-1 result.

These observed primary pose counts match Lucas's reported successful-pose
counts, but represent only {o['diffdock_primary_pose_files'] / qm_denominator:.1%}
(DiffDock) and {o['equibind_primary_pose_files'] / qm_denominator:.1%}
(EquiBind) of the {'verified' if qm_verified is not None else 'reported'} QM-entry count.
Those percentages compare output counts to the QM set, not attempted inputs or verified input success
rates. A complete successful-pose computation is therefore ruled out, while
the completeness of attempted runs remains unknown.

Lucas reported {r['diffdock_input_pairs']:,} DiffDock input pairs,
{r['diffdock_posed']:,} poses, and {r['diffdock_rdkit_skipped']:,} RDKit-unreadable
ligand skips. The latter two sum to {r['diffdock_posed'] + r['diffdock_rdkit_skipped']:,},
leaving {r['diffdock_reported_input_unaccounted']:,} input pairs unexplained.
He reported {r['equibind_posed']:,} EquiBind poses and {r['equibind_failed']:,}
failures, which sum to the {r['qm_ligand_sdfs_written']:,} reported SDFs.
Of those failures, {r['equibind_failed_missing_protein']:,} were attributed to
missing proteins and {r['equibind_failures_other_than_missing_protein']:,}
have no finer verified status.

## Reference and scoring cohort

The RCSB crystal-reference audit admits {o['strict_single_reference_diffdock_ids']:,}
DiffDock and {o['strict_single_reference_equibind_ids']:,} EquiBind delivered
poses with exactly one heavy-element-matched native ligand, a usable direct
coordinate-frame graph check, and a cached reference. Their common strict
cohort has {o['strict_single_reference_paired_ids']:,} IDs. These are
**pre-scoring eligibility** counts, not final OST-scored denominators.
Of those pairs, {o['strict_paired_by_prediction_identity'].get('eligible_exact_graph', 0):,}
have an exact graph/stereochemistry match between the two predictions;
{o['strict_paired_by_prediction_identity'].get('eligible_stereochemistry_review', 0):,}
have differing stereochemical representations and
{o['strict_paired_by_prediction_identity'].get('eligible_representation_review', 0):,}
require other representation review. The exact-between-methods subset is
the conservative primary paired comparison; the larger topology-matched
set is a sensitivity analysis. Neither test alone proves native crystal
stereochemistry.
Multi-copy references selected using a docking pose, composite native
ligands, unmatched ligands, and graph failures are excluded from the primary
paired comparison until independently resolved. The separate score output
must report metric-specific successes and failures; missing values must not
be silently dropped from the denominator.

The static receptor makes pocket Cα RMSD and whole-complex QS comparisons
non-independent of the provided crystal structure. Pose RMSD and lDDT-PLI
are the primary docking comparisons. The original pipeline runs locally
under WSL Ubuntu (`plb` conda environment), with parallel per-target scoring;
Docker is not required. Reference chemistry, coordinate frame, and the exact
Lucas receptor files remain limitations pending his preparation files.

## MD-split context

The official MD lists contain {o.get('published_md_split_ids', 0):,} IDs;
{o.get('delivered_ids_within_md_splits', 0):,} delivered IDs overlap them and
{o.get('delivered_ids_outside_md_splits', 0):,} do not. This is *not* an
estimate of docking coverage because the docking inputs were QM-derived and
Lucas says he did not apply the MD splits.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("misato_output/inventory/target_manifest.csv"))
    parser.add_argument("--identity", type=Path, default=Path("misato_output/inventory/primary_ligand_identity.csv"))
    parser.add_argument("--crystal-audit", type=Path, default=Path("misato_output/crystal_pose_all.csv"))
    parser.add_argument("--split-summary", type=Path, default=Path("misato_output/split_comparison/summary.json"))
    parser.add_argument("--qm-summary", type=Path, default=Path("misato_output/qm_coverage_v3/qm_coverage_summary.json"))
    parser.add_argument("--qm-graph-summary", type=Path,
                        default=Path("misato_output/pose_qm_graph_audit_v1/pose_qm_graph_summary.json"))
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.exists():
        parser.error("--out-dir already exists; choose a new versioned directory")
    split = json.loads(args.split_summary.read_text(encoding="utf-8")) if args.split_summary.exists() else None
    qm = json.loads(args.qm_summary.read_text(encoding="utf-8")) if args.qm_summary.exists() else None
    qm_graph = json.loads(args.qm_graph_summary.read_text(encoding="utf-8")) if args.qm_graph_summary.exists() else None
    result = build_accounting(read_csv(args.manifest), read_csv(args.identity),
                              read_csv(args.crystal_audit), split, qm, qm_graph)
    args.out_dir.mkdir(parents=True)
    (args.out_dir / "cohort_accounting.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    (args.out_dir / "METHODS_ACCOUNTING.md").write_text(methods_text(result), encoding="utf-8")
    print(json.dumps(result["evidence"]["directly_observed_delivered_files_and_reference_audit"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
