"""Produce metric-specific MISATO denominators and paired rank-1 summaries."""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import Counter
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def number(value: str) -> float | None:
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def describe(values: list[float]) -> dict[str, float | int | None]:
    return {"n": len(values),
            "mean": statistics.mean(values) if values else None,
            "median": statistics.median(values) if values else None}


def rmsd_disagreements(scores: list[dict[str, str]],
                       audit: list[dict[str, str]],
                       threshold: float = 0.5) -> list[dict[str, object]]:
    """Keep every aligned-versus-fixed-frame discrepancy, not just examples."""
    audit_by_key = {(row["target_id"], row["method"]): row for row in audit}
    rows = []
    for score in scores:
        key = (score["target_id"], score["method"])
        old = number(audit_by_key[key].get("pose_rmsd_angstrom_provisional", ""))
        formal = number(score.get("bisy_rmsd_angstrom", ""))
        if old is not None and formal is not None and abs(formal - old) > threshold:
            rows.append({"target_id": key[0], "method": key[1],
                         "formal_angstrom": formal,
                         "provisional_angstrom": old,
                         "absolute_difference_angstrom": abs(formal - old)})
    return sorted(rows, key=lambda row: -float(row["absolute_difference_angstrom"]))


def paired_metrics(rows: list[dict[str, object]]) -> dict:
    rmsd_pairs = [(float(row["diffdock_bisy_rmsd_angstrom"]),
                   float(row["equibind_bisy_rmsd_angstrom"]))
                  for row in rows if row["diffdock_bisy_rmsd_angstrom"] != ""
                  and row["equibind_bisy_rmsd_angstrom"] != ""]
    lddt_pairs = [(float(row["diffdock_lddt_pli"]), float(row["equibind_lddt_pli"]))
                  for row in rows if row["diffdock_lddt_pli"] != ""
                  and row["equibind_lddt_pli"] != ""]
    return {
        "eligible_paired_ids_in_score_file": len(rows),
        "both_bisy_rmsd_available": len(rmsd_pairs),
        "both_lddt_pli_available": len(lddt_pairs),
        "bisy_rmsd_diffdock": describe([pair[0] for pair in rmsd_pairs]),
        "bisy_rmsd_equibind": describe([pair[1] for pair in rmsd_pairs]),
        "lddt_pli_diffdock": describe([pair[0] for pair in lddt_pairs]),
        "lddt_pli_equibind": describe([pair[1] for pair in lddt_pairs]),
        "diffdock_lower_rmsd": sum(a < b for a, b in rmsd_pairs),
        "equibind_lower_rmsd": sum(b < a for a, b in rmsd_pairs),
        "rmsd_ties": sum(a == b for a, b in rmsd_pairs),
        "diffdock_higher_lddt": sum(a > b for a, b in lddt_pairs),
        "equibind_higher_lddt": sum(b > a for a, b in lddt_pairs),
        "lddt_ties": sum(a == b for a, b in lddt_pairs),
        "rmsd_le_2_angstrom_diffdock": sum(a <= 2 for a, _ in rmsd_pairs),
        "rmsd_le_2_angstrom_equibind": sum(b <= 2 for _, b in rmsd_pairs),
    }


def summarize(scores: list[dict[str, str]], audit: list[dict[str, str]],
              allow_partial: bool = False,
              identity: list[dict[str, str]] | None = None) -> tuple[dict, list[dict[str, object]]]:
    audit_by_key = {(row["target_id"], row["method"]): row for row in audit}
    score_by_key = {(row["target_id"], row["method"]): row for row in scores}
    if len(audit_by_key) != len(audit) or len(score_by_key) != len(scores):
        raise ValueError("Duplicate target/method rows")
    extra = set(score_by_key) - set(audit_by_key)
    missing = set(audit_by_key) - set(score_by_key)
    if extra or (missing and not allow_partial):
        raise ValueError(f"Score/audit key mismatch: {len(extra)} extra, {len(missing)} missing")
    identity_by_id = {row["target_id"]: row["disposition"] for row in (identity or [])}
    if identity is not None and len(identity_by_id) != len(identity):
        raise ValueError("Duplicate identity-audit target IDs")
    scored = {method: [row for row in scores if row["method"] == method]
              for method in ("diffdock", "equibind")}
    summary: dict = {
        "expected_delivered_pose_rows": len(audit),
        "observed_score_rows": len(scores),
        "missing_score_rows": len(missing),
        "complete_score_file": not missing,
        "method": {},
    }
    rmsd_disagreements_rows = rmsd_disagreements(scores, audit)
    summary["quality_control"] = {
        "formal_vs_provisional_rmsd_disagreement_gt_0_5_angstrom": len(rmsd_disagreements_rows),
        "disagreement_gt_2_angstrom": sum(row["absolute_difference_angstrom"] > 2
                                          for row in rmsd_disagreements_rows),
        "disagreement_gt_10_angstrom": sum(row["absolute_difference_angstrom"] > 10
                                           for row in rmsd_disagreements_rows),
        "disagreement_by_method": dict(sorted(Counter(row["method"]
                                                  for row in rmsd_disagreements_rows).items())),
        "disagreement_by_prediction_identity": dict(sorted(Counter(
            identity_by_id.get(row["target_id"], "not_audited")
            for row in rmsd_disagreements_rows).items())),
        "examples_first_30": rmsd_disagreements_rows[:30],
        "interpretation": "Flag for manual review, not automatic exclusion; the two algorithms are not identical.",
    }
    for method, rows in scored.items():
        eligible = [row for row in rows if row["reference_selection"] == "single_candidate"
                    and row["preaudit_status"] == "provisional_pose_rmsd"]
        rmsd = [value for row in eligible
                if (value := number(row["bisy_rmsd_angstrom"])) is not None]
        lddt = [value for row in eligible
                if (value := number(row["lddt_pli"])) is not None]
        summary["method"][method] = {
            "delivered_rows": len(rows), "strict_reference_eligible_rows": len(eligible),
            "status_counts": dict(sorted(Counter(row["status"] for row in rows).items())),
            "bisy_rmsd_angstrom": describe(rmsd), "lddt_pli": describe(lddt),
            "rmsd_le_2_angstrom": sum(value <= 2 for value in rmsd),
        }
    common = sorted({target for target, method in score_by_key if method == "diffdock"}
                    & {target for target, method in score_by_key if method == "equibind"})
    paired_rows: list[dict[str, object]] = []
    for target_id in common:
        diff = score_by_key[(target_id, "diffdock")]
        equi = score_by_key[(target_id, "equibind")]
        if not all(row["reference_selection"] == "single_candidate"
                   and row["preaudit_status"] == "provisional_pose_rmsd"
                   for row in (diff, equi)):
            continue
        if diff["reference_residue"] != equi["reference_residue"]:
            continue
        dr = number(diff["bisy_rmsd_angstrom"])
        er = number(equi["bisy_rmsd_angstrom"])
        dl = number(diff["lddt_pli"])
        el = number(equi["lddt_pli"])
        ddirect = number(audit_by_key[(target_id, "diffdock")].get("pose_rmsd_angstrom_provisional", ""))
        edirect = number(audit_by_key[(target_id, "equibind")].get("pose_rmsd_angstrom_provisional", ""))
        paired_rows.append({"target_id": target_id, "reference_residue": diff["reference_residue"],
                            "prediction_identity_disposition": identity_by_id.get(target_id, "not_audited"),
                            "diffdock_bisy_rmsd_angstrom": "" if dr is None else dr,
                            "equibind_bisy_rmsd_angstrom": "" if er is None else er,
                            "diffdock_direct_rmsd_angstrom_provisional": "" if ddirect is None else ddirect,
                            "equibind_direct_rmsd_angstrom_provisional": "" if edirect is None else edirect,
                            "diffdock_formal_minus_direct_angstrom": "" if dr is None or ddirect is None else dr - ddirect,
                            "equibind_formal_minus_direct_angstrom": "" if er is None or edirect is None else er - edirect,
                            "diffdock_lddt_pli": "" if dl is None else dl,
                            "equibind_lddt_pli": "" if el is None else el,
                            "diffdock_status": diff["status"], "equibind_status": equi["status"]})
    summary["paired_primary_exact_prediction_graph"] = paired_metrics(
        [row for row in paired_rows
         if row["prediction_identity_disposition"] == "eligible_exact_graph"])
    summary["paired_broad_topology_matched_sensitivity"] = paired_metrics(paired_rows)
    frame_consistent = []
    for row in paired_rows:
        if row["prediction_identity_disposition"] != "eligible_exact_graph":
            continue
        target_id = str(row["target_id"])
        checks = []
        for method in ("diffdock", "equibind"):
            formal = number(score_by_key[(target_id, method)]["bisy_rmsd_angstrom"])
            direct = number(audit_by_key[(target_id, method)]["pose_rmsd_angstrom_provisional"])
            checks.append(formal is not None and direct is not None
                          and abs(formal - direct) <= 0.5)
        if all(checks):
            frame_consistent.append(row)
    summary["paired_exact_graph_frame_consistent_sensitivity"] = paired_metrics(frame_consistent)
    direct_pairs = []
    for target_id in common:
        da = audit_by_key[(target_id, "diffdock")]
        ea = audit_by_key[(target_id, "equibind")]
        if not all(row.get("reference_selection") == "single_candidate"
                   and row.get("reference_status") == "provisional_pose_rmsd"
                   for row in (da, ea)):
            continue
        if da["reference_residues"] != ea["reference_residues"]:
            continue
        dr = number(da.get("pose_rmsd_angstrom_provisional", ""))
        er = number(ea.get("pose_rmsd_angstrom_provisional", ""))
        if dr is not None and er is not None:
            direct_pairs.append((dr, er))
    summary["sensitivity_direct_coordinate_rmsd_not_formal_bisyrmsd"] = {
        "paired_ids": len(direct_pairs),
        "diffdock_angstrom": describe([pair[0] for pair in direct_pairs]),
        "equibind_angstrom": describe([pair[1] for pair in direct_pairs]),
        "diffdock_le_2_angstrom": sum(a <= 2 for a, _ in direct_pairs),
        "equibind_le_2_angstrom": sum(b <= 2 for _, b in direct_pairs),
        "diffdock_lower": sum(a < b for a, b in direct_pairs),
        "equibind_lower": sum(b < a for a, b in direct_pairs),
        "interpretation": "Exploratory sensitivity on the full strict-reference cohort; not the original OST BiSyRMSD metric.",
    }
    exact_direct = [(a, b) for target_id in common
                    if identity_by_id.get(target_id) == "eligible_exact_graph"
                    if (da := audit_by_key[(target_id, "diffdock")])["reference_selection"] == "single_candidate"
                    if (ea := audit_by_key[(target_id, "equibind")])["reference_selection"] == "single_candidate"
                    if da["reference_status"] == ea["reference_status"] == "provisional_pose_rmsd"
                    if da["reference_residues"] == ea["reference_residues"]
                    if (a := number(da["pose_rmsd_angstrom_provisional"])) is not None
                    if (b := number(ea["pose_rmsd_angstrom_provisional"])) is not None]
    summary["sensitivity_direct_coordinate_rmsd_exact_prediction_graph"] = {
        "paired_ids": len(exact_direct),
        "diffdock_angstrom": describe([pair[0] for pair in exact_direct]),
        "equibind_angstrom": describe([pair[1] for pair in exact_direct]),
        "diffdock_le_2_angstrom": sum(a <= 2 for a, _ in exact_direct),
        "equibind_le_2_angstrom": sum(b <= 2 for _, b in exact_direct),
        "diffdock_lower": sum(a < b for a, b in exact_direct),
        "equibind_lower": sum(b < a for a, b in exact_direct),
    }
    return summary, paired_rows


def methods_text(summary: dict) -> str:
    p = summary["paired_primary_exact_prediction_graph"]
    broad = summary["paired_broad_topology_matched_sensitivity"]
    frame = summary["paired_exact_graph_frame_consistent_sensitivity"]
    direct = summary["sensitivity_direct_coordinate_rmsd_exact_prediction_graph"]
    m = summary["method"]
    state = "complete" if summary["complete_score_file"] else "INCOMPLETE — smoke/partial run"
    return f"""# MISATO static-docking score summary ({state})

The source is the delivered primary ligand SDF for each method (`rank1.sdf`
for DiffDock; `lig_equibind_corrected.sdf` for EquiBind), the RCSB
asymmetric-unit crystal CIF, and the unchanged static crystal receptor. The
scorer uses the existing WSL `plb` OpenStructure BiSyRMSD and lDDT-PLI
implementation. It does not use MD frames or select an oracle-best pose.
The reference passed to OST is restricted to the one pre-audited ligand,
and substructure-only matches are disabled. The model complex uses a
temporary heavy-atom-only copy of the SDF so explicit ligand hydrogens do
not force a partial match against the crystal reference. OST BiSyRMSD includes a local
binding-site superposition and may map equivalent protein chains; it is
not identical to fixed-frame direct-coordinate RMSD.

Reference eligibility requires exactly one crystal nonpolymer ligand with
the same heavy-element signature and a successful untyped graph/coordinate
audit. This conservative rule excludes multi-copy references selected by
DiffDock proximity, composite crystal ligands, unmatched ligands, and graph
failures. Metric failures retain their rows and are not imputed as zero.

Expected delivered pose rows: {summary['expected_delivered_pose_rows']:,};
score-file rows: {summary['observed_score_rows']:,}; missing: {summary['missing_score_rows']:,}.
DiffDock delivered/strict-eligible/RMSD-scored/lDDT-scored:
{m['diffdock']['delivered_rows']:,} / {m['diffdock']['strict_reference_eligible_rows']:,} /
{m['diffdock']['bisy_rmsd_angstrom']['n']:,} / {m['diffdock']['lddt_pli']['n']:,}.
EquiBind delivered/strict-eligible/RMSD-scored/lDDT-scored:
{m['equibind']['delivered_rows']:,} / {m['equibind']['strict_reference_eligible_rows']:,} /
{m['equibind']['bisy_rmsd_angstrom']['n']:,} / {m['equibind']['lddt_pli']['n']:,}.

The conservative primary comparison requires the same exact ligand graph in
the two predicted SDFs in addition to the strict crystal-reference rule:
{p['eligible_paired_ids_in_score_file']:,} eligible IDs, with paired BiSyRMSD
available for {p['both_bisy_rmsd_available']:,} and paired lDDT-PLI for
{p['both_lddt_pli_available']:,}. The broader topology-matched sensitivity
cohort contains {broad['eligible_paired_ids_in_score_file']:,} eligible IDs,
with {broad['both_bisy_rmsd_available']:,} paired BiSyRMSDs and
{broad['both_lddt_pli_available']:,} paired lDDT-PLI values. Comparisons use each metric's own
paired denominator, not all delivered IDs. Pocket Cα RMSD and QS are omitted
from the primary comparison because both predicted complexes reuse the same
crystal receptor, making receptor-only structural metrics non-independent.

Formal-versus-fixed-frame BiSyRMSD disagreements greater than 0.5 Å occur
in {summary['quality_control']['formal_vs_provisional_rmsd_disagreement_gt_0_5_angstrom']:,}
method rows. An exact-graph, frame-consistent sensitivity subset contains
{frame['eligible_paired_ids_in_score_file']:,} pairs; it is not substituted
for the full prespecified cohort. The fixed-frame direct-coordinate RMSD
sensitivity analysis contains {direct['paired_ids']:,} exact-graph pairs.

These results are an independent crystal-reference evaluation of the
delivered static docking poses, **not** exact reproduction of Lucas's input
preparation. Exact receptor/ligand preparation and the attempted-run cohort
still require his input files or logs. See cohort accounting separately.

The separately labeled direct-coordinate RMSD sensitivity analysis includes
every strict-reference pair for which that exploratory graph-based metric is
available, even when OpenStructure returns no score. It is not substituted
for BiSyRMSD in the primary comparison.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--crystal-audit", type=Path, default=Path("misato_output/crystal_pose_all.csv"))
    parser.add_argument("--identity", type=Path, default=Path("misato_output/inventory/primary_ligand_identity.csv"))
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--allow-partial", action="store_true", help="Only for smoke tests")
    args = parser.parse_args()
    if args.out_dir.exists():
        parser.error("--out-dir already exists; choose a new versioned directory")
    scores = read_csv(args.scores)
    audit = read_csv(args.crystal_audit)
    summary, pairs = summarize(scores, audit, args.allow_partial, read_csv(args.identity))
    args.out_dir.mkdir(parents=True)
    (args.out_dir / "score_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (args.out_dir / "METHODS_RESULTS.md").write_text(methods_text(summary), encoding="utf-8")
    with (args.out_dir / "strict_paired_scores.csv").open("x", newline="", encoding="utf-8") as handle:
        fieldnames = ["target_id", "reference_residue", "prediction_identity_disposition",
                      "diffdock_bisy_rmsd_angstrom", "equibind_bisy_rmsd_angstrom",
                      "diffdock_direct_rmsd_angstrom_provisional",
                      "equibind_direct_rmsd_angstrom_provisional",
                      "diffdock_formal_minus_direct_angstrom", "equibind_formal_minus_direct_angstrom",
                      "diffdock_lddt_pli", "equibind_lddt_pli",
                      "diffdock_status", "equibind_status"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(pairs)
    with (args.out_dir / "rmsd_qc_disagreements.csv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["target_id", "method", "formal_angstrom",
                                                     "provisional_angstrom",
                                                     "absolute_difference_angstrom"])
        writer.writeheader()
        writer.writerows(rmsd_disagreements(scores, audit))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
