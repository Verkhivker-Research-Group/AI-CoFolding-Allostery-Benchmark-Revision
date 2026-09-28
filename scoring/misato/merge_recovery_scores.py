"""Merge MISATO recovery scores while preserving every v3-scored value.

Formal OST metrics and exploratory fixed-frame fallback RMSD occupy different
columns, so downstream analyses cannot accidentally pool them.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

FIELDS = [
    "target_id", "method", "score_status", "pose_rmsd_angstrom", "pose_rmsd_method",
    "lddt_pli", "direct_rmsd_fallback_angstrom", "metric_source",
    "reference_kind", "assigned_crystal_copy", "n_reference_copies",
    "original_v3_status", "allcopies_status", "recovery_status", "recovery_retry_status", "reason",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def index(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    result = {(row["target_id"], row["method"]): row for row in rows}
    if len(result) != len(rows):
        raise ValueError("Duplicate target/method rows")
    return result


def merge(v3: list[dict[str, str]], allcopies: list[dict[str, str]],
          recovery: list[dict[str, str]], retry: list[dict[str, str]],
          audit: list[dict[str, str]], diagnosis: list[dict[str, str]],
          expected_count: int = 16984,
          recovery_retry: list[dict[str, str]] | None = None) -> list[dict[str, str]]:
    old, broad, multi, retried, provisional, reasons = map(
        index, (v3, allcopies, recovery, retry, audit, diagnosis))
    retry_multi = index(recovery_retry or [])
    if len(old) != expected_count or set(old) != set(broad):
        raise ValueError(f"Expected identical {expected_count:,}-row delivered cohorts")
    if set(multi) - set(old) or set(retried) - set(old):
        raise ValueError("Recovery/retry contains undelivered keys")
    expected_recovery = {key for key, row in provisional.items()
                         if row.get("graph_status") == "exact_full_graph"}
    if set(multi) != expected_recovery:
        raise ValueError(f"Recovery rows incomplete: expected {len(expected_recovery)}, "
                         f"observed {len(multi)}")
    if set(retry_multi) - set(multi):
        raise ValueError("Recovery retry contains non-recovery keys")
    result = []
    for key in sorted(old):
        prior = old[key]
        allcopy = broad[key]
        recovered = multi.get(key, {})
        recovered_retry = retry_multi.get(key, {})
        retry_row = retried.get(key, {})
        entry = {field: "" for field in FIELDS}
        entry.update({"target_id": key[0], "method": key[1],
                      "original_v3_status": prior["status"],
                      "allcopies_status": allcopy["status"],
                      "recovery_status": recovered.get("status", ""),
                      "recovery_retry_status": recovered_retry.get("status", "")})
        chosen = None
        if prior["status"] == "scored_both":
            chosen = prior
            entry["metric_source"] = "v3_strict_original"
            entry["pose_rmsd_method"] = "OST_BiSyRMSD"
        elif allcopy["status"] == "scored_both":
            chosen = allcopy
            entry["metric_source"] = "allcopies_v1"
            entry["pose_rmsd_method"] = "OST_BiSyRMSD"
            entry["assigned_crystal_copy"] = allcopy["assigned_ref_ost_chain"]
            entry["n_reference_copies"] = allcopy["n_target_copies"]
        elif retry_row.get("status") == "scored_both":
            chosen = retry_row
            entry["metric_source"] = "allcopies_retry_300s"
            entry["pose_rmsd_method"] = "OST_BiSyRMSD"
            entry["assigned_crystal_copy"] = retry_row["assigned_ref_ost_chain"]
            entry["n_reference_copies"] = retry_row["n_target_copies"]
        elif recovered.get("status") == "scored_both" or recovered_retry.get("status") == "scored_both":
            chosen = recovered if recovered.get("status") == "scored_both" else recovered_retry
            entry["metric_source"] = ("ost_pseudo_multiresidue" if chosen is recovered
                                      else "ost_pseudo_multiresidue_retry")
            entry["pose_rmsd_method"] = "OST_BiSyRMSD"
            entry["reference_kind"] = chosen["candidate_kind"]
            entry["assigned_crystal_copy"] = chosen["assigned_copy_label"]
            entry["n_reference_copies"] = chosen["n_target_copies"]
        if chosen is not None:
            entry["score_status"] = "scored_both"
            entry["pose_rmsd_angstrom"] = chosen["bisy_rmsd_angstrom"]
            entry["lddt_pli"] = chosen["lddt_pli"]
        else:
            direct = provisional.get(key, {}).get("direct_rmsd_angstrom", "")
            if direct:
                entry["score_status"] = "direct_only_fallback"
                entry["direct_rmsd_fallback_angstrom"] = direct
                entry["metric_source"] = "untyped_graph_fixed_frame_exploratory"
                entry["reference_kind"] = provisional[key]["candidate_kind"]
                entry["assigned_crystal_copy"] = provisional[key]["candidate_label"]
                entry["reason"] = recovered.get("error", "OST pseudo-residue metric unavailable")
            else:
                entry["score_status"] = allcopy["status"]
                entry["reason"] = allcopy.get("error", "") or allcopy.get("preaudit_status", "")
                if recovered.get("error"):
                    entry["reason"] = recovered["error"]
            if key in reasons:
                entry["reason"] = (f"BiSyRMSD={reasons[key]['bisy_reason']}; "
                                   f"lDDT-PLI={reasons[key]['lddt_reason']}")
        result.append(entry)
    for row in result:
        key = (row["target_id"], row["method"])
        if old[key]["status"] == "scored_both" and (
                row["pose_rmsd_angstrom"] != old[key]["bisy_rmsd_angstrom"] or
                row["lddt_pli"] != old[key]["lddt_pli"]):
            raise AssertionError(f"v3 score changed for {key}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v3", type=Path,
                        default=Path("misato_output/misato_score_all_heavy_selected_v3.csv"))
    parser.add_argument("--allcopies", type=Path,
                        default=Path("misato_output/misato_score_allcopies_v1.csv"))
    parser.add_argument("--recovery", type=Path, required=True)
    parser.add_argument("--recovery-retry", type=Path)
    parser.add_argument("--retry", type=Path,
                        default=Path("misato_output/retry_timeouts_v1.csv"))
    parser.add_argument("--audit", type=Path,
                        default=Path("misato_output/multiresidue_recovery_all_v2.csv"))
    parser.add_argument("--diagnosis", type=Path,
                        default=Path("misato_output/allcopies_residual_diagnosis_v1.csv"))
    parser.add_argument("--out-csv", type=Path, required=True)
    args = parser.parse_args()
    if args.out_csv.exists():
        parser.error("Output already exists")
    rows = merge(*(read_csv(path) for path in
                   (args.v3, args.allcopies, args.recovery, args.retry,
                    args.audit, args.diagnosis)),
                 recovery_retry=read_csv(args.recovery_retry) if args.recovery_retry else None)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} delivered poses: {dict(Counter(row['score_status'] for row in rows))}")
    print(f"Sources: {dict(Counter(row['metric_source'] for row in rows))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
