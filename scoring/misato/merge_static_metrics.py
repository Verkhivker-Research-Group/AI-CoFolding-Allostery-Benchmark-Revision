"""Join formal MISATO scores, static-receptor QS, and gnina CNNscore.

Pocket RMSD is intentionally missing: both methods docked into the unchanged
crystal receptor. No blank metric is converted to zero. This is a versioned
intermediate table, not the step-6 evalspreadsheets export.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

try:
    from .export_evalspreadsheets_v2 import diffdock_confidence
    from .merge_recovery_scores import FIELDS as BASE_FIELDS
except ImportError:
    from export_evalspreadsheets_v2 import diffdock_confidence
    from merge_recovery_scores import FIELDS as BASE_FIELDS


EXTRA_FIELDS = ["static-receptor QS", "confidence", "confidence_source",
                "rescore_confidence", "rescore_confidence_source", "gnina_status",
                "gnina_error", "pocket_rmsd_angstrom", "pocket_rmsd_reason",
                "receptor_mode"]
FIELDS = BASE_FIELDS + EXTRA_FIELDS
POCKET_REASON = "not applicable: rigid-receptor docking into the reference crystal structure"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def unique_index(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    result = {(row["target_id"], row["method"]): row for row in rows}
    if len(result) != len(rows):
        raise ValueError("duplicate target/method rows")
    return result


def merge(scores: list[dict[str, str]], qs_rows: list[dict[str, str]],
          gnina_rows: list[dict[str, str]], records: list[dict[str, str]],
          expected_count: int = 16984) -> list[dict[str, str]]:
    base = unique_index(scores)
    qs = unique_index(qs_rows)
    gnina = unique_index(gnina_rows)
    if len(base) != expected_count:
        raise ValueError(f"expected {expected_count} delivered poses; got {len(base)}")
    formally_scored = {key for key, row in base.items() if row["score_status"] == "scored_both"}
    if set(qs) != formally_scored:
        raise ValueError(f"static QS keys differ: {len(qs)} versus {len(formally_scored)} scored")
    if set(gnina) != set(base):
        raise ValueError(f"gnina keys differ: {len(gnina)} versus {len(base)} delivered")
    confidence = diffdock_confidence(records)
    if {key[0] for key in base if key[1] == "diffdock"} - set(confidence):
        raise ValueError("some delivered DiffDock primary poses lack a SHA-256 confidence match")
    result = []
    for key, row in sorted(base.items()):
        out = {field: row.get(field, "") for field in BASE_FIELDS}
        out.update({field: "" for field in EXTRA_FIELDS})
        out["receptor_mode"] = "static_crystal"
        out["pocket_rmsd_reason"] = POCKET_REASON
        if key[1] == "diffdock":
            out["confidence"] = confidence[key[0]]
            out["confidence_source"] = "DiffDock_rank1_SHA256_matched_logit"
        else:
            out["confidence_source"] = "not_reported_by_EquiBind"
        q = qs.get(key)
        if q:
            if q["status"] != "scored":
                raise ValueError(f"static QS failed for formal score {key}: {q['error']}")
            value = float(q["static_receptor_qs"])
            if not 0 <= value <= 1:
                raise ValueError(f"static QS out of range for {key}")
            out["static-receptor QS"] = q["static_receptor_qs"]
        g = gnina[key]
        out["gnina_status"] = g["status"]
        out["gnina_error"] = g.get("error", "")
        if g["status"] == "scored":
            value = float(g["cnnscore"])
            if not 0 <= value <= 1:
                raise ValueError(f"CNNscore out of range for {key}")
            out["rescore_confidence"] = g["cnnscore"]
            out["rescore_confidence_source"] = f"{g['gnina_version']}_CNNscore_score_only"
        result.append(out)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path,
                        default=Path("misato_output/misato_score_recovered_v1.csv"))
    parser.add_argument("--static-qs", type=Path,
                        default=Path("misato_output/static_qs_all_v1.csv"))
    parser.add_argument("--gnina", type=Path, required=True)
    parser.add_argument("--records", type=Path,
                        default=Path("misato_output/inventory/raw_prediction_records.csv"))
    parser.add_argument("--out-csv", type=Path, required=True)
    args = parser.parse_args()
    if args.out_csv.exists():
        parser.error("Output exists; select a new versioned path")
    rows = merge(read_rows(args.scores), read_rows(args.static_qs),
                 read_rows(args.gnina), read_rows(args.records))
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} delivered rows: formal metrics "
          f"{sum(bool(row['pose_rmsd_angstrom'] and row['lddt_pli']) for row in rows)}, "
          f"static QS {sum(bool(row['static-receptor QS']) for row in rows)}, "
          f"gnina CNNscore {sum(bool(row['rescore_confidence']) for row in rows)}")
    print(f"Gnina status: {dict(Counter(row['gnina_status'] for row in rows))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
