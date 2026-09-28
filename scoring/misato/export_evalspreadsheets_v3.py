"""Export the completed MISATO static-docking metrics in legacy CSV layout."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


EXPECTED = {"diffdock": 8419, "equibind": 8565}
METRIC_HEADER = ["id", "pose rmsd", "pocket rmsd", "qs score", "lddt-pli",
                 "confidence", "rescore_confidence", "receptor_mode", "qs metric"]
PROVENANCE_HEADER = ["id", "method", "score_status", "reason", "metric_source",
                     "pose_rmsd_method", "reference_kind", "assigned_crystal_copy",
                     "n_reference_copies", "original_v3_status", "allcopies_status",
                     "recovery_status", "recovery_retry_status", "confidence_source",
                     "rescore_confidence_source", "gnina_status", "gnina_error",
                     "pocket_rmsd_reason", "direct_rmsd_fallback_angstrom"]


def export(rows: list[dict[str, str]], out_dir: Path) -> dict[str, int]:
    counts = Counter(row["method"] for row in rows)
    if dict(counts) != EXPECTED:
        raise ValueError(f"expected complete delivery {EXPECTED}; got {dict(counts)}")
    if len({(row["target_id"].lower(), row["method"]) for row in rows}) != len(rows):
        raise ValueError("duplicate target/method rows")
    required = {"target_id", "method", "score_status", "pose_rmsd_angstrom",
                "lddt_pli", "static-receptor QS", "confidence", "rescore_confidence",
                "receptor_mode", "pocket_rmsd_angstrom", *PROVENANCE_HEADER[2:]}
    if not rows or required - rows[0].keys():
        raise ValueError(f"missing source columns: {sorted(required - rows[0].keys()) if rows else 'empty CSV'}")
    paths = [out_dir / f"{method}mainligand.csv" for method in EXPECTED]
    paths.append(out_dir / "pose_provenance.csv")
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(f"refusing to overwrite: {existing}")
    for row in rows:
        formal = row["score_status"] == "scored_both"
        if row["receptor_mode"] != "static_crystal" or row["pocket_rmsd_angstrom"]:
            raise ValueError(f"unexpected receptor/pocket RMSD for {row['target_id']}")
        if formal != bool(row["pose_rmsd_angstrom"] and row["lddt_pli"] and row["static-receptor QS"]):
            raise ValueError(f"inconsistent formal metric coverage for {row['target_id']}")
        if row["method"] == "equibind" and row["confidence"]:
            raise ValueError("EquiBind native confidence must remain blank")
    out_dir.mkdir(parents=True, exist_ok=True)
    for method in EXPECTED:
        subset = sorted((row for row in rows if row["method"] == method),
                        key=lambda row: row["target_id"].lower())
        with (out_dir / f"{method}mainligand.csv").open("x", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(METRIC_HEADER)
            for row in subset:
                writer.writerow([row["target_id"].lower(), row["pose_rmsd_angstrom"], "",
                                 row["static-receptor QS"], row["lddt_pli"], row["confidence"],
                                 row["rescore_confidence"], row["receptor_mode"],
                                 "static-receptor QS"])
    with (out_dir / "pose_provenance.csv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PROVENANCE_HEADER)
        writer.writeheader()
        for row in sorted(rows, key=lambda row: (row["target_id"].lower(), row["method"])):
            writer.writerow({"id": row["target_id"].lower(), "method": row["method"],
                             **{field: row[field] for field in PROVENANCE_HEADER[2:]}})
    return dict(counts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=Path("misato_output/misato_static_metrics_v1.csv"))
    parser.add_argument("--out-dir", type=Path,
                        default=Path("evalspreadsheets/misato_v3"))
    args = parser.parse_args()
    with args.source.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    print(export(rows, args.out_dir))


if __name__ == "__main__":
    main()
