"""Export MISATO scores in the legacy evalspreadsheets column format."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


HEADER = ["id", "pose rmsd", "pocket rmsd", "qs score", "lddt-pli", "confidence"]
METHODS = ("diffdock", "equibind")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path,
                        default=Path("misato_output/misato_score_all_heavy_selected_v3.csv"))
    parser.add_argument("--out-dir", type=Path, default=Path("evalspreadsheets/misato"))
    args = parser.parse_args()
    with args.scores.open(newline="", encoding="utf-8-sig") as handle:
        source = list(csv.DictReader(handle))
    required = {"target_id", "method", "bisy_rmsd_angstrom", "lddt_pli", "status"}
    if not source or not required.issubset(source[0]):
        parser.error("Score CSV is empty or has the wrong schema")
    counts = Counter(row["method"] for row in source)
    if len(source) != 16984 or counts != {"diffdock": 8419, "equibind": 8565}:
        parser.error(f"Expected the complete scored delivery; got {len(source)} rows: {dict(counts)}")
    if len({(row["target_id"], row["method"]) for row in source}) != len(source):
        parser.error("Duplicate target/method rows")
    destinations = {method: args.out_dir / f"{method}mainligand.csv" for method in METHODS}
    existing = [str(path) for path in destinations.values() if path.exists()]
    if existing:
        parser.error(f"Output exists: {existing}; choose a new directory")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for method in METHODS:
        destination = destinations[method]
        rows = sorted((row for row in source if row["method"] == method),
                      key=lambda row: row["target_id"].lower())
        with destination.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(HEADER)
            for row in rows:
                writer.writerow([row["target_id"].lower(), row["bisy_rmsd_angstrom"],
                                 "", "", row["lddt_pli"], ""])
        print(f"{destination}: {len(rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
