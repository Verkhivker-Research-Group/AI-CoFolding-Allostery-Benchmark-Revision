"""Export only applicable MISATO docking metrics, with no placeholder columns."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


HEADER = ["id", "pose rmsd", "qs score", "lddt-pli", "confidence", "rescore_confidence"]
EXPECTED = {"diffdock": 8419, "equibind": 8565}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=Path("misato_output/misato_static_metrics_v1.csv"))
    parser.add_argument("--out-dir", type=Path,
                        default=Path("evalspreadsheets/misato_v4"))
    args = parser.parse_args()
    with args.source.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    counts = Counter(row["method"] for row in rows)
    if dict(counts) != EXPECTED:
        parser.error(f"expected {EXPECTED}, got {dict(counts)}")
    if len({(row["target_id"].lower(), row["method"]) for row in rows}) != len(rows):
        parser.error("duplicate target/method rows")
    destinations = {method: args.out_dir / f"{method}mainligand.csv" for method in EXPECTED}
    if any(path.exists() for path in destinations.values()):
        parser.error("v4 output already exists; refusing to overwrite")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for method, path in destinations.items():
        subset = sorted((row for row in rows if row["method"] == method),
                        key=lambda row: row["target_id"].lower())
        with path.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(HEADER)
            for row in subset:
                writer.writerow([row["target_id"].lower(), row["pose_rmsd_angstrom"],
                                 row["static-receptor QS"], row["lddt_pli"],
                                 row["confidence"], row["rescore_confidence"]])
        print(f"{path}: {len(subset)} rows")


if __name__ == "__main__":
    main()
