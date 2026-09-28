"""Export MISATO scores in the legacy evalspreadsheets format, with confidence.

Changes from ``export_evalspreadsheets.py``:

* Reads the all-copies score CSV (``score_static_docking_allcopies_wsl.py``).
* DiffDock ``confidence`` is recovered from the delivered files. ``rank1.sdf``
  carries no score in its name, but it is a byte-identical copy of exactly one
  ``rank1_confidence-X.sdf`` in the same target folder; the SHA-256 join in
  ``raw_prediction_records.csv`` identifies it. Values are DiffDock's raw
  confidence (a logit: higher is more confident), not a probability.
* EquiBind produces no confidence score, so its column stays blank.
* ``pocket rmsd`` and ``qs score`` stay blank: with a static crystal receptor
  they are not independent of the input structure (see METHODS_RESULTS.md).
  Blank means not measured, never zero.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

HEADER = ["id", "pose rmsd", "pocket rmsd", "qs score", "lddt-pli", "confidence"]
METHODS = ("diffdock", "equibind")
EXPECTED = {"diffdock": 8419, "equibind": 8565}


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def diffdock_confidence(records: list[dict[str, str]]) -> dict[str, str]:
    """target_id -> confidence of the primary rank1.sdf, via exact SHA-256 match."""
    primary = {}
    scored = defaultdict(list)
    for row in records:
        if row["method"] != "diffdock":
            continue
        if row["is_primary"].strip().lower() == "true":
            primary[row["target_id"]] = row["sha256"]
        elif row["confidence"].strip():
            scored[(row["target_id"], row["sha256"])].append(row["confidence"].strip())
    out = {}
    for target, sha in primary.items():
        values = set(scored.get((target, sha), []))
        if len(values) == 1:
            out[target] = values.pop()
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scores", type=Path, required=True,
                        help="Output of score_static_docking_allcopies_wsl.py")
    parser.add_argument("--records", type=Path,
                        default=Path("misato_output/inventory/raw_prediction_records.csv"))
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    source = read(args.scores)
    required = {"target_id", "method", "bisy_rmsd_angstrom", "lddt_pli", "status"}
    if not source or not required.issubset(source[0]):
        parser.error("Score CSV is empty or has the wrong schema")
    counts = Counter(row["method"] for row in source)
    if dict(counts) != EXPECTED:
        parser.error(f"Expected the complete delivery {EXPECTED}; got {dict(counts)}")
    if len({(row["target_id"], row["method"]) for row in source}) != len(source):
        parser.error("Duplicate target/method rows")
    destinations = {method: args.out_dir / f"{method}mainligand.csv" for method in METHODS}
    existing = [str(path) for path in destinations.values() if path.exists()]
    if existing:
        parser.error(f"Output exists: {existing}; choose a new directory")

    confidence = diffdock_confidence(read(args.records))
    missing_conf = EXPECTED["diffdock"] - sum(1 for row in source if row["method"] == "diffdock"
                                              and row["target_id"] in confidence)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for method in METHODS:
        rows = sorted((row for row in source if row["method"] == method),
                      key=lambda row: row["target_id"].lower())
        with destinations[method].open("x", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(HEADER)
            for row in rows:
                conf = confidence.get(row["target_id"], "") if method == "diffdock" else ""
                writer.writerow([row["target_id"].lower(), row["bisy_rmsd_angstrom"], "", "",
                                 row["lddt_pli"], conf])
        n = len(rows)
        cov = lambda key: sum(1 for row in rows if row[key] != "")
        print(f"{destinations[method]}: {n} rows | pose rmsd {cov('bisy_rmsd_angstrom')} "
              f"({100 * cov('bisy_rmsd_angstrom') / n:.1f}%) | lddt-pli {cov('lddt_pli')} "
              f"({100 * cov('lddt_pli') / n:.1f}%)")
    print(f"DiffDock confidence recovered for {EXPECTED['diffdock'] - missing_conf}/"
          f"{EXPECTED['diffdock']} poses")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
