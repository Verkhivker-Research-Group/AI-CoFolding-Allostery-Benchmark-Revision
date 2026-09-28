"""Create a provenance-preserving manifest for supplied MISATO docking outputs.

Input folders are read-only evidence. This script neither modifies nor copies them.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


DIFFDOCK_CONFIDENCE = re.compile(r"^rank(?P<rank>\d+)_confidence(?P<confidence>-?\d+(?:\.\d+)?)\.sdf$", re.I)
DIFFDOCK_RANK = re.compile(r"^rank(?P<rank>\d+)\.sdf$", re.I)
EQUIBIND_PRIMARY = "lig_equibind_corrected.sdf"
RAW_FIELDS = ["method", "target_id", "record_type", "is_primary", "rank", "confidence", "source_path", "file_name", "file_size_bytes", "sha256", "status", "note"]
TARGET_FIELDS = ["target_id", "equibind_status", "equibind_primary_path", "diffdock_status", "diffdock_primary_path", "diffdock_candidate_sdf_count", "diffdock_noncanonical_sdf_count", "in_shared_id_set", "in_valid_primary_cohort", "analysis_disposition"]


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def directories(root: Path) -> list[Path]:
    if not root.is_dir():
        raise FileNotFoundError(f"Result root does not exist: {root}")
    return sorted((path for path in root.iterdir() if path.is_dir()), key=lambda path: path.name.upper())


def missing(method: str, target_id: str, status: str, note: str) -> dict[str, object]:
    return {"method": method, "target_id": target_id, "record_type": "missing", "is_primary": False, "rank": "", "confidence": "", "source_path": "", "file_name": "", "file_size_bytes": "", "sha256": "", "status": status, "note": note}


def scan_equibind(root: Path) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    records, targets = [], {}
    for directory in directories(root):
        target = directory.name.upper()
        sdfs = sorted(directory.glob("*.sdf"))
        primary = directory / EQUIBIND_PRIMARY
        if not primary.is_file():
            status = "missing_primary" if not sdfs else "unexpected_layout"
            records.append(missing("equibind", target, status, "Expected lig_equibind_corrected.sdf"))
            targets[target] = {"status": status, "primary": ""}
            continue
        status = "candidate" if len(sdfs) == 1 else "unexpected_layout"
        records.append({"method": "equibind", "target_id": target, "record_type": "prediction", "is_primary": True, "rank": 1, "confidence": "", "source_path": str(primary), "file_name": primary.name, "file_size_bytes": primary.stat().st_size, "sha256": digest(primary), "status": status, "note": "" if status == "candidate" else f"{len(sdfs)} SDF files; expected one"})
        targets[target] = {"status": status, "primary": str(primary)}
    return records, targets


def scan_diffdock(root: Path) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    records, targets = [], {}
    for directory in directories(root):
        target = directory.name.upper()
        sdfs = sorted(directory.glob("*.sdf"))
        primary = directory / "rank1.sdf"
        if not primary.is_file():
            records.append(missing("diffdock", target, "missing_primary", "Expected rank1.sdf"))
            targets[target] = {"status": "missing_primary", "primary": "", "count": len(sdfs), "noncanonical": len(sdfs)}
            continue
        noncanonical = 0
        for sdf in sdfs:
            confidence = DIFFDOCK_CONFIDENCE.match(sdf.name)
            rank = DIFFDOCK_RANK.match(sdf.name)
            primary_file = sdf.name.lower() == "rank1.sdf"
            if primary_file:
                record_type, model_rank, score, status, note = "prediction", 1, "", "candidate", "Canonical delivered top-ranked pose"
            elif confidence:
                record_type, model_rank, score, status, note = "ranked_candidate", int(confidence.group("rank")), confidence.group("confidence"), "retained_unselected", "May be a repeat-run candidate; never pooled automatically"
                noncanonical += 1
            elif rank:
                record_type, model_rank, score, status, note = "ranked_candidate", int(rank.group("rank")), "", "retained_unselected", "Non-primary canonical rank"
                noncanonical += 1
            else:
                record_type, model_rank, score, status, note = "unrecognized", "", "", "unexpected_layout", "Filename does not match known DiffDock output pattern"
                noncanonical += 1
            records.append({"method": "diffdock", "target_id": target, "record_type": record_type, "is_primary": primary_file, "rank": model_rank, "confidence": score, "source_path": str(sdf), "file_name": sdf.name, "file_size_bytes": sdf.stat().st_size, "sha256": digest(sdf), "status": status, "note": note})
        targets[target] = {"status": "candidate", "primary": str(primary), "count": len(sdfs), "noncanonical": noncanonical}
    return records, targets


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def target_manifest(equibind: dict[str, dict[str, object]], diffdock: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for target in sorted(set(equibind) | set(diffdock)):
        eq = equibind.get(target, {"status": "not_delivered", "primary": ""})
        dd = diffdock.get(target, {"status": "not_delivered", "primary": "", "count": 0, "noncanonical": 0})
        shared = target in equibind and target in diffdock
        valid = eq["status"] == "candidate" and dd["status"] == "candidate"
        rows.append({"target_id": target, "equibind_status": eq["status"], "equibind_primary_path": eq["primary"], "diffdock_status": dd["status"], "diffdock_primary_path": dd["primary"], "diffdock_candidate_sdf_count": dd.get("count", 0), "diffdock_noncanonical_sdf_count": dd.get("noncanonical", 0), "in_shared_id_set": shared, "in_valid_primary_cohort": valid, "analysis_disposition": "validated_primary_candidate" if valid else "needs_review_or_missing_primary"})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equibind-root", type=Path, required=True)
    parser.add_argument("--diffdock-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    eq_records, eq_targets = scan_equibind(args.equibind_root)
    dd_records, dd_targets = scan_diffdock(args.diffdock_root)
    records, targets = eq_records + dd_records, target_manifest(eq_targets, dd_targets)
    write_csv(args.out_dir / "raw_prediction_records.csv", RAW_FIELDS, records)
    write_csv(args.out_dir / "target_manifest.csv", TARGET_FIELDS, targets)
    summary = {"created_utc": datetime.now(timezone.utc).isoformat(), "sources": {"equibind_root": str(args.equibind_root), "diffdock_root": str(args.diffdock_root)}, "equibind": {"target_directories": len(eq_targets), "record_statuses": dict(sorted(Counter(row["status"] for row in eq_records).items()))}, "diffdock": {"target_directories": len(dd_targets), "record_statuses": dict(sorted(Counter(row["status"] for row in dd_records).items()))}, "cohort": {"all_target_ids": len(targets), "shared_target_ids": sum(bool(row["in_shared_id_set"]) for row in targets), "valid_primary_candidates": sum(bool(row["in_valid_primary_cohort"]) for row in targets)}, "policy": "Only EquiBind lig_equibind_corrected.sdf and DiffDock rank1.sdf are primary candidates. Other DiffDock files are retained but not pooled."}
    (args.out_dir / "inventory_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
