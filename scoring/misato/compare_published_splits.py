"""Compare delivered docking targets with MISATO's published MD split IDs.

This is an MD-split overlap audit, not a reconstruction of Lucas's QM-based
input filtering, an MD-frame analysis, or a docking-accuracy analysis.
Only the three small text lists are downloaded; MD.hdf5 is never accessed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from urllib.request import urlopen


RECORD = "https://zenodo.org/records/7711953"
SPLITS = {
    "train": ("train_MD.txt", "0eff7e1729307f9d773353382fbb1f1e"),
    "val": ("val_MD.txt", "a37f0dfead91accfe089b3555e0d9e33"),
    "test": ("test_MD.txt", "a3f8ada0c6562ff3d8e53653816a36ba"),
}
# Lucas's folders include seven suffixed IDs (for example, 1G42_A).
# Retain those as distinct, auditable IDs; never silently strip the suffix.
ID_PATTERN = re.compile(r"^[A-Z0-9]{4}[A-Z0-9_]{0,12}$")
FIELDS = [
    "target_id", "published_md_split", "equibind_delivered", "diffdock_delivered",
    "equibind_primary_candidate", "diffdock_primary_candidate",
    "paired_primary_candidate", "coverage_category",
]
SUMMARY_FIELDS = [
    "split", "published_ids", "equibind_delivered", "diffdock_delivered",
    "both_delivered", "neither_delivered", "equibind_primary_candidates",
    "diffdock_primary_candidates", "paired_primary_candidates",
]


def checked_list(data: bytes, name: str, expected_md5: str) -> set[str]:
    actual_md5 = hashlib.md5(data).hexdigest()
    if actual_md5 != expected_md5:
        raise ValueError(f"{name}: MD5 {actual_md5} differs from published {expected_md5}")
    ids: set[str] = set()
    for line_number, raw_line in enumerate(data.decode("utf-8-sig").splitlines(), 1):
        target_id = raw_line.strip().upper()
        if not target_id:
            continue
        if not ID_PATTERN.fullmatch(target_id):
            raise ValueError(f"{name}:{line_number}: invalid target ID {target_id!r}")
        if target_id in ids:
            raise ValueError(f"{name}:{line_number}: duplicate target ID {target_id}")
        ids.add(target_id)
    return ids


def load_splits(split_dir: Path | None, out_dir: Path) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    membership: dict[str, str] = {}
    sources: dict[str, dict[str, str]] = {}
    for split, (filename, expected_md5) in SPLITS.items():
        url = f"{RECORD}/files/{filename}?download=1"
        if split_dir is None:
            with urlopen(url, timeout=60) as response:
                data = response.read()
        else:
            data = (split_dir / filename).read_bytes()
        ids = checked_list(data, filename, expected_md5)
        if split_dir is None:
            source_dir = out_dir / "published_lists"
            source_dir.mkdir(parents=True, exist_ok=True)
            (source_dir / filename).write_bytes(data)
        for target_id in ids:
            if target_id in membership:
                raise ValueError(f"{target_id} occurs in both {membership[target_id]} and {split} splits")
            membership[target_id] = split
        sources[split] = {"url": url, "md5": expected_md5, "ids": len(ids)}
    return membership, sources


def boolean(value: str, column: str, target_id: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"{target_id}: {column} must be True or False, got {value!r}")
    return normalized == "true"


def read_manifest(path: Path) -> dict[str, dict[str, str]]:
    required = {"target_id", "equibind_status", "diffdock_status", "in_valid_primary_cohort"}
    targets: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            target_id = row["target_id"].strip().upper()
            if not ID_PATTERN.fullmatch(target_id):
                raise ValueError(f"{path}: invalid target ID {target_id!r}")
            if target_id in targets:
                raise ValueError(f"{path}: duplicate target ID {target_id}")
            both_valid = row["equibind_status"] == row["diffdock_status"] == "candidate"
            if boolean(row["in_valid_primary_cohort"], "in_valid_primary_cohort", target_id) != both_valid:
                raise ValueError(f"{target_id}: primary-cohort flag disagrees with method statuses")
            targets[target_id] = row
    return targets


def compare(membership: dict[str, str], targets: dict[str, dict[str, str]]) -> list[dict[str, str | bool]]:
    rows = []
    for target_id in sorted(membership.keys() | targets.keys()):
        target = targets.get(target_id, {})
        eq_status = target.get("equibind_status", "not_delivered")
        dd_status = target.get("diffdock_status", "not_delivered")
        eq_delivered = eq_status != "not_delivered"
        dd_delivered = dd_status != "not_delivered"
        eq_candidate = eq_status == "candidate"
        dd_candidate = dd_status == "candidate"
        if target_id not in membership:
            category = "delivered_outside_published_md_splits"
        elif eq_candidate and dd_candidate:
            category = "paired_primary_candidate"
        elif eq_delivered and dd_delivered:
            category = "both_delivered_primary_incomplete"
        elif eq_delivered:
            category = "equibind_only_delivered"
        elif dd_delivered:
            category = "diffdock_only_delivered"
        else:
            category = "neither_delivered"
        rows.append({
            "target_id": target_id,
            "published_md_split": membership.get(target_id, ""),
            "equibind_delivered": eq_delivered,
            "diffdock_delivered": dd_delivered,
            "equibind_primary_candidate": eq_candidate,
            "diffdock_primary_candidate": dd_candidate,
            "paired_primary_candidate": eq_candidate and dd_candidate,
            "coverage_category": category,
        })
    return rows


def summarize(rows: list[dict[str, str | bool]]) -> list[dict[str, str | int]]:
    summaries = []
    for split in SPLITS:
        subset = [row for row in rows if row["published_md_split"] == split]
        summaries.append({
            "split": split,
            "published_ids": len(subset),
            "equibind_delivered": sum(row["equibind_delivered"] is True for row in subset),
            "diffdock_delivered": sum(row["diffdock_delivered"] is True for row in subset),
            "both_delivered": sum(row["equibind_delivered"] is True and row["diffdock_delivered"] is True for row in subset),
            "neither_delivered": sum(row["coverage_category"] == "neither_delivered" for row in subset),
            "equibind_primary_candidates": sum(row["equibind_primary_candidate"] is True for row in subset),
            "diffdock_primary_candidates": sum(row["diffdock_primary_candidate"] is True for row in subset),
            "paired_primary_candidates": sum(row["paired_primary_candidate"] is True for row in subset),
        })
    return summaries


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", type=Path, default=Path("misato_output/inventory/target_manifest.csv"))
    parser.add_argument("--split-dir", type=Path, help="Use locally saved official train_MD.txt, val_MD.txt, and test_MD.txt instead of downloading them")
    parser.add_argument("--out-dir", type=Path, default=Path("misato_output/split_comparison"))
    args = parser.parse_args()
    targets = read_manifest(args.target_manifest)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    membership, sources = load_splits(args.split_dir, args.out_dir)
    rows = compare(membership, targets)
    split_rows = summarize(rows)
    outside_rows = [row for row in rows if not row["published_md_split"]]
    write_csv(args.out_dir / "target_comparison.csv", FIELDS, rows)
    write_csv(args.out_dir / "split_summary.csv", SUMMARY_FIELDS, split_rows)
    summary = {
        "published_record": RECORD,
        "published_lists": sources,
        "delivered_manifest": str(args.target_manifest.resolve()),
        "published_id_count": len(membership),
        "delivered_id_count": len(targets),
        "outside_published_md_splits": {
            "delivered_ids": len(outside_rows),
            "equibind_delivered": sum(row["equibind_delivered"] is True for row in outside_rows),
            "diffdock_delivered": sum(row["diffdock_delivered"] is True for row in outside_rows),
            "both_delivered": sum(row["equibind_delivered"] is True and row["diffdock_delivered"] is True for row in outside_rows),
            "paired_primary_candidates": sum(row["paired_primary_candidate"] is True for row in outside_rows),
        },
        "categories": dict(sorted(Counter(str(row["coverage_category"]) for row in rows).items())),
        "splits": split_rows,
        "interpretation": "MD-split overlap only; Lucas reports that these splits were not used to select his static-receptor/QM-ligand docking inputs.",
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
