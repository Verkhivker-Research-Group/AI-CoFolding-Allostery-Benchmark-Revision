"""Fetch RCSB asymmetric-unit crystal CIFs for delivered MISATO poses.

This is separate from the published benchmark's reference cache. Existing
crystal files are never replaced, and an existing report is never overwritten.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import time
from pathlib import Path

import requests


PDB_ID = re.compile(r"^[A-Z0-9]{4}")
FIELDS = ["pdb_id", "status", "path", "bytes", "sha256", "http_status", "note"]


def requested_ids(manifest: Path) -> list[str]:
    ids: set[str] = set()
    with manifest.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"target_id", "equibind_status", "diffdock_status"}
        if required - set(reader.fieldnames or []):
            raise ValueError(f"{manifest} is not a MISATO target manifest")
        for row in reader:
            if row["equibind_status"] != "candidate" and row["diffdock_status"] != "candidate":
                continue
            match = PDB_ID.match(row["target_id"].strip().upper())
            if not match:
                raise ValueError(f"Invalid target ID: {row['target_id']!r}")
            ids.add(match.group(0))
    return sorted(ids)


def valid_cif(data: bytes, pdb_id: str) -> bool:
    if len(data) < 1000:
        return False
    # RCSB's entry CIF starts with a data block and has atom-site records.
    start = data[:200].lstrip()
    return start.startswith(b"data_") and b"_atom_site." in data


def fetch_one(pdb_id: str, directory: Path, session: requests.Session, timeout: int) -> dict[str, object]:
    destination = directory / f"{pdb_id}.cif"
    if destination.exists():
        data = destination.read_bytes()
        status = "cached_valid" if valid_cif(data, pdb_id) else "cached_invalid_needs_review"
        return {"pdb_id": pdb_id, "status": status, "path": str(destination),
                "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), "http_status": "", "note": ""}
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        response = session.get(url, timeout=timeout)
    except requests.RequestException as error:
        return {"pdb_id": pdb_id, "status": "network_error", "path": "", "bytes": "",
                "sha256": "", "http_status": "", "note": f"{type(error).__name__}: {error}"}
    if response.status_code != 200:
        return {"pdb_id": pdb_id, "status": "http_error", "path": "", "bytes": "",
                "sha256": "", "http_status": response.status_code, "note": ""}
    data = response.content
    if not valid_cif(data, pdb_id):
        return {"pdb_id": pdb_id, "status": "invalid_response", "path": "", "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(), "http_status": 200, "note": "Not a plausible RCSB coordinate CIF"}
    # Exclusive creation protects existing files if another process wins a race.
    try:
        with destination.open("xb") as handle:
            handle.write(data)
    except FileExistsError:
        return fetch_one(pdb_id, directory, session, timeout)
    return {"pdb_id": pdb_id, "status": "downloaded", "path": str(destination),
            "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), "http_status": 200, "note": ""}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", type=Path, default=Path("misato_output/inventory/target_manifest.csv"))
    parser.add_argument("--cache-dir", type=Path, default=Path("misato_output/rcsb_asymmetric_unit"))
    parser.add_argument("--report", type=Path, default=Path("misato_output/rcsb_asymmetric_unit/fetch_report.csv"))
    parser.add_argument("--limit", type=int, default=0, help="First N unique PDB IDs; 0 means all")
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--sleep", type=float, default=0.05)
    args = parser.parse_args()
    if args.limit < 0 or args.timeout <= 0 or args.sleep < 0:
        parser.error("--limit/--sleep must be nonnegative and --timeout positive")
    if args.report.exists():
        parser.error(f"Report already exists; choose a new --report: {args.report}")
    ids = requested_ids(args.target_manifest)
    if args.limit:
        ids = ids[:args.limit]
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with requests.Session() as session, args.report.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for index, pdb_id in enumerate(ids, 1):
            row = fetch_one(pdb_id, args.cache_dir, session, args.timeout)
            writer.writerow(row)
            handle.flush()
            print(f"[{index}/{len(ids)}] {pdb_id}: {row['status']}", flush=True)
            if row["status"] != "cached_valid" and args.sleep:
                time.sleep(args.sleep)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
