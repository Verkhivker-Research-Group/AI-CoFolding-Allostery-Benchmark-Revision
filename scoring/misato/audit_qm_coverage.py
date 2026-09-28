"""Verify public MISATO QM.hdf5 and compare its IDs with delivered poses."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


PUBLISHED_QM_MD5 = "1199f0af1eac684da3a6c2ddd0f321df"
PUBLISHED_SOURCE = "https://zenodo.org/records/7711953/files/QM.hdf5"


def md5_file(path: Path) -> str:
    digest = hashlib.md5()  # Published Zenodo checksum; used only for integrity.
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def audit(qm_path: Path, manifest: list[dict[str, str]],
          identity: list[dict[str, str]] | None = None,
          unpaired: list[dict[str, str]] | None = None) -> tuple[dict, list[dict]]:
    import h5py

    checksum = md5_file(qm_path)
    if checksum != PUBLISHED_QM_MD5:
        raise ValueError(f"QM.hdf5 MD5 mismatch: {checksum} != {PUBLISHED_QM_MD5}")
    manifest_by_id = {row["target_id"].upper(): row for row in manifest}
    if len(manifest_by_id) != len(manifest):
        raise ValueError("Manifest contains duplicate target IDs")
    identity_by_id = {row["target_id"].upper(): row for row in (identity or [])}
    if identity is not None and len(identity_by_id) != len(identity):
        raise ValueError("Identity audit contains duplicate IDs")
    unpaired_by_key = {(row["target_id"].upper(), row["method"]): row
                       for row in (unpaired or [])}
    if unpaired is not None and len(unpaired_by_key) != len(unpaired):
        raise ValueError("Unpaired signature audit contains duplicate target/method pairs")
    qm_signatures: dict[str, str] = {}
    qm_signature_errors: Counter[str] = Counter()
    with h5py.File(qm_path, "r") as source:
        source_ids = {key.upper() for key in source.keys() if isinstance(source[key], h5py.Group)}
        root_items = len(source)
        for target_id in manifest_by_id.keys() & source_ids:
            try:
                numbers = [int(item.decode("ascii"))
                           for item in source[target_id]["atom_properties"]["atom_names"][:]]
                qm_signatures[target_id] = ".".join(str(number) for number in sorted(
                    number for number in numbers if number > 1))
            except (KeyError, ValueError, UnicodeDecodeError, TypeError) as error:
                qm_signature_errors[type(error).__name__] += 1
    if len(source_ids) != root_items:
        raise ValueError("QM root contains non-group keys or case-colliding IDs")
    delivered = set(manifest_by_id)
    pose = {method: {target_id for target_id, row in manifest_by_id.items()
                     if row[f"{method}_status"] == "candidate"}
            for method in ("diffdock", "equibind")}
    both = pose["diffdock"] & pose["equibind"]
    rows = []
    chemical_counts = {method: Counter() for method in pose}
    for target_id in sorted(source_ids | delivered):
        row = manifest_by_id.get(target_id, {})
        qm_signature = qm_signatures.get(target_id, "")
        chemical = {}
        for method in pose:
            if target_id not in pose[method]:
                state = "no_delivered_primary_pose"
            elif not qm_signature:
                state = "qm_atom_signature_unavailable"
            else:
                pose_signature = (identity_by_id.get(target_id, {}).get(f"{method}_element_signature", "")
                                  or unpaired_by_key.get((target_id, method), {}).get("element_signature", ""))
                state = ("match" if pose_signature == qm_signature else
                         "pose_signature_unavailable" if not pose_signature else "mismatch")
            chemical[f"{method}_element_vs_qm"] = state
            if target_id in pose[method]:
                chemical_counts[method][state] += 1
        rows.append({"target_id": target_id, "in_qm_hdf5": target_id in source_ids,
                     "delivered_directory": target_id in delivered,
                     "diffdock_primary_pose": target_id in pose["diffdock"],
                     "equibind_primary_pose": target_id in pose["equibind"],
                     "paired_primary_poses": target_id in both,
                     "qm_heavy_element_signature": qm_signature,
                     "qm_heavy_atom_count": len(qm_signature.split(".")) if qm_signature else "",
                     **chemical,
                     "equibind_manifest_status": row.get("equibind_status", "not_delivered"),
                     "diffdock_manifest_status": row.get("diffdock_status", "not_delivered")})
    summary = {
        "qm_path": str(qm_path.resolve()), "published_source": PUBLISHED_SOURCE,
        "published_md5": PUBLISHED_QM_MD5, "verified_md5": checksum,
        "qm_group_ids": len(source_ids), "qm_root_items": root_items,
        "delivered_unique_ids": len(delivered),
        "delivered_ids_in_qm": len(delivered & source_ids),
        "delivered_ids_outside_qm": len(delivered - source_ids),
        "qm_ids_without_delivered_directory": len(source_ids - delivered),
        "diffdock_primary_in_qm": len(pose["diffdock"] & source_ids),
        "equibind_primary_in_qm": len(pose["equibind"] & source_ids),
        "paired_primary_in_qm": len(both & source_ids),
        "diffdock_primary_outside_qm": len(pose["diffdock"] - source_ids),
        "equibind_primary_outside_qm": len(pose["equibind"] - source_ids),
        "paired_primary_outside_qm": len(both - source_ids),
        "qm_atom_signature_errors": dict(qm_signature_errors),
        "pose_heavy_element_match_qm_by_method": {
            method: dict(sorted(counts.items())) for method, counts in chemical_counts.items()},
        "interpretation": (
            "Verifies source-set membership and delivered successful-pose coverage. "
            "Does not reveal which missing QM IDs were submitted, failed, or never attempted."
        ),
    }
    return summary, rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qm", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=Path("misato_output/inventory/target_manifest.csv"))
    parser.add_argument("--identity", type=Path, default=Path("misato_output/inventory/primary_ligand_identity.csv"))
    parser.add_argument("--unpaired-signatures", type=Path,
                        default=Path("misato_output/inventory/unpaired_pose_signatures_v1.csv"))
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.exists():
        parser.error("--out-dir exists; choose a new versioned directory")
    with args.manifest.open(newline="", encoding="utf-8-sig") as handle:
        manifest = list(csv.DictReader(handle))
    with args.identity.open(newline="", encoding="utf-8-sig") as handle:
        identity = list(csv.DictReader(handle))
    if args.unpaired_signatures.exists():
        with args.unpaired_signatures.open(newline="", encoding="utf-8-sig") as handle:
            unpaired = list(csv.DictReader(handle))
    else:
        unpaired = []
    summary, rows = audit(args.qm, manifest, identity, unpaired)
    args.out_dir.mkdir(parents=True)
    (args.out_dir / "qm_coverage_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    with (args.out_dir / "qm_id_coverage.csv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
