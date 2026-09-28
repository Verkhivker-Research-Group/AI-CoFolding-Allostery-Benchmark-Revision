"""Export a portable, path-free index of Lucas's delivered primary pose IDs.

The source manifests and complete score CSV stay local under misato_output/.
This smaller index is suitable for tracking in the analysis repository.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "target_id", "diffdock_primary_delivered", "equibind_primary_delivered",
    "diffdock_strict_reference_eligible", "equibind_strict_reference_eligible",
    "prediction_identity_disposition", "diffdock_score_status", "equibind_score_status",
    "diffdock_bisy_rmsd_available", "equibind_bisy_rmsd_available",
    "diffdock_lddt_pli_available", "equibind_lddt_pli_available",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def index_unique(rows: list[dict[str, str]], keys: tuple[str, ...]) -> dict[tuple[str, ...], dict[str, str]]:
    indexed = {tuple(row[key] for key in keys): row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError(f"Duplicate rows by {keys}")
    return indexed


def export(manifest: list[dict[str, str]], audit: list[dict[str, str]],
           scores: list[dict[str, str]], identity: list[dict[str, str]]) -> list[dict[str, str]]:
    manifests = index_unique(manifest, ("target_id",))
    audits = index_unique(audit, ("target_id", "method"))
    scored = index_unique(scores, ("target_id", "method"))
    identities = index_unique(identity, ("target_id",))
    delivered_keys = {(target_id, method) for (target_id,), row in manifests.items()
                      for method in ("diffdock", "equibind")
                      if row[f"{method}_status"] == "candidate"}
    if set(audits) != delivered_keys or set(scored) != delivered_keys:
        raise ValueError("Crystal audit and score file must exactly cover delivered primary poses")
    result = []
    for (target_id,), manifest_row in sorted(manifests.items()):
        item = {"target_id": target_id,
                "prediction_identity_disposition": identities.get((target_id,), {}).get(
                    "disposition", "not_audited")}
        for method in ("diffdock", "equibind"):
            key = (target_id, method)
            present = key in delivered_keys
            item[f"{method}_primary_delivered"] = str(present).lower()
            item[f"{method}_strict_reference_eligible"] = str(
                present and audits[key]["reference_selection"] == "single_candidate"
                and audits[key]["reference_status"] == "provisional_pose_rmsd"
                and bool(audits[key]["reference_residues"])).lower()
            item[f"{method}_score_status"] = scored[key]["status"] if present else "not_delivered"
            item[f"{method}_bisy_rmsd_available"] = str(
                present and bool(scored[key]["bisy_rmsd_angstrom"])).lower()
            item[f"{method}_lddt_pli_available"] = str(
                present and bool(scored[key]["lddt_pli"])).lower()
        result.append(item)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path,
                        default=Path("misato_output/inventory/target_manifest.csv"))
    parser.add_argument("--crystal-audit", type=Path,
                        default=Path("misato_output/crystal_pose_all.csv"))
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--identity", type=Path,
                        default=Path("misato_output/inventory/primary_ligand_identity.csv"))
    parser.add_argument("--out-csv", type=Path, required=True)
    args = parser.parse_args()
    if args.out_csv.exists():
        parser.error("Output already exists; choose a new file")
    rows = export(read_csv(args.manifest), read_csv(args.crystal_audit),
                  read_csv(args.scores), read_csv(args.identity))
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} target IDs to {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
