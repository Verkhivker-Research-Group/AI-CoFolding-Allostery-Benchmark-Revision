"""Read heavy-element signatures for delivered poses absent from paired audit."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from rdkit import Chem, RDLogger


RDLogger.DisableLog("rdApp.*")
FIELDS = ["target_id", "method", "pose_path", "status", "element_signature", "heavy_atoms", "error"]


def signature(path: Path) -> tuple[str, int]:
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
    mol = next((item for item in supplier if item is not None and item.GetNumConformers()), None)
    if mol is None:
        raise ValueError("unreadable_sdf")
    elements = sorted(atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    return ".".join(map(str, elements)), len(elements)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("misato_output/inventory/target_manifest.csv"))
    parser.add_argument("--out-csv", type=Path, required=True)
    args = parser.parse_args()
    if args.out_csv.exists():
        parser.error("--out-csv exists; choose a new versioned filename")
    with args.manifest.open(newline="", encoding="utf-8-sig") as handle:
        manifest = list(csv.DictReader(handle))
    rows = []
    for target in manifest:
        present = [method for method in ("diffdock", "equibind")
                   if target[f"{method}_status"] == "candidate"]
        if len(present) != 1:
            continue
        method = present[0]
        path = Path(target[f"{method}_primary_path"])
        row = {"target_id": target["target_id"], "method": method, "pose_path": str(path),
               "status": "", "element_signature": "", "heavy_atoms": "", "error": ""}
        try:
            row["element_signature"], row["heavy_atoms"] = signature(path)
            row["status"] = "ok"
        except Exception as error:
            row["status"] = "failed"
            row["error"] = f"{type(error).__name__}: {error}"
        rows.append(row)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Unpaired primary poses: {len(rows)}; readable: {sum(row['status'] == 'ok' for row in rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
