"""
00_download_ref_cifs.py
Download reference CIF structures from RCSB for all ASD entries
needed to score the Chai and AF3 runs.

Usage:
    python 00_download_ref_cifs.py

Writes:  asd_score/ref_cifs/<PDBID>.cif  (667 files total)
"""
from __future__ import annotations

import csv
import os
import time
import urllib.request
import urllib.error
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
ASD_MANIFEST   = Path("C:/Users/Ryan/Downloads/ASD_Release_202309_AS.txt")
CHAI_PRED_DIR  = Path("C:/Users/Ryan/Downloads/chai/outputs")
AF3_PRED_DIR   = Path("C:/Users/Ryan/Downloads/af3/output_final")
OUT_DIR        = Path("C:/Users/Ryan/AI-CoFolding-Allostery-Benchmark/asd_score/ref_cifs")

RCSB_CIF_URL   = "https://files.rcsb.org/download/{pdb_id}.cif"
SLEEP_S        = 0.05   # polite rate-limit between requests


def collect_needed_pdb_ids() -> set[str]:
    """
    Chai folders: 1B86_D_DG2 → PDB = 1B86
    AF3 folders:  O00411 (UniProt) → all allosteric_pdb entries in ASD manifest
    """
    # Chai: PDB ID is the first token of each folder name
    chai_pdbs = {
        d.split("_")[0].upper()
        for d in os.listdir(CHAI_PRED_DIR)
        if (CHAI_PRED_DIR / d).is_dir() and "_" in d
    }

    # AF3: map UniProt → allosteric PDB entries via ASD manifest
    af3_uniprots = {
        d for d in os.listdir(AF3_PRED_DIR)
        if (AF3_PRED_DIR / d).is_dir()
    }
    af3_pdbs: set[str] = set()
    with ASD_MANIFEST.open(encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            uid = row.get("pdb_uniprot", "").strip()
            pdb = row.get("allosteric_pdb", "").strip().upper()
            if uid in af3_uniprots and pdb:
                af3_pdbs.add(pdb)

    all_pdbs = chai_pdbs | af3_pdbs
    print(f"Chai PDB IDs : {len(chai_pdbs)}")
    print(f"AF3  PDB IDs : {len(af3_pdbs)}")
    print(f"Total unique : {len(all_pdbs)}")
    return all_pdbs


def download_cif(pdb_id: str, out_dir: Path) -> bool:
    dest = out_dir / f"{pdb_id}.cif"
    if dest.exists() and dest.stat().st_size > 1000:
        return True  # already downloaded
    url = RCSB_CIF_URL.format(pdb_id=pdb_id)
    try:
        urllib.request.urlretrieve(url, dest)
        return True
    except urllib.error.HTTPError as e:
        print(f"  [HTTP {e.code}] {pdb_id}")
        return False
    except Exception as e:
        print(f"  [ERROR] {pdb_id}: {e}")
        return False


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pdb_ids = sorted(collect_needed_pdb_ids())

    # Check what we already have
    existing = {p.stem.upper() for p in OUT_DIR.glob("*.cif")}
    to_download = [p for p in pdb_ids if p not in existing]
    print(f"Already in ref_cifs/ : {len(existing)}")
    print(f"To download          : {len(to_download)}")

    ok = skipped = failed = 0
    for i, pdb_id in enumerate(to_download, 1):
        if download_cif(pdb_id, OUT_DIR):
            ok += 1
            if i % 50 == 0:
                print(f"  [{i}/{len(to_download)}] downloaded {ok} OK, {failed} failed")
        else:
            failed += 1
        time.sleep(SLEEP_S)

    print(f"\nDone. Downloaded: {ok}  Failed: {failed}  Already present: {len(existing)}")
    print(f"Total CIFs in {OUT_DIR}: {len(list(OUT_DIR.glob('*.cif')))}")


if __name__ == "__main__":
    main()
