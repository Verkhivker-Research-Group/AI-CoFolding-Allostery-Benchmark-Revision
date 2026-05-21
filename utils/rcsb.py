"""RCSB helper functions.

Small wrappers for downloading PDB/CIF files and fetching chem comp metadata.

Notes:
- PDB IDs are case-insensitive; we standardize to lower-case for URLs.
- These utilities are used by pipeline scripts.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import requests


@dataclass(frozen=True)
class DownloadResult:
    pdb_id: str
    fmt: str
    url: str
    path: Path
    ok: bool
    status_code: Optional[int] = None
    error: Optional[str] = None


def download_structure(
    pdb_id: str,
    fmt: str,
    out_dir: str | Path,
    *,
    overwrite: bool = False,
    timeout: int = 60,
    sleep_s: float = 0.0,
) -> DownloadResult:
    """Download a structure file from RCSB.

    Args:
        pdb_id: 4-character PDB ID.
        fmt: 'pdb' or 'cif'.
        out_dir: Output directory.
        overwrite: If False and file exists, returns ok=True without re-downloading.
        timeout: HTTP timeout seconds.
        sleep_s: Optional sleep between requests.

    Returns:
        DownloadResult with status.
    """
    fmt = fmt.lower().lstrip('.')
    if fmt not in {"pdb", "cif"}:
        raise ValueError("fmt must be 'pdb' or 'cif'")

    pdb_id_clean = pdb_id.strip().lower()
    url = f"https://files.rcsb.org/download/{pdb_id_clean}.{fmt}"

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pdb_id_clean}.{fmt}"

    if out_path.exists() and not overwrite:
        return DownloadResult(pdb_id=pdb_id_clean, fmt=fmt, url=url, path=out_path, ok=True, status_code=None, error=None)

    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            return DownloadResult(pdb_id=pdb_id_clean, fmt=fmt, url=url, path=out_path, ok=False, status_code=r.status_code, error=f"HTTP {r.status_code}")

        out_path.write_text(r.text)
        if sleep_s:
            time.sleep(sleep_s)
        return DownloadResult(pdb_id=pdb_id_clean, fmt=fmt, url=url, path=out_path, ok=True, status_code=200, error=None)
    except Exception as e:
        return DownloadResult(pdb_id=pdb_id_clean, fmt=fmt, url=url, path=out_path, ok=False, status_code=None, error=str(e))


def fetch_chemcomp_smiles(comp_id: str, *, timeout: int = 30) -> Optional[str]:
    """Fetch a ligand (chem comp) SMILES string from RCSB, if available."""
    cid = comp_id.strip().lower()
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{cid}"
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            return None
        data = r.json()
        descriptors = data.get("pdbx_chem_comp_descriptor", [])
        for d in descriptors:
            if str(d.get("type", "")).upper() == "SMILES":
                return d.get("descriptor")
        return None
    except Exception:
        return None
