"""Phase 3: Automated reference retrieval.

For each ``pdb_id`` we need the experimental reference structure.
Lookup order:

1. ``<refs_dir>/<PDB_ID>.cif.gz`` (biological assembly 1)
2. ``<refs_dir>/<PDB_ID>.cif``
3. ``<refs_dir>/<pdb_id>.cif[.gz]`` (lowercase)
4. Download from RCSB: ``https://files.rcsb.org/download/<PDB_ID>-assembly1.cif.gz``

Downloads are cached into ``refs_dir`` so concurrent workers don't re-fetch.
"""
from __future__ import annotations

import gzip
import logging
import os
import tempfile
import threading
from pathlib import Path

import requests

log = logging.getLogger(__name__)

RCSB_ASSEMBLY_URL = "https://files.rcsb.org/download/{pdb}-assembly1.cif.gz"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb}.cif.gz"

# Per-pdb file locks to prevent two workers downloading the same reference.
_download_locks: dict[str, threading.Lock] = {}
_locks_lock = threading.Lock()


def _lock_for(pdb_id: str) -> threading.Lock:
    with _locks_lock:
        lk = _download_locks.get(pdb_id)
        if lk is None:
            lk = _download_locks[pdb_id] = threading.Lock()
        return lk


def _search_local(refs_dir: Path, pdb_id: str) -> Path | None:
    for variant in (pdb_id.upper(), pdb_id.lower()):
        for name in (f"{variant}.cif", f"{variant}.cif.gz",
                     f"{variant}-assembly1.cif", f"{variant}-assembly1.cif.gz"):
            p = refs_dir / name
            if p.exists():
                return p
    return None


def _atomic_write(dest: Path, data: bytes) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Atomic: write to a sibling temp file, then rename. Prevents partial
    # files being read by concurrent workers.
    fd, tmp_name = tempfile.mkstemp(dir=str(dest.parent), prefix=".tmp_",
                                     suffix=dest.suffix)
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        os.replace(tmp_name, dest)
    except Exception:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise


def _download(pdb_id: str, refs_dir: Path, timeout: float = 30.0) -> Path:
    pdb_up = pdb_id.upper()
    lk = _lock_for(pdb_up)
    with lk:
        # Re-check after acquiring the lock in case another thread won the race.
        existing = _search_local(refs_dir, pdb_id)
        if existing is not None:
            return existing

        last_err: Exception | None = None
        for url in (RCSB_ASSEMBLY_URL.format(pdb=pdb_up),
                    RCSB_CIF_URL.format(pdb=pdb_up)):
            try:
                resp = requests.get(url, timeout=timeout)
                if resp.status_code == 200:
                    dest = refs_dir / url.rsplit("/", 1)[-1]
                    _atomic_write(dest, resp.content)
                    log.info("cached reference %s -> %s", pdb_up, dest.name)
                    return dest
                last_err = RuntimeError(f"HTTP {resp.status_code} for {url}")
            except Exception as e:
                last_err = e
                log.warning("download failed %s: %s", url, e)
        raise RuntimeError(f"no reference available for {pdb_up}: {last_err}")


def get_reference(pdb_id: str, refs_dir: Path,
                  allow_download: bool = True) -> tuple[Path, str]:
    """Return (path, source) where source is ``"local"`` or ``"rcsb"``."""
    refs_dir = Path(refs_dir)
    refs_dir.mkdir(parents=True, exist_ok=True)

    local = _search_local(refs_dir, pdb_id)
    if local is not None:
        return local, "local"
    if not allow_download:
        raise FileNotFoundError(
            f"reference for {pdb_id} not in {refs_dir} and downloads disabled")
    return _download(pdb_id, refs_dir), "rcsb"


def read_reference_text(path: Path) -> str:
    """Return the mmCIF text content, transparently decompressing .gz."""
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            return fh.read()
    return path.read_text()
