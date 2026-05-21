"""ID / label normalization utilities.

These functions centralize the ad-hoc ID rules that were scattered throughout conc.ipynb.
They are intentionally conservative: they only apply transformations that appeared in the
notebook (remove underscores, case-normalize, handle '*seed*' style strings, etc.).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional, Tuple


_SEED_RE = re.compile(r"seed", re.IGNORECASE)


def normalize_pdb_id(raw: Optional[str]) -> Optional[str]:
    """Normalize a PDB/complex identifier to a stable join key.

    Notebook behaviors replicated:
    - if the string contains 'seed' => take first 5 characters
    - else remove underscores and lowercase
    """
    if raw is None:
        return None
    s = str(raw)
    if s.strip() == "" or s.lower() == "nan":
        return None
    if _SEED_RE.search(s):
        return s[:5].lower()
    return s.replace("_", "").lower()


def normalize_join_key_upper_no_underscore(raw: Optional[str]) -> Optional[str]:
    """Uppercase join key used in several notebook merges."""
    if raw is None:
        return None
    s = str(raw)
    if s.strip() == "" or s.lower() == "nan":
        return None
    return s.replace("_", "").upper()


def extract_prefix(raw: Optional[str], n: int = 5) -> Optional[str]:
    """Extract the first N characters (used heavily for 5-char model codes)."""
    if raw is None:
        return None
    s = str(raw)
    if s.strip() == "" or s.lower() == "nan":
        return None
    return s[:n]


@dataclass(frozen=True)
class SourceDataset:
    source: str
    dataset: str


def parse_source_and_dataset(source_cell: str) -> SourceDataset:
    """Map a filename/source string to (source, dataset).

    Notebook behavior replicated:
    - Source prefix mapping: boltz*, af3*, protenix* -> canonical
    - dataset inferred by substring 'main' or 'allo'
    """
    src_lower = str(source_cell).lower()

    if src_lower.startswith("boltz"):
        source = "boltz"
    elif src_lower.startswith("af3"):
        source = "af3"
    elif src_lower.startswith("protenix"):
        source = "protenix"
    elif src_lower.startswith("chai"):
        # notebook used both 'chai' and 'chai1' inconsistently
        source = "chai"
    else:
        source = src_lower

    if "main" in src_lower:
        dataset = "main"
    elif "allo" in src_lower or "allosteric" in src_lower:
        dataset = "allo"
    else:
        dataset = ""

    return SourceDataset(source=source, dataset=dataset)
