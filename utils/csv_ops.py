"""Small CSV dataframe helpers used across aggregation scripts."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd


def read_csv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path)


def safe_mkdir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def write_csv(df: pd.DataFrame, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def concat_csvs(paths: Iterable[str | Path]) -> pd.DataFrame:
    dfs: List[pd.DataFrame] = []
    for p in paths:
        dfs.append(read_csv(p))
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def best_row_per_group(df: pd.DataFrame, group_col: str, score_col: str, keep: str = "min") -> pd.DataFrame:
    """Return one row per group based on min/max of score_col."""
    if df.empty:
        return df
    if keep not in {"min", "max"}:
        raise ValueError("keep must be 'min' or 'max'")
    df2 = df.copy()
    df2[score_col] = pd.to_numeric(df2[score_col], errors="coerce")
    if keep == "min":
        idx = df2.groupby(group_col)[score_col].idxmin()
    else:
        idx = df2.groupby(group_col)[score_col].idxmax()
    return df2.loc[idx].copy()
