"""plb_bench — Protein-Ligand Benchmarking pipeline for AF3/Proteinx/Chai/DynamicBind via OpenStructure."""
from .schema import ModelRecord, ScoreRecord
from .registry import PARSERS, register_parser
from .normalize import (
    DISCOVERERS, register_discoverer,
    normalize_producer, normalize_all, normalize_one,
    NormalizedEntry,
)
from .score_tree import ScoreConfig, run as score_normalized_tree

__version__ = "0.2.0"
__all__ = [
    "ModelRecord", "ScoreRecord",
    "PARSERS", "register_parser",
    "DISCOVERERS", "register_discoverer",
    "normalize_producer", "normalize_all", "normalize_one",
    "NormalizedEntry",
    "ScoreConfig", "score_normalized_tree",
]
