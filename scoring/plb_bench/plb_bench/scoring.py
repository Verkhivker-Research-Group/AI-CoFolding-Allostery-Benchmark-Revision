"""Phase 4: OpenStructure execution.

We compute three metrics per (pdb_id, producer, model_idx):

* **BiSyRMSD** — symmetry-corrected ligand RMSD (``LigandScorer.rmsd_details``)
* **lDDT-PLI** — Local Distance Difference Test on the protein-ligand
  interface (``LigandScorer.lddt_pli_details``)
* **QS-global** — global quaternary-structure score
  (``ost.mol.alg.qsscore.QSScorer``)

The OST API surface changed multiple times between 2.4 and 2.7. We use the
2.7+ ``scoring.LigandScorer`` entry point and the ``qsscore`` module, with
guards so an API drift in one scorer doesn't sink the whole row.
"""
from __future__ import annotations

import io
import logging
import tempfile
from pathlib import Path
from typing import Any

from .sanitizer import SanitizedCIF
from .references import read_reference_text

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# OST loading helpers
# ---------------------------------------------------------------------------
def _load_mmcif_text(cif_text: str):
    """Load an mmCIF string into an OST EntityHandle.

    OST's io.LoadMMCIF only accepts a filename, so we write to a temp file.
    We prefer ``process=True`` so entity/compound info is filled in, but fall
    back to a stripped load if the "proper" path trips on residual quirks.
    """
    import ost
    from ost import io as ost_io

    with tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as tmp:
        tmp.write(cif_text)
        tmp_path = tmp.name

    try:
        try:
            ent, seqres, info = ost_io.LoadMMCIF(
                tmp_path, seqres=True, info=True, fault_tolerant=True,
            )
            return ent, seqres, info
        except Exception as e:
            log.warning("strict MMCIF load failed (%s); retrying fault-tolerant", e)
            ent = ost_io.LoadMMCIF(tmp_path, fault_tolerant=True)
            return ent, None, None
    finally:
        try:
            Path(tmp_path).unlink()
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Ligand extraction helper
# ---------------------------------------------------------------------------
_SOLVENT = {"HOH", "DOD", "WAT", "EDO", "GOL", "PEG", "PO4", "SO4", "ACT"}

# OST chain types that are polymer or solvent (not ligands).
# Built lazily so we don't import ost at module level.
_POLY_CHAIN_TYPES: set | None = None

def _polymer_chain_types() -> set:
    global _POLY_CHAIN_TYPES
    if _POLY_CHAIN_TYPES is not None:
        return _POLY_CHAIN_TYPES
    from ost import mol
    _POLY_CHAIN_TYPES = set()
    for name in (
        "CHAINTYPE_POLY_PEPTIDE_L", "CHAINTYPE_POLY_PEPTIDE_D",
        "CHAINTYPE_POLY_DN",        "CHAINTYPE_POLY_RN",
        "CHAINTYPE_POLY",           "CHAINTYPE_POLY_SAC_D",
        "CHAINTYPE_POLY_SAC_L",     "CHAINTYPE_WATER",
    ):
        ct = getattr(mol, name, None)
        if ct is not None:
            _POLY_CHAIN_TYPES.add(ct)
    return _POLY_CHAIN_TYPES


def _get_ligands(ent):
    """Return non-solvent ligand ResidueHandles using OST chain_type."""
    poly = _polymer_chain_types()
    ligs = []
    for chain in ent.chains:
        if chain.chain_type not in poly:
            for res in chain.residues:
                if res.name not in _SOLVENT:
                    ligs.append(res)
    return ligs


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------
def _compute_ligand_metrics(model_ent, ref_ent) -> dict[str, float | None]:
    """Run OST's LigandScorer, return {bisy_rmsd, lddt_pli}."""
    from ost.mol.alg import ligand_scoring_scrmsd, ligand_scoring_lddtpli

    out: dict[str, float | None] = {"bisy_rmsd": None, "lddt_pli": None}

    model_ligs  = _get_ligands(model_ent)
    target_ligs = _get_ligands(ref_ent)

    if not model_ligs or not target_ligs:
        log.warning("No ligands found in model or reference — skipping ligand metrics")
        return out

    def _extract_float(v, *keys):
        """Pull a float out of a score value — handles dicts, objects, and plain numbers."""
        if v is None:
            return None
        if isinstance(v, dict):
            # Try named keys first
            for k in keys:
                try:
                    inner = v[k]
                    if inner is not None:
                        return float(inner)
                except (KeyError, TypeError, ValueError):
                    continue
            # Fallback: convert first non-None value to float
            for inner in v.values():
                if inner is not None:
                    try:
                        return float(inner)
                    except (TypeError, ValueError):
                        continue
            return None
        # Plain number or object with __float__
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    # BiSyRMSD via symmetry-corrected RMSD scorer
    try:
        sc = ligand_scoring_scrmsd.SCRMSDScorer(
            model=model_ent, target=ref_ent,
            model_ligands=model_ligs, target_ligands=target_ligs,
            substructure_match=True,
        )
        scores = sc.score
        if scores:
            vals = [_extract_float(v, "rmsd", "score", "bisy_rmsd")
                    for v in scores.values()]
            vals = [v for v in vals if v is not None]
            if vals:
                out["bisy_rmsd"] = float(min(vals))
    except Exception as e:
        log.warning("BiSyRMSD failed: %s", e)

    # lDDT-PLI
    try:
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(
            model=model_ent, target=ref_ent,
            model_ligands=model_ligs, target_ligands=target_ligs,
            substructure_match=True,
        )
        scores = sc.score
        if scores:
            vals = [_extract_float(v, "lddt_pli", "lddt", "score")
                    for v in scores.values()]
            vals = [v for v in vals if v is not None]
            if vals:
                out["lddt_pli"] = float(max(vals))
    except Exception as e:
        log.warning("lDDT-PLI failed: %s", e)

    return out


def _compute_qs_global(model_ent, ref_ent) -> float | None:
    """Compute QS-global using the original QSEntity/QSScorer API.

    The original pipeline (4_score/run_all_metrics.py) used QSEntity +
    QSScorer directly, which includes ALL chains — both protein AND ligand.
    This means QS-global measures protein-ligand interface contact similarity,
    producing continuous values (0.8–1.0 range) that reflect how well the
    model reproduced the binding contacts.

    Newer OST API variants (ChainMapper-based, OST 2.7+) only map polymer
    chains and strip the ligand, reducing QS-global to a trivial 1.0 for all
    monomer predictions (since there are no protein-protein interfaces in
    single-chain predictions). We therefore try the original QSEntity API
    first, and fall back to the newer variants only if it is unavailable.
    """
    # 0. Original API: QSEntity + QSScorer directly — includes ligand chain.
    #    This matches 4_score/run_all_metrics.py and produces continuous scores.
    try:
        from ost.mol.alg.qsscore import QSEntity, QSScorer as _QSScorer
        mdl_q = QSEntity(model_ent)
        ref_q = QSEntity(ref_ent)
        result = _QSScorer(ref_q, mdl_q).Score()
        v = getattr(result, "qs_global", None)
        if v is not None:
            return float(v)
    except Exception:
        pass

    # 1. High-level ost.mol.alg.scoring.Scorer (OST 2.9+)
    try:
        from ost.mol.alg.scoring import Scorer
        sc = Scorer(model_ent, ref_ent)
        v = sc.qs_global
        if v is not None:
            return float(v)
    except Exception:
        pass

    # 2. QSScorer via ChainMapper (OST 2.7 / 2.8) — polymer chains only,
    #    returns 1.0 for monomer predictions (ligand excluded).
    try:
        from ost.mol.alg.qsscore import QSScorer
        from ost.mol.alg import chain_mapping as _cm
        mapper = _cm.ChainMapper(ref_ent, n_max_naive=1)
        mapping = mapper.GetMapping(model_ent)
        sc = QSScorer(model_ent, mapping.alns, ref_ent)
        return float(sc.global_score)
    except Exception:
        pass

    # 3. QSScorer.FromEntities class method
    try:
        from ost.mol.alg.qsscore import QSScorer
        if hasattr(QSScorer, "FromEntities"):
            return float(QSScorer.FromEntities(model_ent, ref_ent).global_score)
    except Exception:
        pass

    # 4. QSScorer.FromMappingResult class method
    try:
        from ost.mol.alg.qsscore import QSScorer
        if hasattr(QSScorer, "FromMappingResult"):
            mapping = QSScorer.GetMappingResult(model_ent, ref_ent)
            return float(QSScorer.FromMappingResult(mapping).global_score)
    except Exception:
        pass

    log.warning("QS-global: all API attempts exhausted")
    return None


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def score_model(
    sanitized_model: SanitizedCIF,
    reference_cif_path: Path,
) -> dict[str, Any]:
    """Score one sanitized model against its reference mmCIF.

    Returns a dict with keys: bisy_rmsd, lddt_pli, qs_global.
    Missing / failed metrics are ``None``. This function *never* raises —
    callers rely on a best-effort contract.
    """
    result: dict[str, Any] = {
        "bisy_rmsd": None, "lddt_pli": None, "qs_global": None,
    }
    try:
        ref_text = read_reference_text(reference_cif_path)
    except Exception as e:
        log.error("cannot read reference %s: %s", reference_cif_path, e)
        return result

    try:
        model_ent, _, _ = _load_mmcif_text(sanitized_model.text)
    except Exception as e:
        log.error("cannot load sanitized model into OST: %s", e)
        return result

    try:
        ref_ent, _, _ = _load_mmcif_text(ref_text)
    except Exception as e:
        log.error("cannot load reference into OST: %s", e)
        return result

    # Ligand metrics
    lig = _compute_ligand_metrics(model_ent, ref_ent)
    result["bisy_rmsd"] = lig["bisy_rmsd"]
    result["lddt_pli"] = lig["lddt_pli"]

    # QS-global
    result["qs_global"] = _compute_qs_global(model_ent, ref_ent)

    return result
