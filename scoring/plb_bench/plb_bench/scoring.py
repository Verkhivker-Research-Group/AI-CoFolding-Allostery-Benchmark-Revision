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
import re
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
    """Return non-solvent ligand ResidueHandles using OST chain_type.

    Falls back to a residue-name scan for ``LIG`` when the chain_type
    approach returns nothing — this handles DynamicBind CIFs where the
    ligand chain may be assigned an unexpected chain type at load time.
    """
    poly = _polymer_chain_types()
    ligs = []
    for chain in ent.chains:
        if chain.chain_type not in poly:
            for res in chain.residues:
                if res.name not in _SOLVENT:
                    ligs.append(res)
    # Fallback: scan ALL chains for LIG residues (DynamicBind chain_type fix)
    if not ligs:
        for chain in ent.chains:
            for res in chain.residues:
                if res.name == "LIG" and res.atom_count > 1:
                    ligs.append(res)
    return ligs


# ---------------------------------------------------------------------------
# CIF bond-wiring helpers (DynamicBind LIG fix)
# ---------------------------------------------------------------------------
_CIF_BOND_ORDERS = {"SING": 1, "DOUB": 2, "TRIP": 3, "AROM": 1}

def _parse_cif_lig_bonds(cif_text: str) -> list[tuple[str, str, int]]:
    """Parse ``_chem_comp_bond`` records for the ``LIG`` residue out of a CIF string.

    Returns a list of (atom1_name, atom2_name, bond_order_int) triples.
    The sanitizer writes these rows as::

        LIG  C1  N2  SING
        LIG  C1  C3  DOUB
        ...

    We match exactly that format (comp_id must be ``LIG``).
    """
    pattern = re.compile(
        r'^LIG\s+(\S+)\s+(\S+)\s+(SING|DOUB|TRIP|AROM)',
        re.MULTILINE,
    )
    bonds = []
    for m in pattern.finditer(cif_text):
        a1, a2, order_str = m.group(1), m.group(2), m.group(3)
        bonds.append((a1, a2, _CIF_BOND_ORDERS.get(order_str, 1)))
    return bonds


def _wire_lig_bonds_into_entity(model_ent, cif_text: str) -> None:
    """Wire ``_chem_comp_bond`` bonds from *cif_text* into ``model_ent`` in place.

    OST's ``LoadMMCIF`` reads ``_chem_comp_bond`` into ``MMCifInfo`` but does
    **not** create ``BondHandle`` objects for unknown residue names like ``LIG``.
    This function parses those records and connects atoms directly in the
    already-loaded entity, keeping ligand handles in the *same* entity object
    (avoiding the cross-entity bug that silently breaks ``SCRMSDScorer``).

    Uses ``STANDARD_EDIT`` so bond handles are immediately visible to scorers.
    No-op when there are no ``_chem_comp_bond`` rows for ``LIG``.
    """
    bonds = _parse_cif_lig_bonds(cif_text)
    if not bonds:
        return

    # Collect all LIG residues in the entity (usually just one)
    lig_residues = [
        res for ch in model_ent.chains
        for res in ch.residues
        if res.name == "LIG"
    ]
    if not lig_residues:
        return

    try:
        import ost.mol as _omol
        edi = model_ent.EditXCS(_omol.STANDARD_EDIT)
        wired = 0
        for lig_res in lig_residues:
            atom_map = {a.name: a for a in lig_res.atoms}
            for a1_name, a2_name, bo in bonds:
                a1 = atom_map.get(a1_name)
                a2 = atom_map.get(a2_name)
                if a1 is None or a2 is None:
                    continue
                try:
                    edi.Connect(a1, a2, bo)
                    wired += 1
                except Exception:
                    pass
        # Commit the editor (destructor also commits, but be explicit)
        del edi
        log.debug("_wire_lig_bonds_into_entity: wired %d bonds", wired)
    except Exception as e:
        log.warning("_wire_lig_bonds_into_entity failed: %s", e)


_IONS_AND_COFACTORS = _SOLVENT | {
    # Common mono-atomic ions and tiny cofactors that should not be treated as
    # the primary ligand when scanning a reference structure.
    "MG", "ZN", "CA", "NA", "CL", "K",  "MN", "FE", "CU", "CO", "NI",
    "CD", "HG", "PB", "AL", "CS", "RB", "IN", "TL", "SR", "BA", "LI",
    "BR", "F",  "I",  "NH4", "NO3", "IOD", "FLC",
}


def _get_ref_lig_name(ref_ent) -> str | None:
    """Return the compound name of the primary ligand in *ref_ent*.

    Skips water, common ions, and single-atom residues so that metals or
    buffer molecules do not masquerade as the ligand of interest.
    Used to rename the reference ligand → ``LIG`` for QS-global scoring so
    OST's chain mapper can pair model (always ``LIG``) and reference chains.
    """
    poly = _polymer_chain_types()
    for chain in ref_ent.chains:
        if chain.chain_type not in poly:
            for res in chain.residues:
                if res.name not in _IONS_AND_COFACTORS and res.atom_count > 1:
                    return res.name
    return None


def _load_sdf_as_ligand(sdf_path: str) -> list:
    """Load an SDF file into a list of ResidueHandles with explicit bond connectivity.

    Used for DynamicBind (and any producer that stores the ligand as a separate
    SDF rather than embedding it in the model CIF).  OST's LoadMMCIF cannot
    assign bond connectivity to unknown residue names like ``LIG``; loading the
    SDF directly via RDKit and building an OST entity with explicit bonds gives
    ``SCRMSDScorer`` the molecular graph it needs for substructure matching.

    Returns an empty list on any failure so callers fall back to the CIF path.
    """
    # ── Attempt 1: OST native SDF reader (available in some OST builds) ──────
    try:
        from ost import io as _ost_io
        if hasattr(_ost_io, "LoadSDF"):
            ent = _ost_io.LoadSDF(sdf_path)
            ligs = [res for ch in ent.chains for res in ch.residues
                    if res.name not in _SOLVENT and res.atom_count > 0]
            if ligs:
                log.debug("_load_sdf_as_ligand: OST LoadSDF ok, %d lig(s)", len(ligs))
                return ligs
    except Exception as e:
        log.debug("OST LoadSDF failed (%s); trying RDKit", e)

    # ── Attempt 2: RDKit → OST entity builder ────────────────────────────────
    try:
        from rdkit import Chem
        import ost.mol as _omol

        suppl = Chem.SDMolSupplier(sdf_path, removeHs=True, sanitize=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            suppl2 = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
            mol = next((m for m in suppl2 if m is not None), None)
        if mol is None or not mol.GetNumConformers():
            log.warning("_load_sdf_as_ligand: RDKit returned None for %s", sdf_path)
            return []

        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass  # continue with unsanitized mol

        conf = mol.GetConformer(0)
        ent = _omol.CreateEntity()
        edi = ent.EditXCS(_omol.BUFFERED_EDIT)
        ch  = edi.InsertChain("L")
        res = edi.AppendResidue(ch, "LIG")

        atom_handles: dict[int, object] = {}
        for i, atom in enumerate(mol.GetAtoms()):
            pos  = conf.GetAtomPosition(i)
            name = f"{atom.GetSymbol()}{i + 1}"[:4]
            ah   = edi.InsertAtom(res, name,
                                  _omol.Vec3(pos.x, pos.y, pos.z))
            try:
                ah.element = atom.GetSymbol()
            except Exception:
                pass
            atom_handles[i] = ah

        _BO = {
            Chem.rdchem.BondType.SINGLE:   1,
            Chem.rdchem.BondType.DOUBLE:   2,
            Chem.rdchem.BondType.TRIPLE:   3,
            Chem.rdchem.BondType.AROMATIC: 1,
        }
        for bond in mol.GetBonds():
            bo = _BO.get(bond.GetBondType(), 1)
            try:
                edi.Connect(atom_handles[bond.GetBeginAtomIdx()],
                            atom_handles[bond.GetEndAtomIdx()], bo)
            except Exception:
                pass

        ligs = [r for c in ent.chains for r in c.residues]
        log.debug("_load_sdf_as_ligand: RDKit builder ok, %d atom(s) in lig",
                  ligs[0].atom_count if ligs else 0)
        return ligs

    except Exception as e:
        log.warning("_load_sdf_as_ligand failed for %s: %s", sdf_path, e)
        return []


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------
def _compute_ligand_metrics(model_ent, ref_ent,
                            sdf_path: str | None = None) -> dict[str, float | None]:
    """Run OST's LigandScorer, return {bisy_rmsd, lddt_pli}.

    Model-ligand strategy (two candidates tried in order):
      A. RDKit mol loaded directly from the SDF file (OST 2.9+ accepts
         ``rdkit.Chem.Mol`` as ``model_ligands``). Gives correct bond graph,
         correct heavy-atom count (removeHs=True), and avoids the cross-entity
         issue entirely. This is the primary path for DynamicBind.
      B. ResidueHandle extracted from model_ent via ``_get_ligands`` (with
         bonds wired by ``_wire_lig_bonds_into_entity``). Used as fallback
         when ``sdf_path`` is absent or the RDKit approach fails.

    For each candidate we try ``substructure_match=True`` first, then ``False``
    as a fallback.
    """
    from ost.mol.alg import ligand_scoring_scrmsd, ligand_scoring_lddtpli

    out: dict[str, float | None] = {"bisy_rmsd": None, "lddt_pli": None}

    # ── Build ordered list of model-ligand candidates ─────────────────────────
    mdl_lig_candidates: list = []

    # Candidate A: RDKit mol from SDF (OST 2.9+ natively accepts Chem.Mol)
    if sdf_path:
        _sdf_p = Path(sdf_path)
        if _sdf_p.exists():
            try:
                from rdkit import Chem as _Chem
                _rdmol = None
                for _rm, _san in ((True, True), (True, False), (False, False)):
                    _suppl = _Chem.SDMolSupplier(str(_sdf_p), removeHs=_rm, sanitize=_san)
                    _rdmol = next((_m for _m in _suppl
                                   if _m is not None and _m.GetNumConformers() > 0), None)
                    if _rdmol is not None:
                        if not _san:
                            try:
                                _Chem.SanitizeMol(_rdmol)
                            except Exception:
                                pass
                        break
                if _rdmol is not None:
                    mdl_lig_candidates.append([_rdmol])
                    log.info("ligand: RDKit mol from SDF  (%d heavy atoms)",
                             _rdmol.GetNumHeavyAtoms())
                else:
                    log.warning("ligand: RDKit could not load SDF %s", sdf_path)
            except Exception as _e:
                log.warning("ligand: SDF→RDKit failed (%s)", _e)
        else:
            log.warning("ligand: sdf_path not found on disk: %s", sdf_path)

    # Candidate B: entity handles (bonds pre-wired)
    _ent_ligs = _get_ligands(model_ent)
    if _ent_ligs:
        mdl_lig_candidates.append(_ent_ligs)
        log.info("ligand: entity handles  (%d residue(s) found)", len(_ent_ligs))
    else:
        log.warning("ligand: _get_ligands returned nothing for model entity")

    target_ligs = _get_ligands(ref_ent)
    if not target_ligs:
        log.warning("ligand: no target ligands in reference — skipping metrics")
        return out

    if not mdl_lig_candidates:
        log.warning("ligand: no model ligand candidates — skipping metrics")
        return out

    # ── Shared score extractor ────────────────────────────────────────────────
    def _extract_float(v, *keys):
        if v is None:
            return None
        if isinstance(v, dict):
            for k in keys:
                try:
                    inner = v[k]
                    if inner is not None:
                        return float(inner)
                except (KeyError, TypeError, ValueError):
                    continue
            for inner in v.values():
                if inner is not None:
                    try:
                        return float(inner)
                    except (TypeError, ValueError):
                        continue
            return None
        # OST result objects: try attribute access before float()
        for attr in ("rmsd", "bisy_rmsd", "lddt_pli", "lddt", "score", "value"):
            try:
                av = getattr(v, attr, None)
                if av is not None:
                    return float(av)
            except (TypeError, ValueError):
                continue
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    # ── BiSyRMSD ─────────────────────────────────────────────────────────────
    for _cand in mdl_lig_candidates:
        _cand_label = type(_cand[0]).__name__
        for _sm in (True, False):
            try:
                sc = ligand_scoring_scrmsd.SCRMSDScorer(
                    model=model_ent, target=ref_ent,
                    model_ligands=_cand, target_ligands=target_ligs,
                    substructure_match=_sm,
                )
                scores = sc.score
                if scores:
                    vals = [_extract_float(v, "rmsd", "score", "bisy_rmsd")
                            for v in scores.values()]
                    vals = [v for v in vals if v is not None]
                    if vals:
                        out["bisy_rmsd"] = float(min(vals))
                        log.info("bisy_rmsd=%.3f  (cand=%s, sm=%s)",
                                 out["bisy_rmsd"], _cand_label, _sm)
                        break
                    else:
                        log.info("BiSyRMSD: scorer ran but scores dict empty "
                                 "(cand=%s, sm=%s)", _cand_label, _sm)
                else:
                    log.info("BiSyRMSD: sc.score is empty/None "
                             "(cand=%s, sm=%s)", _cand_label, _sm)
            except Exception as _e:
                log.info("BiSyRMSD exception (cand=%s, sm=%s): %s",
                         _cand_label, _sm, _e)
        if out["bisy_rmsd"] is not None:
            break

    # ── lDDT-PLI ─────────────────────────────────────────────────────────────
    for _cand in mdl_lig_candidates:
        _cand_label = type(_cand[0]).__name__
        for _sm in (True, False):
            try:
                sc = ligand_scoring_lddtpli.LDDTPLIScorer(
                    model=model_ent, target=ref_ent,
                    model_ligands=_cand, target_ligands=target_ligs,
                    substructure_match=_sm,
                )
                scores = sc.score
                if scores:
                    vals = [_extract_float(v, "lddt_pli", "lddt", "score")
                            for v in scores.values()]
                    vals = [v for v in vals if v is not None]
                    if vals:
                        out["lddt_pli"] = float(max(vals))
                        log.info("lddt_pli=%.3f  (cand=%s, sm=%s)",
                                 out["lddt_pli"], _cand_label, _sm)
                        break
                    else:
                        log.info("lDDT-PLI: scores dict empty "
                                 "(cand=%s, sm=%s)", _cand_label, _sm)
                else:
                    log.info("lDDT-PLI: sc.score is empty/None "
                             "(cand=%s, sm=%s)", _cand_label, _sm)
            except Exception as _e:
                log.info("lDDT-PLI exception (cand=%s, sm=%s): %s",
                         _cand_label, _sm, _e)
        if out["lddt_pli"] is not None:
            break

    return out


def _compute_qs_global(model_ent, ref_ent,
                       ref_cif_text: str | None = None) -> float | None:
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

    DynamicBind QS fix
    ------------------
    DynamicBind models use ``LIG`` as the ligand residue name.  The reference
    uses the real CCD code (e.g. ``DTP``).  OST's QSScorer cannot pair these
    chains because the names differ → binary 0/1 scores.

    Previous approach (rename model LIG → DTP) failed: OST recognises the
    real CCD name, applies its compound library, and conflicts with the
    SDF-derived atom names → the re-loaded entity is corrupted or the load
    throws, so we silently fell back to the original (unmapped) entity.

    Correct approach: rename the REFERENCE ligand (``DTP`` → ``LIG``) instead.
    ``LIG`` is unknown to OST's compound library so the renamed reference CIF
    loads cleanly as a generic non-polymer residue.  Now both model and
    reference have ``LIG`` → QSScorer pairs the ligand chains correctly →
    protein-ligand interface is evaluated → continuous QS scores.
    """
    # ── DynamicBind fix: rename ref compound → LIG so chains can be paired ──
    qs_ref_ent = ref_ent  # default: use reference as-is
    model_has_lig = any(
        res.name == "LIG"
        for ch in model_ent.chains
        for res in ch.residues
    )
    if model_has_lig and ref_cif_text is not None:
        ref_lig_name = _get_ref_lig_name(ref_ent)
        if ref_lig_name and ref_lig_name != "LIG":
            try:
                renamed_ref = re.sub(
                    r'\b' + re.escape(ref_lig_name) + r'\b',
                    "LIG",
                    ref_cif_text,
                )
                qs_ref_ent, _, _ = _load_mmcif_text(renamed_ref)
                log.info("QS: renamed ref %s→LIG for chain pairing", ref_lig_name)
            except Exception as e:
                log.warning("QS ref rename %s→LIG failed (%s); using original ref",
                            ref_lig_name, e)
                qs_ref_ent = ref_ent

    # 0. Original API: QSEntity + QSScorer directly — includes ligand chain.
    #    This matches 4_score/run_all_metrics.py and produces continuous scores.
    try:
        from ost.mol.alg.qsscore import QSEntity, QSScorer as _QSScorer
        mdl_q = QSEntity(model_ent)
        ref_q = QSEntity(qs_ref_ent)
        result = _QSScorer(ref_q, mdl_q).Score()
        v = getattr(result, "qs_global", None)
        if v is not None:
            return float(v)
    except Exception:
        pass

    # 1. High-level ost.mol.alg.scoring.Scorer (OST 2.9+)
    try:
        from ost.mol.alg.scoring import Scorer
        sc = Scorer(model_ent, qs_ref_ent)
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
        mapper = _cm.ChainMapper(qs_ref_ent, n_max_naive=1)
        mapping = mapper.GetMapping(model_ent)
        sc = QSScorer(model_ent, mapping.alns, qs_ref_ent)
        return float(sc.global_score)
    except Exception:
        pass

    # 3. QSScorer.FromEntities class method
    try:
        from ost.mol.alg.qsscore import QSScorer
        if hasattr(QSScorer, "FromEntities"):
            return float(QSScorer.FromEntities(model_ent, qs_ref_ent).global_score)
    except Exception:
        pass

    # 4. QSScorer.FromMappingResult class method
    try:
        from ost.mol.alg.qsscore import QSScorer
        if hasattr(QSScorer, "FromMappingResult"):
            mapping = QSScorer.GetMappingResult(model_ent, qs_ref_ent)
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
    ligand_sdf_path: str | None = None,
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

    # Wire _chem_comp_bond bonds from the sanitized CIF text into model_ent.
    # OST ignores these bonds for unknown residue names (e.g. LIG) during load,
    # so we patch them in manually using atom-name matching on STANDARD_EDIT.
    # No-op for producers that don't embed _chem_comp_bond records.
    _wire_lig_bonds_into_entity(model_ent, sanitized_model.text)

    # Ligand metrics
    lig = _compute_ligand_metrics(model_ent, ref_ent, sdf_path=ligand_sdf_path)
    result["bisy_rmsd"] = lig["bisy_rmsd"]
    result["lddt_pli"] = lig["lddt_pli"]

    # QS-global — pass ref CIF text so the QS path can rename the reference
    # ligand (e.g. DTP) → LIG, enabling correct chain pairing with DynamicBind
    # models that always use LIG as the residue name.
    result["qs_global"] = _compute_qs_global(
        model_ent, ref_ent, ref_cif_text=ref_text)

    return result
