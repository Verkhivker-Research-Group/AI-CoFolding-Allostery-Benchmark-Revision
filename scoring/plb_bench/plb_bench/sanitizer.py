"""Phase 2: Robust CIF Pre-processing Layer.

OpenStructure's mmCIF reader is strict: it needs coherent ``_entity``,
``_entity_poly``, ``_struct_asym`` and ``_chem_comp`` loops. Generative
models routinely skip these. Instead of patching every model's output,
we rebuild the structure in-memory via Biopython + gemmi and emit a
minimal-but-complete mmCIF that OST accepts.

Design choices:

* Everything runs via ``io.StringIO`` / ``tempfile.SpooledTemporaryFile`` —
  no disk round-trips on the hot path.
* The sanitizer is format-agnostic: PDB, mmCIF, and separate ligand SDFs
  all converge to a single clean mmCIF text blob.
* We *never* try to preserve every annotation. We preserve coordinates,
  residue identity, chain IDs, elements, and the protein/ligand split —
  which is all the scoring needs.
"""
from __future__ import annotations

import io
import logging
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Lazy imports — we want the sanitizer importable even on machines where
# only a subset of (biopython, gemmi, rdkit) is installed.
# ---------------------------------------------------------------------------
def _import_biopython():
    from Bio.PDB import PDBParser, MMCIFParser   # type: ignore
    from Bio.PDB.mmcifio import MMCIFIO           # type: ignore
    return PDBParser, MMCIFParser, MMCIFIO


def _import_gemmi():
    import gemmi                                  # type: ignore
    return gemmi


def _import_rdkit():
    from rdkit import Chem                        # type: ignore
    return Chem


# ---------------------------------------------------------------------------
# Canonical output
# ---------------------------------------------------------------------------
@dataclass
class SanitizedCIF:
    """Result of sanitization — kept in-memory."""
    text: str                        # full mmCIF as a string
    pdb_id: str
    n_chains: int
    n_ligand_residues: int
    ligand_chain_ids: tuple[str, ...]

    def write(self, path: Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.text)
        return path

    def as_buffer(self) -> io.StringIO:
        return io.StringIO(self.text)


# ---------------------------------------------------------------------------
# Core rebuild via gemmi (preferred — fastest, writes a proper mmCIF doc)
# ---------------------------------------------------------------------------
def _rebuild_with_gemmi(
    structure_bytes: bytes,
    fmt: str,                         # "pdb" | "cif"
    pdb_id: str,
    extra_ligand_sdf: bytes | None = None,
) -> SanitizedCIF:
    gemmi = _import_gemmi()

    # Read into a gemmi.Structure regardless of input format
    with tempfile.NamedTemporaryFile(suffix=f".{fmt}", delete=True) as tmp:
        tmp.write(structure_bytes)
        tmp.flush()
        if fmt == "pdb":
            st = gemmi.read_pdb(tmp.name)
        else:
            st = gemmi.read_structure(tmp.name)
    st.name = pdb_id.upper()

    # CIFs from generative models often omit fields gemmi's high-level reader
    # needs (pdbx_PDB_model_num, _cell, _entity_poly_seq). When that happens
    # read_structure silently returns 0 models. Recover manually from the raw
    # _atom_site loop — this is the central "broken CIF" case.
    if fmt == "cif" and (len(st) == 0 or _structure_is_empty(st)):
        log.info("gemmi read returned empty structure — rebuilding from raw _atom_site")
        st = _rebuild_from_atom_site(structure_bytes, pdb_id)

    # Append an SDF ligand as its own chain if provided
    if extra_ligand_sdf is not None:
        _append_sdf_as_chain(st, extra_ligand_sdf)

    # Make sure required metadata loops exist. gemmi's setup_entities rebuilds
    # entity / entity_poly / struct_asym from the coordinate records — this is
    # the single most important call in the whole sanitizer.
    st.setup_entities()
    st.assign_label_seq_id()

    # Audit pass: collect ligand chains (non-polymer with HETATM residues)
    ligand_chain_ids: list[str] = []
    n_ligand_res = 0
    for model in st:
        for chain in model:
            has_poly = False
            lig_here = 0
            for res in chain:
                info = gemmi.find_tabulated_residue(res.name)
                if info is not None and info.is_amino_acid():
                    has_poly = True
                elif res.het_flag == "H" or (info is not None and not info.is_water()
                                             and not info.is_nucleic_acid()
                                             and not info.is_amino_acid()):
                    lig_here += 1
            if lig_here and not has_poly:
                ligand_chain_ids.append(chain.name)
                n_ligand_res += lig_here
        break  # first model only

    # Serialize to mmCIF text in-memory. gemmi's Document exposes
    # as_string() which returns the full serialized mmCIF.
    doc = st.make_mmcif_document()
    cif_text = doc.as_string()

    n_chains = len(st[0]) if len(st) else 0
    return SanitizedCIF(
        text=cif_text, pdb_id=pdb_id.upper(),
        n_chains=n_chains, n_ligand_residues=n_ligand_res,
        ligand_chain_ids=tuple(ligand_chain_ids),
    )


def _structure_is_empty(st) -> bool:
    """A gemmi Structure with no atoms is functionally empty."""
    for model in st:
        for chain in model:
            for _res in chain:
                return False
        break
    return True


def _rebuild_from_atom_site(cif_bytes: bytes, pdb_id: str):
    """Rebuild a gemmi.Structure directly from the ``_atom_site`` loop.

    Used when the high-level reader gave up because of missing sibling
    metadata (``pdbx_PDB_model_num``, ``_cell``, etc.). We only need atoms,
    residues, and chains — the downstream ``setup_entities()`` call rebuilds
    the rest.

    Each column is fetched independently via ``find_loop`` so that any single
    missing optional tag does not abort the whole recovery. The mandatory
    minimum is ``Cartn_x/y/z``; everything else has a sensible default.
    """
    gemmi = _import_gemmi()
    doc = gemmi.cif.read_string(cif_bytes.decode("utf-8", errors="replace"))
    block = doc.sole_block()

    def col(tag: str) -> list[str]:
        lp = block.find_loop(f"_atom_site.{tag}")
        return list(lp) if lp else []

    x = col("Cartn_x")
    y = col("Cartn_y")
    z = col("Cartn_z")
    n = len(x)
    if n == 0 or len(y) != n or len(z) != n:
        raise ValueError("no _atom_site loop found")

    def padded(tag: str, default: str) -> list[str]:
        v = col(tag)
        if len(v) == n:
            return v
        return [default] * n

    group_col = padded("group_PDB", "ATOM")
    serial_col = padded("id", "0")
    elem_col = padded("type_symbol", "")
    atom_name_col = padded("label_atom_id", "")
    comp_col = padded("label_comp_id", "UNK")
    label_asym_col = padded("label_asym_id", "A")
    auth_asym_col = padded("auth_asym_id", "")
    label_seq_col = padded("label_seq_id", "")
    auth_seq_col = padded("auth_seq_id", "")
    occ_col = padded("occupancy", "1.00")
    b_col = padded("B_iso_or_equiv", "20.00")
    model_col = padded("pdbx_PDB_model_num", "1")

    def clean(v: str) -> str:
        v = v.strip().strip("'\"")
        if v in ("?", "."):
            return ""
        return v

    st = gemmi.Structure()
    st.name = pdb_id.upper()
    model = None
    last_model_num = None
    chain_by_id: dict[str, "gemmi.Chain"] = {}
    last_res_key: tuple | None = None
    residue = None

    for i in range(n):
        try:
            model_num = int(clean(model_col[i]) or "1")
        except ValueError:
            model_num = 1

        if model is None or model_num != last_model_num:
            _new_model = gemmi.Model(str(model_num))
            st.add_model(_new_model)
            model = st[-1]                # re-bind to the copy gemmi stored
            chain_by_id = {}
            last_model_num = model_num
            last_res_key = None
            residue = None

        group = (clean(group_col[i]) or "ATOM").upper()
        element = clean(elem_col[i])
        atom_name = clean(atom_name_col[i]) or element or "X"
        comp = clean(comp_col[i]) or "UNK"
        asym = clean(label_asym_col[i]) or clean(auth_asym_col[i]) or "A"
        seq_raw = clean(auth_seq_col[i]) or clean(label_seq_col[i]) or "1"
        try:
            seq_num = int(seq_raw)
        except ValueError:
            seq_num = 1

        try:
            cx = float(clean(x[i]) or "0")
            cy = float(clean(y[i]) or "0")
            cz = float(clean(z[i]) or "0")
        except ValueError:
            continue  # skip unparseable rows rather than aborting the whole rebuild

        try:
            occ = float(clean(occ_col[i]) or "1.0")
        except ValueError:
            occ = 1.0
        try:
            b = float(clean(b_col[i]) or "20.0")
        except ValueError:
            b = 20.0

        chain = chain_by_id.get(asym)
        if chain is None:
            _new_chain = gemmi.Chain(asym)
            model.add_chain(_new_chain)
            chain = model[-1]             # re-bind to stored copy
            chain_by_id[asym] = chain
            last_res_key = None

        res_key = (asym, comp, seq_num)
        if res_key != last_res_key:
            _new_res = gemmi.Residue()
            _new_res.name = comp
            _new_res.seqid = gemmi.SeqId(seq_num, " ")
            _new_res.het_flag = "H" if group == "HETATM" else "A"
            chain.add_residue(_new_res)
            residue = chain[-1]           # re-bind to stored copy
            last_res_key = res_key

        at = gemmi.Atom()
        at.name = atom_name[:4]
        if element:
            try:
                at.element = gemmi.Element(element)
            except Exception:
                at.element = gemmi.Element("X")
        at.pos = gemmi.Position(cx, cy, cz)
        at.occ = occ
        at.b_iso = b
        try:
            at.serial = int(clean(serial_col[i]) or "0")
        except ValueError:
            pass
        residue.add_atom(at)               # Atom does not need re-bind; we never mutate it after

    if len(st) == 0 or _structure_is_empty(st):
        raise ValueError("no atoms recovered from _atom_site")
    return st


def _append_sdf_as_chain(st, sdf_bytes: bytes, chain_id: str = "L",
                          res_name: str = "LIG", res_num: int = 1) -> None:
    """Append an SDF molecule to a gemmi Structure as a HETATM chain.

    Uses RDKit to read atoms + coordinates, then materializes a gemmi chain.
    """
    gemmi = _import_gemmi()
    Chem = _import_rdkit()

    mol = Chem.MolFromMolBlock(sdf_bytes.decode("utf-8", errors="ignore"),
                               sanitize=False, removeHs=False)
    if mol is None:
        log.warning("could not parse SDF — ligand not appended")
        return
    conf = mol.GetConformer()

    # Avoid colliding with an existing chain id
    existing = {c.name for c in st[0]} if len(st) else set()
    cid = chain_id
    i = 0
    while cid in existing:
        i += 1
        cid = f"{chain_id}{i}"

    model = st[0]

    # Add empty chain first, then re-bind and add an empty residue, then re-bind
    # again and add atoms. This ordering matters: in gemmi, add_* copies the
    # passed object, so we must always mutate via the live reference returned
    # by indexing the parent.
    _empty_chain = gemmi.Chain(cid)
    model.add_chain(_empty_chain)
    chain = model[-1]

    _empty_res = gemmi.Residue()
    _empty_res.name = res_name
    _empty_res.seqid = gemmi.SeqId(res_num, " ")
    _empty_res.het_flag = "H"
    chain.add_residue(_empty_res)
    residue = chain[-1]

    for atom_idx in range(mol.GetNumAtoms()):
        a = mol.GetAtomWithIdx(atom_idx)
        pos = conf.GetAtomPosition(atom_idx)
        at = gemmi.Atom()
        at.name = f"{a.GetSymbol()}{atom_idx + 1}"[:4]
        at.element = gemmi.Element(a.GetSymbol())
        at.pos = gemmi.Position(pos.x, pos.y, pos.z)
        at.occ = 1.0
        at.b_iso = 20.0
        residue.add_atom(at)


# ---------------------------------------------------------------------------
# Fallback rebuild via Biopython (used if gemmi import fails)
# ---------------------------------------------------------------------------
def _rebuild_with_biopython(
    structure_bytes: bytes, fmt: str, pdb_id: str,
    extra_ligand_sdf: bytes | None = None,
) -> SanitizedCIF:
    PDBParser, MMCIFParser, MMCIFIO = _import_biopython()

    with tempfile.NamedTemporaryFile(suffix=f".{fmt}", delete=True) as tmp:
        tmp.write(structure_bytes)
        tmp.flush()
        parser = PDBParser(QUIET=True) if fmt == "pdb" else MMCIFParser(QUIET=True)
        st = parser.get_structure(pdb_id, tmp.name)

    # Biopython doesn't rebuild entity loops. We accept that — OST will
    # auto-infer from coordinates when label_asym_id is present in a
    # minimal form, which MMCIFIO writes.
    io_writer = MMCIFIO()
    io_writer.set_structure(st)
    buf = io.StringIO()
    io_writer.save(buf)
    text = buf.getvalue()

    ligand_chain_ids: list[str] = []
    n_lig = 0
    for model in st:
        for chain in model:
            lig_here = [res for res in chain
                        if res.id[0].strip() and res.id[0] != "W"]
            if lig_here and len(lig_here) == len(list(chain)):
                ligand_chain_ids.append(chain.id)
                n_lig += len(lig_here)
        break

    n_chains = len(list(next(iter(st)))) if len(list(st)) else 0
    return SanitizedCIF(
        text=text, pdb_id=pdb_id.upper(),
        n_chains=n_chains, n_ligand_residues=n_lig,
        ligand_chain_ids=tuple(ligand_chain_ids),
    )


# ---------------------------------------------------------------------------
# Public entrypoint
# ---------------------------------------------------------------------------
def sanitize(
    structure_path: Path,
    pdb_id: str,
    ligand_path: Path | None = None,
) -> SanitizedCIF:
    """Convert any supported input into a clean, OST-ready mmCIF (in memory).

    Parameters
    ----------
    structure_path : Path
        The predicted structure. May be .pdb, .cif, or .mmcif.
    pdb_id : str
        Used as the mmCIF data_ block name and entry.id.
    ligand_path : Path | None
        Optional separate SDF (DynamicBind emits ligand as SDF alongside a
        receptor PDB). When present, it's merged as its own HETATM chain.

    Raises
    ------
    ValueError
        If no backend (gemmi or biopython) can parse the structure.
    """
    structure_path = Path(structure_path)
    if not structure_path.exists():
        raise FileNotFoundError(structure_path)

    ext = structure_path.suffix.lower().lstrip(".")
    if ext in ("mmcif", "cif"):
        fmt = "cif"
    elif ext in ("pdb", "ent"):
        fmt = "pdb"
    else:
        raise ValueError(f"unsupported structure extension: .{ext}")

    structure_bytes = structure_path.read_bytes()
    sdf_bytes = ligand_path.read_bytes() if ligand_path else None

    # Try gemmi first — it rebuilds entity metadata properly.
    try:
        return _rebuild_with_gemmi(structure_bytes, fmt, pdb_id, sdf_bytes)
    except ImportError:
        log.info("gemmi not available; falling back to biopython")
    except Exception as e:
        log.warning("gemmi sanitize failed (%s); falling back to biopython", e)

    return _rebuild_with_biopython(structure_bytes, fmt, pdb_id, sdf_bytes)
