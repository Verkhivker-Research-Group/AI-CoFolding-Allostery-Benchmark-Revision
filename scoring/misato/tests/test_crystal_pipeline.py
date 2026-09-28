"""Small offline tests; no existing benchmark files are accessed or modified."""
import csv
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fetch_crystal_cifs import requested_ids, valid_cif  # noqa: E402

try:
    import numpy as np
    from rdkit import Chem
    from audit_crystal_poses import crystal_candidates, direct_pose_rmsd
except ImportError:
    np = None
    Chem = None


class CrystalFetchTests(unittest.TestCase):
    def test_manifest_ids_are_deduplicated_without_losing_suffix_provenance(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "manifest.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["target_id", "equibind_status", "diffdock_status"])
                writer.writeheader()
                writer.writerow({"target_id": "1G42_A", "equibind_status": "candidate", "diffdock_status": "not_delivered"})
                writer.writerow({"target_id": "1G42_B", "equibind_status": "candidate", "diffdock_status": "candidate"})
                writer.writerow({"target_id": "10GS", "equibind_status": "not_delivered", "diffdock_status": "candidate"})
            self.assertEqual(requested_ids(path), ["10GS", "1G42"])

    def test_response_must_look_like_coordinate_cif(self):
        self.assertFalse(valid_cif(b"<html>" * 300, "10GS"))
        self.assertTrue(valid_cif(b"data_10GS\n" + b"_atom_site.\n" + b"x" * 1100, "10GS"))


@unittest.skipIf(Chem is None, "RDKit/numpy not installed")
class CrystalPoseTests(unittest.TestCase):
    @staticmethod
    def molecule():
        mol = Chem.MolFromSmiles("CCO")
        conformer = Chem.Conformer(3)
        for index, coords in enumerate(((0, 0, 0), (1, 0, 0), (2, 1, 0))):
            conformer.SetAtomPosition(index, coords)
        mol.AddConformer(conformer)
        return mol

    def test_rmsd_is_symmetry_aware_but_not_spatially_aligned(self):
        native = self.molecule()
        reordered = Chem.RenumberAtoms(native, [1, 0, 2])
        rmsd, mappings, issue = direct_pose_rmsd(reordered, native)
        self.assertAlmostEqual(rmsd, 0.0)
        self.assertGreaterEqual(mappings, 1)
        self.assertEqual(issue, "")
        shifted = Chem.Mol(reordered)
        for atom in range(shifted.GetNumAtoms()):
            position = shifted.GetConformer().GetAtomPosition(atom)
            shifted.GetConformer().SetAtomPosition(atom, (position.x + 3, position.y, position.z))
        rmsd, _, _ = direct_pose_rmsd(shifted, native)
        self.assertAlmostEqual(rmsd, 3.0)

    def test_composite_is_identified_not_force_scored(self):
        residues = [{"label": "A:AAA:1", "signature": Counter({6: 1})},
                    {"label": "A:BBB:2", "signature": Counter({8: 1})}]
        candidates, composite = crystal_candidates(residues, Counter({6: 1, 8: 1}))
        self.assertTrue(composite)
        self.assertEqual(len(candidates), 1)


if __name__ == "__main__":
    unittest.main()
