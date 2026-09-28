"""Small deterministic checks for cohort denominators and paired scoring."""
from __future__ import annotations

import unittest
from unittest.mock import patch

from scoring.misato.build_methods_accounting import build_accounting
from scoring.misato.export_delivered_cohort import export
from scoring.misato.score_static_docking_wsl import initial_row, resolve_path
from scoring.misato.selected_ligand_metrics import choose_native, write_heavy_only_pose
from scoring.misato.summarize_static_docking import rmsd_disagreements, summarize


class MethodsAccountingTests(unittest.TestCase):
    def test_reported_and_observed_never_mix(self):
        manifest = [
            {"target_id": "AAAA", "diffdock_status": "candidate", "equibind_status": "candidate"},
            {"target_id": "BBBB", "diffdock_status": "not_delivered", "equibind_status": "candidate"},
        ]
        identity = [{"target_id": "AAAA", "disposition": "eligible_exact_graph"}]
        crystal = [
            {"target_id": "AAAA", "method": method, "reference_status": "provisional_pose_rmsd",
             "reference_selection": "single_candidate", "reference_residues": "A:LIG:1"}
            for method in ("diffdock", "equibind")
        ]
        # Production function intentionally verifies Lucas's observed pose totals;
        # override those only to exercise a tiny synthetic cohort.
        with patch("scoring.misato.build_methods_accounting.REPORTED", {
            "qm_hdf5_entries": 4, "qm_ligand_sdfs_written": 3,
            "qm_ligand_sanitization_failed": 1, "rcsb_proteins_saved": 3,
            "rcsb_proteins_missing": 1, "diffdock_input_pairs": 2,
            "diffdock_posed": 1, "diffdock_rdkit_skipped": 1,
            "equibind_posed": 2, "equibind_failed": 1,
            "equibind_failed_missing_protein": 1,
        }):
            result = build_accounting(manifest, identity, crystal, None)
        observed = result["evidence"]["directly_observed_delivered_files_and_reference_audit"]
        self.assertEqual(observed["paired_primary_pose_ids"], 1)
        self.assertEqual(observed["strict_single_reference_paired_ids"], 1)
        self.assertEqual(observed["equibind_only_primary_pose_ids"], 1)


class StaticScoringTests(unittest.TestCase):
    def test_portable_cohort_index_requires_complete_scores(self):
        manifest = [{"target_id": "AAAA", "diffdock_status": "candidate",
                     "equibind_status": "not_delivered", "diffdock_primary_path": r"C:\pose.sdf"}]
        audit = [{"target_id": "AAAA", "method": "diffdock",
                  "reference_selection": "single_candidate",
                  "reference_status": "provisional_pose_rmsd", "reference_residues": "A:LIG:1"}]
        score = [{"target_id": "AAAA", "method": "diffdock", "status": "scored_both",
                  "bisy_rmsd_angstrom": "1", "lddt_pli": "0.5"}]
        rows = export(manifest, audit, score, [])
        self.assertEqual(rows[0]["diffdock_primary_delivered"], "true")
        self.assertEqual(rows[0]["equibind_score_status"], "not_delivered")
        self.assertNotIn("C:", str(rows))
        with self.assertRaisesRegex(ValueError, "exactly cover"):
            export(manifest, audit, [], [])

    def test_native_selection_rejects_partial_same_named_copy(self):
        from rdkit import Chem
        from types import SimpleNamespace

        pose = Chem.MolFromSmiles("CCO")
        def residue(name, elements):
            return SimpleNamespace(name=name,
                                   atoms=[SimpleNamespace(element=element) for element in elements])
        partial = residue("LIG", ["C", "O"])
        exact = residue("LIG", ["C", "C", "O"])
        other = residue("COF", ["C", "C", "O"])
        self.assertIs(choose_native([partial, exact, other], "LIG", pose), exact)
        with self.assertRaisesRegex(ValueError, "selected_native_residue_count_2"):
            choose_native([exact, residue("LIG", ["C", "C", "O"])], "LIG", pose)

    def test_heavy_only_temporary_pose_preserves_heavy_coordinates(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from tempfile import TemporaryDirectory
        from pathlib import Path

        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        AllChem.EmbedMolecule(mol, randomSeed=5)
        with TemporaryDirectory() as temp:
            source, heavy_path = Path(temp) / "source.sdf", Path(temp) / "heavy.sdf"
            writer = Chem.SDWriter(str(source))
            writer.write(mol)
            writer.close()
            self.assertEqual(write_heavy_only_pose(source, heavy_path), 3)
            heavy = Chem.SDMolSupplier(str(heavy_path), removeHs=False)[0]
            self.assertEqual(heavy.GetNumAtoms(), 3)
            for index in range(3):
                original = mol.GetConformer().GetAtomPosition(index)
                actual = heavy.GetConformer().GetAtomPosition(index)
                self.assertAlmostEqual(original.x, actual.x, places=3)

    def test_windows_paths_are_resolved_for_wsl(self):
        from pathlib import Path
        root = Path("/mnt/c/repo")
        self.assertEqual(resolve_path(r"C:\Users\Ryan\file.sdf", root),
                         Path("/mnt/c/Users/Ryan/file.sdf"))
        self.assertEqual(resolve_path(r"misato_output\ref.cif", root),
                         root / "misato_output/ref.cif")

    def test_multicopy_reference_is_excluded(self):
        row = {"target_id": "AAAA", "method": "diffdock", "pose_path": "pose.sdf",
               "reference_path": "ref.cif", "reference_residues": "A:LIG:1",
               "reference_selection": "copy_by_diffdock_proximity",
               "reference_status": "provisional_pose_rmsd"}
        self.assertEqual(initial_row(row)["status"], "excluded_reference_or_graph_check")

    def test_missing_scores_do_not_enter_paired_denominator(self):
        audit = [{"target_id": "AAAA", "method": method} for method in ("diffdock", "equibind")]
        common = {"target_id": "AAAA", "reference_selection": "single_candidate",
                  "preaudit_status": "provisional_pose_rmsd", "reference_residue": "A:LIG:1",
                  "lddt_pli": "0.5", "status": "partial_or_no_metrics"}
        scores = [{**common, "method": "diffdock", "bisy_rmsd_angstrom": "1.0"},
                  {**common, "method": "equibind", "bisy_rmsd_angstrom": ""}]
        result, _ = summarize(scores, audit)
        broad = result["paired_broad_topology_matched_sensitivity"]
        self.assertEqual(broad["eligible_paired_ids_in_score_file"], 1)
        self.assertEqual(broad["both_bisy_rmsd_available"], 0)
        self.assertEqual(broad["both_lddt_pli_available"], 1)
        self.assertEqual(result["paired_primary_exact_prediction_graph"]["eligible_paired_ids_in_score_file"], 0)

    def test_every_large_frame_disagreement_is_retained(self):
        audit = [{"target_id": target, "method": "diffdock",
                  "pose_rmsd_angstrom_provisional": old}
                 for target, old in (("AAAA", "10"), ("BBBB", "20"), ("CCCC", "2"))]
        scores = [{"target_id": target, "method": "diffdock",
                   "bisy_rmsd_angstrom": new}
                  for target, new in (("AAAA", "1"), ("BBBB", "2"), ("CCCC", "1.7"))]
        rows = rmsd_disagreements(scores, audit)
        self.assertEqual([row["target_id"] for row in rows], ["BBBB", "AAAA"])


if __name__ == "__main__":
    unittest.main()
