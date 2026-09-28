"""Regression tests for conservative multi-residue recovery and merging."""
from __future__ import annotations

import unittest
from types import SimpleNamespace

from rdkit import Chem

from scoring.misato.merge_recovery_scores import merge
from scoring.misato.score_multiresidue_recovery_wsl import (
    graph_valid_copies,
    same_indexed_heavy_graph,
)


class MultiresidueRecoveryTests(unittest.TestCase):
    def test_atom_order_must_match_before_reusing_native_mapping(self):
        first = Chem.MolFromSmiles("CCO")
        reordered = Chem.RenumberAtoms(first, [2, 1, 0])
        self.assertTrue(same_indexed_heavy_graph(first, Chem.Mol(first)))
        self.assertFalse(same_indexed_heavy_graph(first, reordered))

    def test_disconnected_native_cannot_match_connected_pose(self):
        pose = Chem.MolFromSmiles("CC")
        candidate = lambda distance: {
            "kind": "nonpolymer_pair", "label": "A+B", "subchains": [],
            "atoms": [SimpleNamespace(element=SimpleNamespace(name="C"),
                                      pos=SimpleNamespace(x=x, y=0.0, z=0.0))
                      for x in (0.0, distance)],
        }
        self.assertEqual(graph_valid_copies(pose, [candidate(10.0)]), [])
        self.assertEqual(len(graph_valid_copies(pose, [candidate(1.5)])), 1)

    def test_merge_preserves_v3_and_separates_direct_fallback(self):
        v3 = [
            {"target_id": "AAAA", "method": "diffdock", "status": "scored_both",
             "bisy_rmsd_angstrom": "1.000000", "lddt_pli": "0.9"},
            {"target_id": "BBBB", "method": "diffdock", "status": "excluded_reference_or_graph_check"},
            {"target_id": "CCCC", "method": "diffdock", "status": "excluded_reference_or_graph_check"},
        ]
        allcopies = [
            {"target_id": "AAAA", "method": "diffdock", "status": "scored_both",
             "bisy_rmsd_angstrom": "1.1", "lddt_pli": "0.8"},
            {"target_id": "BBBB", "method": "diffdock", "status": "excluded_reference_or_graph_check",
             "preaudit_status": "composite_reference_review"},
            {"target_id": "CCCC", "method": "diffdock", "status": "excluded_reference_or_graph_check",
             "preaudit_status": "composite_reference_review"},
        ]
        recovery = [
            {"target_id": "BBBB", "method": "diffdock", "status": "scored_both",
             "bisy_rmsd_angstrom": "2.0", "lddt_pli": "0.5", "candidate_kind": "short_polymer",
             "assigned_copy_label": "B:B", "n_target_copies": "2"},
            {"target_id": "CCCC", "method": "diffdock", "status": "failed_scoring",
             "error": "test failure"},
        ]
        audit = [{"target_id": "BBBB", "method": "diffdock",
                  "graph_status": "exact_full_graph"},
                 {"target_id": "CCCC", "method": "diffdock", "direct_rmsd_angstrom": "3.0",
                  "candidate_kind": "branched_glycan", "candidate_label": "C:C",
                  "graph_status": "exact_full_graph"}]
        rows = merge(v3, allcopies, recovery, [], audit, [], expected_count=3)
        self.assertEqual(rows[0]["pose_rmsd_angstrom"], "1.000000")
        self.assertEqual(rows[1]["pose_rmsd_angstrom"], "2.0")
        self.assertEqual(rows[2]["pose_rmsd_angstrom"], "")
        self.assertEqual(rows[2]["direct_rmsd_fallback_angstrom"], "3.0")
        self.assertEqual(rows[2]["lddt_pli"], "")
        with self.assertRaisesRegex(ValueError, "Recovery rows incomplete"):
            merge(v3, allcopies, recovery[:1], [], audit, [], expected_count=3)
        retried = [{"target_id": "CCCC", "method": "diffdock", "status": "scored_both",
                    "bisy_rmsd_angstrom": "2.5", "lddt_pli": "0.4",
                    "candidate_kind": "branched_glycan", "assigned_copy_label": "C:C",
                    "n_target_copies": "1"}]
        retried_rows = merge(v3, allcopies, recovery, [], audit, [], expected_count=3,
                             recovery_retry=retried)
        self.assertEqual(retried_rows[2]["pose_rmsd_angstrom"], "2.5")
        self.assertEqual(retried_rows[2]["metric_source"], "ost_pseudo_multiresidue_retry")


if __name__ == "__main__":
    unittest.main()
