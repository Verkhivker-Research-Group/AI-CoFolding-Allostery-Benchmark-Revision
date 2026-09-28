"""Static-QS, gnina parsing, and missing-pocket-RMSD regression tests."""
from __future__ import annotations

import unittest

import numpy as np

from scoring.misato.compute_static_qs_wsl import contact_qs, residue_min_distances
from scoring.misato.merge_static_metrics import merge
from scoring.misato.rescore_gnina_wsl import parse_cnnscores


class StaticMetricsTests(unittest.TestCase):
    def test_static_qs_controls_and_continuous_score(self):
        receptor = np.array([[0., 0., 0.], [1., 0., 0.], [10., 0., 0.]])
        residue_index = np.array([0, 0, 1])
        native = residue_min_distances(receptor, residue_index, 2,
                                       np.array([[0., 1., 0.]]), contact_d=5)
        moved = residue_min_distances(receptor, residue_index, 2,
                                      np.array([[100., 1., 0.]]), contact_d=5)
        partial = residue_min_distances(receptor, residue_index, 2,
                                        np.array([[5.5, 1., 0.]]), contact_d=5)
        self.assertEqual(contact_qs(native, native, contact_d=5), 1.0)
        self.assertEqual(contact_qs(native, moved, contact_d=5), 0.0)
        self.assertGreater(contact_qs(native, partial, contact_d=5), 0.0)
        self.assertLess(contact_qs(native, partial, contact_d=5), 1.0)

    def test_gnina_score_parsing_requires_one_per_pose(self):
        output = "Affinity: -6.0\nCNNscore: 0.2 \nCNNscore: 0.85\n"
        self.assertEqual(parse_cnnscores(output, 2), [0.2, 0.85])
        with self.assertRaisesRegex(ValueError, "expected_1_CNNscores_got_2"):
            parse_cnnscores(output, 1)
        with self.assertRaisesRegex(ValueError, "out_of_range"):
            parse_cnnscores("CNNscore: 1.1\n", 1)

    def test_confidence_columns_and_pocket_rmsd_missing(self):
        scores = [
            {"target_id": "AAAA", "method": "diffdock", "score_status": "scored_both",
             "pose_rmsd_angstrom": "1.25", "lddt_pli": "0.9"},
            {"target_id": "AAAA", "method": "equibind", "score_status": "scored_both",
             "pose_rmsd_angstrom": "2.25", "lddt_pli": "0.5"},
        ]
        qs = [{"target_id": "AAAA", "method": method, "status": "scored",
               "static_receptor_qs": value} for method, value in
              (("diffdock", "0.8"), ("equibind", "0.4"))]
        gnina = [{"target_id": "AAAA", "method": method, "status": "scored",
                  "cnnscore": value, "gnina_version": "gnina v1.3.3"} for method, value in
                 (("diffdock", "0.7"), ("equibind", "0.6"))]
        records = [
            {"target_id": "AAAA", "method": "diffdock", "is_primary": "True",
             "sha256": "abc", "confidence": ""},
            {"target_id": "AAAA", "method": "diffdock", "is_primary": "False",
             "sha256": "abc", "confidence": "-0.71"},
        ]
        result = merge(scores, qs, gnina, records, expected_count=2)
        by_method = {row["method"]: row for row in result}
        self.assertEqual(by_method["diffdock"]["confidence"], "-0.71")
        self.assertEqual(by_method["equibind"]["confidence"], "")
        self.assertEqual(by_method["equibind"]["rescore_confidence"], "0.6")
        self.assertEqual(by_method["diffdock"]["static-receptor QS"], "0.8")
        self.assertTrue(all(row["pocket_rmsd_angstrom"] == "" for row in result))
        self.assertTrue(all(row["receptor_mode"] == "static_crystal" for row in result))
        with self.assertRaisesRegex(ValueError, "static QS keys differ"):
            merge(scores, qs[:1], gnina, records, expected_count=2)


if __name__ == "__main__":
    unittest.main()
