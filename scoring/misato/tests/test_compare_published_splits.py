"""Unit tests for the MISATO split-ID audit (no network or MD archive needed)."""
import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_published_splits import checked_list, compare, read_manifest, summarize  # noqa: E402


class PublishedSplitComparisonTests(unittest.TestCase):
    def test_list_validation_and_duplicate_rejection(self):
        data = b"10GS\n11GS\n"
        self.assertEqual(checked_list(data, "test", hashlib.md5(data).hexdigest()), {"10GS", "11GS"})
        with self.assertRaisesRegex(ValueError, "MD5"):
            checked_list(data, "test", "0" * 32)
        duplicate = b"10GS\n10GS\n"
        with self.assertRaisesRegex(ValueError, "duplicate"):
            checked_list(duplicate, "test", hashlib.md5(duplicate).hexdigest())

    def test_manifest_preserves_suffixed_ids_and_validates_flags(self):
        fields = ["target_id", "equibind_status", "diffdock_status", "in_valid_primary_cohort"]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "manifest.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow({"target_id": "1G42_A", "equibind_status": "candidate", "diffdock_status": "candidate", "in_valid_primary_cohort": "True"})
            self.assertIn("1G42_A", read_manifest(path))
            content = path.read_text(encoding="utf-8").replace("True", "False")
            path.write_text(content, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "disagrees"):
                read_manifest(path)

    def test_comparison_keeps_missing_and_outside_ids_distinct(self):
        membership = {"10GS": "train", "11GS": "val", "16PK": "test"}
        targets = {
            "10GS": {"equibind_status": "candidate", "diffdock_status": "candidate"},
            "11GS": {"equibind_status": "candidate", "diffdock_status": "not_delivered"},
            "1G42_A": {"equibind_status": "candidate", "diffdock_status": "candidate"},
        }
        rows = {row["target_id"]: row for row in compare(membership, targets)}
        self.assertEqual(rows["10GS"]["coverage_category"], "paired_primary_candidate")
        self.assertEqual(rows["11GS"]["coverage_category"], "equibind_only_delivered")
        self.assertEqual(rows["16PK"]["coverage_category"], "neither_delivered")
        self.assertEqual(rows["1G42_A"]["coverage_category"], "delivered_outside_published_md_splits")
        summaries = {row["split"]: row for row in summarize(list(rows.values()))}
        self.assertEqual(summaries["train"]["paired_primary_candidates"], 1)
        self.assertEqual(summaries["test"]["neither_delivered"], 1)


if __name__ == "__main__":
    unittest.main()
