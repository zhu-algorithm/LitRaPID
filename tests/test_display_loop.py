import csv
import tempfile
import unittest
from pathlib import Path

from litrapid.display_loop import BACKENDS, design_next_round, displayability, export_panel, feedback_update, macrocycle_topology


class DisplayLoopTests(unittest.TestCase):
    def test_noncanonical_is_synthesis_only(self):
        score, flags, status = displayability("ACD[BETA]FG", BACKENDS["yeast_aga2"])
        self.assertEqual(status, "SYNTHESIS_ONLY")
        self.assertEqual(score, 0.0)
        self.assertIn("contains_noncanonical_or_invalid_residue", flags)

    def test_odd_cysteine_is_flagged_for_p3(self):
        _, flags, _ = displayability("ACDEFGH", BACKENDS["m13_p3"])
        self.assertIn("odd_cysteine_count_for_p3", flags)

    def test_gpi_macrocycle_backend_and_topology(self):
        topology, score, positions = macrocycle_topology("CADEFGHIC")
        self.assertEqual(topology, "one_ring")
        self.assertGreater(score, 0.8)
        self.assertEqual(positions, [1, 9])
        display_score, flags, status = displayability("CADEFGHIC", BACKENDS["yeast_gpi_macrocycle"])
        self.assertGreater(display_score, 0.5)
        self.assertEqual(status, "DIRECT_DISPLAY")
        self.assertNotIn("gpi_macrocycle_requires_at_least_two_cysteines", flags)

    def test_gpi_backend_requires_cysteines(self):
        _, flags, status = displayability("ASDFGHKLM", BACKENDS["yeast_gpi_macrocycle"])
        self.assertIn("gpi_macrocycle_requires_at_least_two_cysteines", flags)
        self.assertEqual(status, "DISPLAY_REVIEW")

    def test_export_and_feedback(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            panel = Path(tmp) / "panel.csv"
            next_round = Path(tmp) / "next.csv"
            export_panel(root / "examples" / "candidates.csv", "yeast_aga2", panel, 96)
            feedback_update(panel, root / "examples" / "display_results.csv", next_round, 3)
            with next_round.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertTrue(rows)
            self.assertIn("posterior_score", rows[0])
            self.assertIn("binding_display_ratio", rows[0])
            self.assertLessEqual(len(rows), 3)

            designed = Path(tmp) / "designed.csv"
            design_next_round(next_round, "yeast_aga2", designed, 10, 2, 4)
            with designed.open(encoding="utf-8") as handle:
                variants = list(csv.DictReader(handle))
            self.assertTrue(variants)
            self.assertIn("parent_id", variants[0])
            self.assertIn("structure_hypothesis", variants[0])
            self.assertTrue(all(v["dna_insert"] for v in variants if v["display_status"] != "SYNTHESIS_ONLY"))


if __name__ == "__main__":
    unittest.main()
