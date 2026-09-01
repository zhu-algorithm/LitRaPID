import unittest

from litrapid.platform_server import api_design, api_export, api_feedback, api_meta


class PlatformServerTests(unittest.TestCase):
    def setUp(self):
        self.candidates = [
            {"candidate_id": "C1", "canonical_sequence": "CADEFGHIC", "predicted_affinity_score": "0.8"},
            {"candidate_id": "C2", "canonical_sequence": "CWQKPGTLC", "predicted_affinity_score": "0.7"},
        ]

    def test_meta_lists_three_backends(self):
        meta = api_meta()
        self.assertEqual(len(meta["backends"]), 3)
        self.assertEqual(meta["literature_route"]["backend"], "yeast_gpi_macrocycle")

    def test_full_api_cycle(self):
        panel = api_export({"backend": "yeast_gpi_macrocycle", "candidates": self.candidates})["rows"]
        self.assertEqual(len(panel), 2)
        results = [
            {"candidate_id": "C1", "round": "1", "ngs_input_count": "100", "ngs_selected_count": "800", "target_binding_mfi": "900", "display_mfi": "450"},
            {"candidate_id": "C2", "round": "1", "ngs_input_count": "100", "ngs_selected_count": "200", "target_binding_mfi": "500", "display_mfi": "500"},
        ]
        ranked = api_feedback({"panel": panel, "results": results})["rows"]
        self.assertIn("posterior_score", ranked[0])
        designed = api_design({"backend": "yeast_gpi_macrocycle", "ranked": ranked, "limit": 6, "parents": 1, "variants_per_parent": 4})["rows"]
        self.assertTrue(designed)
        self.assertIn("dna_insert", designed[0])


if __name__ == "__main__":
    unittest.main()
