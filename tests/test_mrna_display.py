import unittest

from litrapid.mrna_display import SimulationParameters, simulate_mrna_display
from litrapid.platform_server import api_simulate_mrna


CANDIDATES = [
    {"candidate_id": "HIGH", "canonical_sequence": "CAAAAAAAC", "kd_nm": 1},
    {"candidate_id": "MID", "canonical_sequence": "CGGGGGGGC", "kd_nm": 100},
    {"candidate_id": "LOW", "canonical_sequence": "CSSSSSSSC", "kd_nm": 10_000},
]


class MrnaDisplayTests(unittest.TestCase):
    def test_affinity_drives_multiround_enrichment(self):
        result = simulate_mrna_display(CANDIDATES, SimulationParameters(rounds=5, ngs_reads=10_000, seed=7))
        self.assertEqual(result["final_ranking"][0]["candidate_id"], "HIGH")
        high = [r for r in result["trajectories"] if r["candidate_id"] == "HIGH"]
        self.assertGreater(high[-1]["ngs_frequency"], high[0]["ngs_frequency"])

    def test_seed_is_reproducible_and_frequencies_normalize(self):
        params = SimulationParameters(rounds=3, ngs_reads=2_000, seed=99)
        first = simulate_mrna_display(CANDIDATES, params)
        second = simulate_mrna_display(CANDIDATES, params)
        self.assertEqual(first, second)
        total = sum(row["ngs_frequency"] for row in first["final_ranking"])
        self.assertAlmostEqual(total, 1.0, places=7)

    def test_api_accepts_parameter_overrides(self):
        result = api_simulate_mrna({"candidates": CANDIDATES, "parameters": {"rounds": 2, "seed": 4}})
        self.assertEqual(len(result["round_summaries"]), 2)


if __name__ == "__main__":
    unittest.main()
