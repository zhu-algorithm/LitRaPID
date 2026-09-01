import unittest
from litrapid.digital_twin import infer_latent_fitness, optimize_protocol
from litrapid.mrna_display import SimulationParameters, simulate_mrna_display

class DigitalTwinTests(unittest.TestCase):
    def setUp(self):
        self.c=[{"candidate_id":"A","canonical_sequence":"CAAAAAAC","kd_nm":1},{"candidate_id":"B","canonical_sequence":"CGGGGGGC","kd_nm":100}]
    def test_bias_aware_inference(self):
        sim=simulate_mrna_display(self.c,SimulationParameters(rounds=3,ngs_reads=1000,seed=2,pcr_bias_sigma=.5))
        self.assertEqual({x["candidate_id"] for x in infer_latent_fitness(sim["trajectories"])},{"A","B"})
    def test_protocol_optimizer_returns_pareto_set(self):
        out=optimize_protocol(self.c,SimulationParameters(rounds=2,ngs_reads=500,seed=3),[5,20],[.4,.8],[.1,.5])
        self.assertTrue(out); self.assertTrue(all("effective_diversity" in x for x in out))
