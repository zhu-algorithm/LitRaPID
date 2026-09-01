import unittest
from litrapid.pdl1_bridge import adapt_pdl1_report,validate_pdl1_in_silico_display
from litrapid.mrna_display import SimulationParameters

def report(target="PD-L1"):
    return {"epitope_design":{"profile_id":"pdl1_test","target":target},"candidates":[
      {"id":"CP-1","sequence":"FWYPG","cyclization":"head-to-tail (conceptual)","priority_score":.9,"binding":{"pd_l1_binding_proxy":.9,"selectivity_score":.7},"epitope":{"antibody_epitope_compatibility_proxy":.8},"admet":{"synthesizability":.8}},
      {"id":"CP-2","sequence":"ADEKGC","cyclization":"head-to-tail (conceptual)","priority_score":.6,"binding":{"pd_l1_binding_proxy":.5,"selectivity_score":.4},"epitope":{"antibody_epitope_compatibility_proxy":.5},"admet":{"synthesizability":.7}}]}
class BridgeTests(unittest.TestCase):
    def test_chemistry_conversion_is_explicit(self):
        rows=adapt_pdl1_report(report())["candidates"]
        self.assertEqual(rows[0]["canonical_sequence"],"FWYPGC"); self.assertIn("appended",rows[0]["chemistry_conversion"])
    def test_wrong_target_profile_is_rejected(self):
        with self.assertRaises(ValueError): adapt_pdl1_report(report("PD-1"))
    def test_end_to_end_validation(self):
        out=validate_pdl1_in_silico_display(report(),SimulationParameters(rounds=2,ngs_reads=500,seed=8))
        self.assertEqual(len(out["simulation"]["round_summaries"]),2)
    def test_versioned_upstream_contract(self):
        payload={"schema":"litrapid.pdl1-report.v1","pdl1_report":report()}
        self.assertEqual(adapt_pdl1_report(payload)["target"],"PD-L1")
