import unittest

from litrapid.biosynthesis import (
    PLATFORMS,
    STOP_FREQUENCY,
    LibrarySpec,
    biosynthesis_profile,
    compare_platforms,
    disulfide_closure,
    display_probability_table,
    elongation_fidelity,
    intein_splicing,
    simulate_biosynthesis,
    stop_free_fraction,
    thioether_closure,
)


class TestGeneticCode(unittest.TestCase):
    def test_stop_frequencies_match_enumeration(self):
        # B excludes adenine from the third position, so NNB carries only amber.
        self.assertAlmostEqual(STOP_FREQUENCY["NNK"], 1 / 32)
        self.assertAlmostEqual(STOP_FREQUENCY["NNS"], 1 / 32)
        self.assertAlmostEqual(STOP_FREQUENCY["NNN"], 3 / 64)
        self.assertAlmostEqual(STOP_FREQUENCY["NNB"], 1 / 48)
        self.assertEqual(STOP_FREQUENCY["trimer"], 0.0)

    def test_stop_free_fraction_of_a_12mer_nnk_library(self):
        self.assertAlmostEqual(stop_free_fraction("NNK", 12), (31 / 32) ** 12, places=10)
        self.assertEqual(stop_free_fraction("trimer", 12), 1.0)


class TestChemistry(unittest.TestCase):
    def test_thioether_requires_a_cysteine(self):
        self.assertEqual(thioether_closure("AAAAAAAA")["efficiency"], 0.0)
        self.assertGreater(thioether_closure("AAAAAAAC")["efficiency"], 0.5)

    def test_thioether_ring_size_has_an_optimum(self):
        # closure through a 9-residue ring beats both a 2-mer and a 25-mer
        small = thioether_closure("AC")["efficiency"]
        good = thioether_closure("AAAAAAAAC")["efficiency"]
        large = thioether_closure("A" * 24 + "C")["efficiency"]
        self.assertGreater(good, small)
        self.assertGreater(good, large)

    def test_second_cysteine_flags_an_isomer_mixture(self):
        result = thioether_closure("AAACAAAC")
        self.assertTrue(result["isomer_mixture"])
        self.assertEqual(result["cysteines"], 2)

    def test_disulfide_needs_a_pair(self):
        self.assertEqual(disulfide_closure("AAACAAAA")["efficiency"], 0.0)
        self.assertGreater(disulfide_closure("ACAAAAAC")["efficiency"], 0.0)
        # three thiols give a pairing mixture, so the yield of any one ring falls
        self.assertLess(disulfide_closure("ACAAACAC")["efficiency"],
                        disulfide_closure("ACAAAAAC")["efficiency"])

    def test_intein_needs_a_plus_one_nucleophile(self):
        with_cys = intein_splicing("CGSRLF")["efficiency"]
        without = intein_splicing("AGSRLF")["efficiency"]
        self.assertGreater(with_cys, 0.0)
        self.assertLess(without, with_cys / 10)

    def test_proline_after_the_catalytic_cysteine_blocks_splicing(self):
        self.assertLess(intein_splicing("CPSRLF")["efficiency"],
                        0.3 * intein_splicing("CGSRLF")["efficiency"])


class TestTranslation(unittest.TestCase):
    def test_consecutive_n_methyls_cost_more_than_isolated_ones(self):
        isolated = elongation_fidelity("AmAmAAAA")
        adjacent = elongation_fidelity("AAmmAAAA")
        self.assertEqual("AmAmAAAA".count("m"), "AAmmAAAA".count("m"))
        self.assertLess(adjacent, isolated)


class TestProfile(unittest.TestCase):
    def test_cell_based_platforms_reject_noncanonical_residues(self):
        profile = biosynthesis_profile("AmAAAAAC", "phage_m13")
        self.assertEqual(profile["display_probability"], 0.0)
        self.assertEqual(profile["limiting_step"], "genetic_code")

    def test_mrna_display_accepts_noncanonical_residues(self):
        profile = biosynthesis_profile("yAmAAAAC", "mrna_rapid")
        self.assertGreater(profile["display_probability"], 0.0)

    def test_every_platform_reports_a_limiting_step(self):
        for key in PLATFORMS:
            profile = biosynthesis_profile("CGSRLFAC", key)
            self.assertIn("limiting_step", profile)
            self.assertTrue(profile["limiting_reason"])
            self.assertGreaterEqual(profile["display_probability"], 0.0)
            self.assertLessEqual(profile["display_probability"], 1.0)

    def test_dibasic_motifs_are_penalised_only_on_bacterial_routes(self):
        clean = biosynthesis_profile("ACGSDQSAEC", "phage_m13")["display_probability"]
        dibasic = biosynthesis_profile("ACGSKRSAEC", "phage_m13")["display_probability"]
        self.assertLess(dibasic, clean)

    def test_unknown_residue_is_rejected(self):
        with self.assertRaises(ValueError):
            biosynthesis_profile("ACGZZ", "phage_m13")


class TestFunnel(unittest.TestCase):
    def setUp(self):
        self.spec = LibrarySpec(length=12, sample_size=500, seed=7)

    def test_molecules_and_distinct_decrease_monotonically(self):
        for key in PLATFORMS:
            result = simulate_biosynthesis(key, self.spec)
            molecules = [row["molecules"] for row in result["funnel"]]
            distinct = [row["distinct"] for row in result["funnel"]]
            for earlier, later in zip(molecules, molecules[1:]):
                self.assertLessEqual(later, earlier + 1e-9, key)
            for earlier, later in zip(distinct, distinct[1:]):
                self.assertLessEqual(later, earlier + 1e-9, key)

    def test_distinct_never_exceeds_molecules(self):
        for key in PLATFORMS:
            result = simulate_biosynthesis(key, self.spec)
            for row in result["funnel"]:
                self.assertLessEqual(row["distinct"], row["molecules"] + 1e-9,
                                     f"{key}:{row['token']}")

    def test_transformation_ceiling_binds_on_cell_based_platforms(self):
        result = simulate_biosynthesis("phage_m13", self.spec)
        transformation = next(r for r in result["funnel"] if r["token"] == "XFORM")
        self.assertTrue(transformation["ceiling_applied"])
        self.assertLessEqual(transformation["distinct"],
                             PLATFORMS["phage_m13"].step("XFORM").ceiling + 1e-9)

    def test_cell_free_route_keeps_more_diversity_than_cell_based_routes(self):
        mrna = simulate_biosynthesis("mrna_rapid", self.spec)["final_distinct"]
        for key in ("phage_m13", "bacteria_siclopps", "bacteria_ecpx", "yeast_aga2"):
            self.assertGreater(mrna, simulate_biosynthesis(key, self.spec)["final_distinct"],
                               key)

    def test_reproducible_under_a_fixed_seed(self):
        first = simulate_biosynthesis("phage_m13", self.spec)["final_distinct"]
        second = simulate_biosynthesis("phage_m13", self.spec)["final_distinct"]
        self.assertEqual(first, second)

    def test_residue_bias_is_reported_for_every_canonical_residue(self):
        result = simulate_biosynthesis("phage_m13", self.spec)
        self.assertEqual(len(result["residue_bias_log2"]), 20)

    def test_invalid_library_spec_is_rejected(self):
        with self.assertRaises(ValueError):
            LibrarySpec(length=2).validate()
        with self.assertRaises(ValueError):
            LibrarySpec(codon_design="NNQ").validate()


class TestComparison(unittest.TestCase):
    def test_comparison_covers_every_platform(self):
        report = compare_platforms(LibrarySpec(length=12, sample_size=300, seed=3))
        self.assertEqual(len(report["comparison"]), len(PLATFORMS))
        self.assertEqual(report["comparison"][0]["platform"], "mrna_rapid")

    def test_display_probability_table_is_sequence_resolved(self):
        # a single cysteine cannot be oxidised into a ring, so the phage route
        # cannot make it at all; the two-cysteine sequences differ from each
        # other through export, proteolysis and toxicity, not through chemistry
        sequences = ["ACSGDQSAEC", "CKRKRKRKRC", "AAAAAAAAAC"]
        table = display_probability_table(sequences, "phage_m13")
        self.assertEqual(len(table), 3)
        self.assertEqual(table["AAAAAAAAAC"], 0.0)
        self.assertGreater(table["ACSGDQSAEC"], 0.0)
        self.assertGreater(table["ACSGDQSAEC"], table["CKRKRKRKRC"])


class TestSelectionIntegration(unittest.TestCase):
    """The scalar `fusion_probability` must become sequence-resolved on request."""

    CANDIDATES = [
        # identical affinity, different manufacturability: a 9-residue thioether
        # ring versus a sequence whose only cysteine sits at the N-terminus,
        # leaving a strained two-atom ring
        {"candidate_id": "GOOD", "peptide_sequence": "YSWDGFTLC", "kd_nm": 10.0},
        {"candidate_id": "STRAINED", "peptide_sequence": "CSWDGFTLY", "kd_nm": 10.0},
    ]

    def test_default_behaviour_is_unchanged(self):
        from litrapid.mrna_display import SimulationParameters, simulate_mrna_display
        result = simulate_mrna_display(self.CANDIDATES,
                                       SimulationParameters(rounds=2, ngs_reads=1000))
        fusion = {row["candidate_id"]: row["fusion_probability"]
                  for row in result["trajectories"] if row["round"] == 1}
        self.assertEqual(fusion["GOOD"], fusion["STRAINED"])

    def test_mechanistic_mode_separates_equal_affinity_candidates(self):
        from litrapid.mrna_display import SimulationParameters, simulate_mrna_display
        result = simulate_mrna_display(
            self.CANDIDATES,
            SimulationParameters(rounds=2, ngs_reads=1000,
                                 biosynthesis_platform="mrna_rapid"))
        fusion = {row["candidate_id"]: row["fusion_probability"]
                  for row in result["trajectories"] if row["round"] == 1}
        self.assertGreater(fusion["GOOD"], fusion["STRAINED"])
        self.assertIn("mrna_rapid", result["model"])

    def test_unknown_platform_is_rejected(self):
        from litrapid.mrna_display import SimulationParameters
        with self.assertRaises(ValueError):
            SimulationParameters(biosynthesis_platform="nanopore_display").validate()


if __name__ == "__main__":
    unittest.main()
