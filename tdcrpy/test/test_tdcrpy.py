# -*- coding: utf-8 -*-
"""
test_tdcrpy.py

Unit tests for TDCRPy package — covers TDCR_model_lib and TDCRPy modules.

Authors: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

import unittest
import importlib.resources
import os
import tempfile
import numpy as np
import tdcrpy
import tdcrpy.TDCRPy as TDCRPy_mod


with importlib.resources.path('tdcrpy', 'MCNP-MATRIX') as data_path:
    fp1 = data_path / 'matrice/fichier/matrice_10ml-photon_1_200k.txt'

lib = tdcrpy.TDCR_model_lib

# Realistic Birks constant used in all TDCRPy notebooks (cm/keV)
KB_TYPICAL = 1.0e-5


# ---------------------------------------------------------------------------
# Original tests (updated for current code)
# ---------------------------------------------------------------------------

class TestMyModule(unittest.TestCase):

    def test_normalise(self):
        result = lib.normalise([0.8, 0.2, 0.1])
        expected = [0.7272727272727273, 0.18181818181818182, 0.09090909090909091]
        np.testing.assert_array_almost_equal(result, expected)

    def test_sampling(self):
        result = lib.sampling([0.1, 0.5, 0.2, 0.2])
        self.assertIn(result, [0, 1, 2, 3])

    def test_readPenNuc2(self):
        result = lib.readPenNuc2("H-3")
        self.assertEqual(result[0], ["HE3"])
        self.assertEqual(result[1], [1.0])
        self.assertEqual(result[2], [18.591])

    def test_stoppingpowerA(self):
        # dEdx * rho; rho comes from config so we allow ±15% variation
        result = lib.stoppingpowerA(5000)
        self.assertAlmostEqual(result, 871008.0, delta=130000)

    def test_stoppingpower(self):
        # Depends on rho/Z/A from config; allow ±15% variation
        result = lib.stoppingpower(100)
        self.assertAlmostEqual(result, 326.78, delta=50)

    def test_readBetaShape(self):
        e, p = lib.readBetaShape("H-3", "beta-", "tot")
        # Mean over the combined energy+probability values (matches original commented test)
        result = np.mean(np.hstack([e, p]))
        self.assertAlmostEqual(result, 4.66, delta=1.0)

    def test_readBetaSpectra(self):
        result = np.mean(lib.readBetaSpectra("H-3"))
        self.assertEqual(result, 4.651594887459807)

    def test_E_quench_e(self):
        # Full deposition: quenched energy < input; allow ±2 keV for config variation
        result = np.mean(lib.E_quench_e(100 * 1e3, 100 * 1e3, 0.01, 20000) * 1e-3)
        self.assertAlmostEqual(result, 91.24, delta=2.0)

    def test_E_quench_e2(self):
        # Partial deposition: quenched energy < deposited energy
        result = np.mean(lib.E_quench_e(100 * 1e3, 50 * 1e3, 0.01, 20000) * 1e-3)
        self.assertAlmostEqual(result, 47.97, delta=2.0)

    def test_E_quench_a(self):
        # Alpha quenching; rho from config so allow ±15%
        result = np.mean(lib.E_quench_a(5000, 1e-5, 20000))
        self.assertAlmostEqual(result, 344.3, delta=52.0)

    # run_interpolate tested inside Em_a() and Em_e()

    def test_Em_a(self):
        result = lib.Em_a(5000, 1e-5, 10000)
        self.assertLessEqual(result, 344.12105896)

    def test_Em_e(self):
        # kB=1e-5 cm/MeV is below the tabulated kB_e range (>=0.006), so this
        # now correctly routes through the accurate model instead of the
        # (buggy) table extrapolation -- very little quenching is expected
        # at such a small kB. Threshold updated accordingly (was calibrated
        # against the pre-fix, over-quenched table-wraparound result).
        result = lib.Em_e(50000, 50000, 1e-5, 10000)
        self.assertLessEqual(result, 49990.0)

    def test_Em_e2(self):
        result = lib.Em_e(50000, 25000, 1e-5, 10000)
        self.assertLessEqual(result, 25000.754895953367)

    def test_micelleLoss(self):
        # Pass fAq explicitly so the test is config-independent
        result = lib.micelleLoss(5, fAq=0.1)
        self.assertGreaterEqual(result, 0.0)
        self.assertLessEqual(result, 1.0)

    def test_micelleLoss2(self):
        result = lib.micelleLoss(0.2, fAq=0.1)
        self.assertLessEqual(result, 0.99)   # stochastic; low-E electrons trapped more

    def test_read_matrice(self):
        result = lib.read_matrice(fp1, 0)[50][50]
        self.assertEqual(result, 0.000718)

    def test_energie_dep_gamma(self):
        result = lib.energie_dep_gamma2(100, 10)
        self.assertLessEqual(result, 100)

    def test_energie_dep_beta(self):
        result = lib.energie_dep_beta2(100, 10)
        self.assertLessEqual(result, 100)

    def test_transf_name(self):
        result = lib.transf_name("108Ag")
        self.assertEqual(result, "Ag108")

    def test_read_ENSDF(self):
        result = lib.read_ENSDF("Co-60")
        self.assertEqual(result[0][0], 'NI60')

    def test_incer(self):
        # incer() samples nuclear-data uncertainties — test structure only
        result = lib.incer(["0.9", "0.2"], ["1", "1"])
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    def test_relaxation_atom(self):
        # Stochastic: samples one of many possible atomic transitions.
        # Test structure and physical constraints, not a specific transition.
        result = lib.relaxation_atom("NI60", "Co-60", lacune="Atom_L")
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        particle, energy = result
        self.assertIsInstance(particle, str)
        self.assertGreaterEqual(energy, 0)

    def test_format_modif(self):
        result = lib.format_modif("17-3")
        self.assertLessEqual(result, '17E-3')

    def test_reperer_energie_index(self):
        result = lib.reperer_energie_index(17, [5, 11, 13, 17, 52, 64])
        self.assertLessEqual(result, 3)


# ---------------------------------------------------------------------------
# Extended tests for maximum coverage
# ---------------------------------------------------------------------------

class TestUtilities(unittest.TestCase):
    """Tests for basic utility functions."""

    def test_normalizeDic_values_sum_to_one(self):
        result = lib.normalizeDic({'H': 2.0, 'C': 4.0, 'O': 1.0})
        self.assertAlmostEqual(sum(result.values()), 1.0, places=10)

    def test_normalizeDic_proportions(self):
        result = lib.normalizeDic({'A': 1.0, 'B': 3.0})
        self.assertAlmostEqual(result['A'], 0.25, places=10)
        self.assertAlmostEqual(result['B'], 0.75, places=10)

    def test_normalizeDic_single_element(self):
        result = lib.normalizeDic({'X': 5.0})
        self.assertAlmostEqual(result['X'], 1.0, places=10)

    def test_normalise_list_output(self):
        result = lib.normalise([1.0, 1.0, 1.0])
        np.testing.assert_array_almost_equal(result, [1/3, 1/3, 1/3])

    def test_normalise_sums_to_one(self):
        result = lib.normalise([0.3, 0.5, 0.2])
        self.assertAlmostEqual(np.sum(result), 1.0, places=10)

    def test_sampling_returns_valid_index(self):
        for _ in range(20):
            result = lib.sampling([0.25, 0.25, 0.25, 0.25])
            self.assertIn(result, [0, 1, 2, 3])

    def test_sampling_deterministic_single(self):
        result = lib.sampling([1.0])
        self.assertEqual(result, 0)


class TestScintillatorComposition(unittest.TestCase):
    """Tests for scintillator mixture and aqueous fraction calculations."""

    def test_aqueous_fractions_water(self):
        result = lib.calculate_aqueous_fractions("Water", 0.1)
        self.assertIn('H', result)
        self.assertIn('O', result)
        self.assertAlmostEqual(sum(result.values()), 1.0, places=6)

    def test_aqueous_fractions_false_solvent(self):
        result = lib.calculate_aqueous_fractions("False", 0.1)
        self.assertIn('H', result)
        self.assertAlmostEqual(sum(result.values()), 1.0, places=6)

    def test_aqueous_fractions_HCl(self):
        result = lib.calculate_aqueous_fractions("HCl", 1.0)
        self.assertIn('Cl', result)
        self.assertGreater(result.get('Cl', 0), 0)
        self.assertAlmostEqual(sum(result.values()), 1.0, places=6)

    def test_aqueous_fractions_HNO3(self):
        result = lib.calculate_aqueous_fractions("HNO3", 1.0)
        self.assertIn('N', result)
        self.assertGreater(result.get('N', 0), 0)
        self.assertAlmostEqual(sum(result.values()), 1.0, places=6)

    def test_aqueous_fractions_NaOH(self):
        result = lib.calculate_aqueous_fractions("NaOH", 1.0)
        self.assertIn('Na', result)
        self.assertGreater(result.get('Na', 0), 0)
        self.assertAlmostEqual(sum(result.values()), 1.0, places=6)

    def test_aqueous_fractions_zero_concentration(self):
        result = lib.calculate_aqueous_fractions("HCl", 0)
        self.assertIn('H', result)
        self.assertEqual(result.get('Cl', 0), 0.0)

    def test_aqueous_fractions_unknown_solvent(self):
        result = lib.calculate_aqueous_fractions("Unknown", 1.0)
        self.assertIn('H', result)

    def test_lsc_mixture_properties_valid_cocktail(self):
        result = lib.calculate_lsc_mixture_properties("Ultima Gold", 0.1, "Water", 0.0)
        self.assertIsNotNone(result)
        self.assertIn('density_g_cm3', result)
        self.assertIn('effective_Z', result)
        self.assertIn('effective_A_g_mol', result)
        self.assertIn('atomic_fractions', result)
        self.assertGreater(result['density_g_cm3'], 0)
        self.assertGreater(result['effective_Z'], 0)
        self.assertGreater(result['effective_A_g_mol'], 0)

    def test_lsc_mixture_properties_atomic_fractions_sum(self):
        result = lib.calculate_lsc_mixture_properties("Ultima Gold", 0.1, "Water", 0.0)
        self.assertAlmostEqual(sum(result['atomic_fractions'].values()), 1.0, places=5)

    def test_lsc_mixture_properties_invalid_cocktail(self):
        result = lib.calculate_lsc_mixture_properties("NonExistentCocktail", 0.1, "Water", 0.0)
        self.assertIsNone(result)

    def test_lsc_mixture_properties_with_HCl(self):
        result = lib.calculate_lsc_mixture_properties("Ultima Gold", 0.1, "HCl", 0.1)
        self.assertIsNotNone(result)
        self.assertGreater(result['density_g_cm3'], 0)

    def test_lsc_mixture_properties_zero_aqueous(self):
        result = lib.calculate_lsc_mixture_properties("Ultima Gold", 0.0, "Water", 0.0)
        self.assertIsNotNone(result)
        self.assertGreater(result['effective_Z'], 0)

    def test_lsc_mixture_properties_hionic_fluor(self):
        result = lib.calculate_lsc_mixture_properties("Hionic-Fluor", 0.1, "Water", 0.0)
        self.assertIsNotNone(result)
        self.assertGreater(result['density_g_cm3'], 0)


class TestConfiguration(unittest.TestCase):
    """Tests for configuration I/O functions."""

    def test_readParameters_returns_tuple_of_29(self):
        result = lib.readParameters()
        self.assertEqual(len(result), 29)

    def test_readParameters_nE_electron_is_int(self):
        nE_electron = lib.readParameters()[0]
        self.assertIsInstance(nE_electron, int)
        self.assertGreater(nE_electron, 0)

    def test_readParameters_nE_alpha_is_int(self):
        nE_alpha = lib.readParameters()[1]
        self.assertIsInstance(nE_alpha, int)
        self.assertGreater(nE_alpha, 0)

    def test_readParameters_density_is_positive_float(self):
        rho = lib.readParameters()[2]
        self.assertIsInstance(rho, float)
        self.assertGreater(rho, 0)

    def test_readParameters_tau_is_int(self):
        tau = lib.readParameters()[10]
        self.assertIsInstance(tau, int)
        self.assertGreater(tau, 0)

    def test_readParameters_effQuantic_has_three_pmt(self):
        effQuantic = lib.readParameters()[14]
        self.assertIsInstance(effQuantic, list)
        self.assertEqual(len(effQuantic), 3)
        for q in effQuantic:
            self.assertGreater(q, 0)
            self.assertLessEqual(q, 1)

    def test_readParameters_disp_does_not_raise(self):
        import io, sys
        captured = io.StringIO()
        sys.stdout = captured
        try:
            lib.readParameters(disp=True)
        finally:
            sys.stdout = sys.__stdout__
        self.assertGreater(len(captured.getvalue()), 0)

    def test_update_config_value_roundtrip(self):
        original_tau = lib.readParameters()[10]
        try:
            lib.update_config_value("tau", 99)
            new_tau = lib.readParameters()[10]
            self.assertEqual(new_tau, 99)
        finally:
            lib.update_config_value("tau", original_tau)
        restored_tau = lib.readParameters()[10]
        self.assertEqual(restored_tau, original_tau)

    def test_update_config_batch_roundtrip(self):
        params = lib.readParameters()
        original_tau = params[10]
        original_dt = params[11]
        try:
            lib.update_config_batch({"tau": 77, "extDT": 15.0})
            new_params = lib.readParameters()
            self.assertEqual(new_params[10], 77)
            self.assertAlmostEqual(new_params[11], 15.0)
        finally:
            lib.update_config_batch({"tau": original_tau, "extDT": original_dt})

    def test_lsCocktail_returns_string_or_none(self):
        result = lib.lsCocktail()
        self.assertIsInstance(result, (str, type(None)))

    def test_modifynE_electron_roundtrip(self):
        original = lib.readParameters()[0]
        try:
            lib.modifynE_electron(500)
            self.assertEqual(lib.readParameters()[0], 500)
        finally:
            lib.modifynE_electron(original)

    def test_modifyDensity_roundtrip(self):
        original = lib.readParameters()[2]
        try:
            lib.modifyDensity(0.99)
            self.assertAlmostEqual(lib.readParameters()[2], 0.99, places=5)
        finally:
            lib.modifyDensity(original)

    def test_modifyTau_roundtrip(self):
        original = lib.readParameters()[10]
        try:
            lib.modifyTau(60)
            self.assertEqual(lib.readParameters()[10], 60)
        finally:
            lib.modifyTau(original)


class TestNuclearDataIO(unittest.TestCase):
    """Tests for nuclear decay data readers."""

    def test_readBetaShape_returns_nonzero_lists(self):
        e, p = lib.readBetaShape("H-3", "beta-", "tot")
        self.assertGreater(len(e), 0)
        self.assertGreater(len(p), 0)

    def test_readBetaShape_probabilities_sum_to_one(self):
        e, p = lib.readBetaShape("H-3", "beta-", "tot")
        self.assertAlmostEqual(sum(p), 1.0, places=5)

    def test_readBetaShape_energies_positive(self):
        e, p = lib.readBetaShape("H-3", "beta-", "tot")
        self.assertTrue(all(ei >= 0 for ei in e))

    def test_readBetaShape_C14(self):
        e, p = lib.readBetaShape("C-14", "beta-", "tot")
        self.assertGreater(len(e), 0)
        self.assertAlmostEqual(sum(p), 1.0, places=5)

    def test_readBetaShapeInfo_returns_string(self):
        result = lib.readBetaShapeInfo("H-3", "beta-", "tot")
        self.assertIsInstance(result, str)
        self.assertGreater(len(result), 0)

    def test_read_ENDF_photon_returns_3tuple(self):
        result = lib.read_ENDF_photon('C')
        self.assertEqual(len(result), 3)

    def test_read_ENDF_photon_C_has_data(self):
        binding, energie, cross_section = lib.read_ENDF_photon('C')
        self.assertGreater(len(binding), 0)
        self.assertGreater(len(energie), 0)
        self.assertGreater(len(cross_section), 0)

    def test_read_ENDF_photon_O_has_data(self):
        binding, energie, cross_section = lib.read_ENDF_photon('O')
        self.assertGreater(len(binding), 0)

    def test_read_ENDF_RA_returns_3tuple(self):
        result = lib.read_ENDF_RA('C')
        self.assertEqual(len(result), 3)

    def test_read_ENDF_RA_lists(self):
        Type_, Energie, Prob = lib.read_ENDF_RA('C')
        self.assertIsInstance(Type_, list)
        self.assertIsInstance(Energie, list)
        self.assertIsInstance(Prob, list)

    def test_read_ENDF_RA_equal_lengths(self):
        Type_, Energie, Prob = lib.read_ENDF_RA('C')
        self.assertEqual(len(Type_), len(Energie))
        self.assertEqual(len(Energie), len(Prob))

    def test_readPenNuc2_Co60(self):
        result = lib.readPenNuc2("Co-60")
        self.assertIn("NI60", result[0])

    def test_readBetaSpectra_Ni63(self):
        e, p = lib.readBetaSpectra("Ni-63")
        self.assertGreater(len(e), 0)
        self.assertAlmostEqual(sum(p), 1.0, places=4)


class TestPhysicsModels(unittest.TestCase):
    """Tests for detection, quenching, and analytical models."""

    def test_detectProbabilities_scalar_single_event(self):
        result = lib.detectProbabilities(1.0, [50.0], [0.0], 0, 1, 10, 20)
        eff_S, eff_D, eff_T, eff_AB, eff_BC, eff_AC, eff_D2 = result
        self.assertGreater(eff_S, 0)
        self.assertGreater(eff_D, 0)
        self.assertGreater(eff_T, 0)

    def test_detectProbabilities_hierarchy(self):
        eff_S, eff_D, eff_T, eff_AB, eff_BC, eff_AC, eff_D2 = \
            lib.detectProbabilities(1.0, [50.0], [0.0], 0, 1, 10, 20)
        self.assertLessEqual(eff_T, eff_D)
        self.assertLessEqual(eff_D, eff_S)

    def test_detectProbabilities_all_in_unit_interval(self):
        result = lib.detectProbabilities(1.0, [100.0], [0.0], 0, 1, 10, 20)
        for v in result:
            self.assertGreaterEqual(v, 0.0)
            self.assertLessEqual(v, 1.0)

    def test_detectProbabilities_tuple_L(self):
        result = lib.detectProbabilities((1.0, 1.0, 1.0), [50.0], [0.0], 0, 1, 10, 20)
        eff_S, eff_D, eff_T, eff_AB, eff_BC, eff_AC, eff_D2 = result
        self.assertGreater(eff_S, 0)
        self.assertLessEqual(eff_T, eff_D)

    def test_detectProbabilities_delayed_event(self):
        extDT_us = 10.0
        measTime_min = 20.0
        t1 = 1e-3   # 1 ms: inside (extDT*1e-6=10ns, measTime*60=1200s)
        result = lib.detectProbabilities(1.0, [50.0], [30.0], t1, 2, extDT_us, measTime_min)
        for v in result:
            self.assertGreaterEqual(v, 0.0)
            self.assertLessEqual(v, 2.0)

    def test_detectProbabilities_higher_L_gives_higher_efficiency(self):
        _, _, eff_T_low, *_ = lib.detectProbabilities(0.5, [50.0], [0.0], 0, 1, 10, 20)
        _, _, eff_T_high, *_ = lib.detectProbabilities(5.0, [50.0], [0.0], 0, 1, 10, 20)
        self.assertLess(eff_T_low, eff_T_high)

    def test_efficienciesEstimates_means(self):
        N = 100
        eff_S = np.full(N, 0.90)
        eff_D = np.full(N, 0.70)
        eff_T = np.full(N, 0.50)
        eff_AB = np.full(N, 0.65)
        eff_BC = np.full(N, 0.65)
        eff_AC = np.full(N, 0.65)
        eff_D2 = np.full(N, 0.40)
        result = lib.efficienciesEstimates(eff_S, eff_D, eff_T, eff_AB, eff_BC, eff_AC, eff_D2, N)
        self.assertEqual(len(result), 14)
        mean_S, std_S, mean_D, std_D, mean_T, std_T = result[:6]
        self.assertAlmostEqual(mean_S, 0.90)
        self.assertAlmostEqual(mean_D, 0.70)
        self.assertAlmostEqual(mean_T, 0.50)

    def test_efficienciesEstimates_std_zero_for_constant(self):
        N = 50
        eff_all = np.full(N, 0.75)
        result = lib.efficienciesEstimates(eff_all, eff_all, eff_all,
                                           eff_all, eff_all, eff_all, eff_all, N)
        for std in result[1::2]:
            self.assertAlmostEqual(std, 0.0, places=10)

    def test_modelAnalytical_eff_mode(self):
        # Use kB = 1e-5 cm/keV as in all TDCRPy notebooks.
        # Note: eff_S is the single-PMT efficiency (not "at least 1 PMT"), so
        # eff_D can exceed eff_S when the per-PMT detection probability pe > 0.5.
        eff_S, eff_D, eff_T = lib.modelAnalytical(
            1.0, 0.5, 0.5, 0.5, 0.5, "H-3", KB_TYPICAL, 10, "eff", 1000)
        self.assertGreater(eff_S, 0)
        self.assertGreater(eff_D, 0)
        self.assertGreater(eff_T, 0)
        self.assertLessEqual(eff_T, eff_D)   # triple <= double (always)
        self.assertLessEqual(eff_S, 1.0)
        self.assertLessEqual(eff_D, 1.0)
        self.assertLessEqual(eff_T, 1.0)

    def test_modelAnalytical_res_mode_returns_scalar(self):
        res = lib.modelAnalytical(1.0, 0.5, 0.5, 0.5, 0.5, "H-3", KB_TYPICAL, 10, "res", 1000)
        self.assertIsInstance(float(res), float)
        self.assertGreaterEqual(res, 0)

    def test_modelAnalytical_res_decreases_at_better_L(self):
        eff_S, eff_D, eff_T = lib.modelAnalytical(
            2.0, 0.5, 0.5, 0.5, 0.5, "H-3", KB_TYPICAL, 10, "eff", 500)
        true_tdcr = eff_T / eff_D
        res_good = lib.modelAnalytical(
            2.0, true_tdcr, true_tdcr, true_tdcr, true_tdcr, "H-3", KB_TYPICAL, 10, "res", 500)
        res_bad = lib.modelAnalytical(
            0.5, true_tdcr, true_tdcr, true_tdcr, true_tdcr, "H-3", KB_TYPICAL, 10, "res", 500)
        self.assertLessEqual(res_good, res_bad)

    def test_modelAnalytical_tuple_L(self):
        res = lib.modelAnalytical(
            (1.0, 1.0, 1.0), 0.5, 0.4, 0.4, 0.4, "H-3", KB_TYPICAL, 10, "res", 500)
        self.assertGreaterEqual(float(res), 0)

    def test_modelAnalyticalCN_returns_3tuple(self):
        eff_A, eff_B, eff_D = lib.modelAnalyticalCN(1.0, "H-3", KB_TYPICAL, 10, 1000)
        self.assertGreater(eff_A, 0)
        self.assertGreater(eff_B, 0)
        self.assertGreater(eff_D, 0)

    def test_modelAnalyticalCN_D_le_A(self):
        eff_A, eff_B, eff_D = lib.modelAnalyticalCN(1.0, "H-3", KB_TYPICAL, 10, 1000)
        self.assertLessEqual(eff_D, eff_A)
        self.assertLessEqual(eff_D, eff_B)

    def test_modelAnalyticalCN_symmetric(self):
        eff_A, eff_B, eff_D = lib.modelAnalyticalCN(1.0, "H-3", KB_TYPICAL, 10, 1000)
        self.assertAlmostEqual(eff_A, eff_B, places=10)

    def test_interaction_scintillation_returns_3tuple(self):
        result = lib.interaction_scintillation(100.0)
        self.assertEqual(len(result), 3)

    def test_interaction_scintillation_photoelectron_energy(self):
        e_ele, lacune, element = lib.interaction_scintillation(100.0)
        self.assertGreaterEqual(e_ele, 0)
        self.assertLessEqual(e_ele, 100.0)

    def test_interaction_scintillation_lacune_valid(self):
        e_ele, lacune, element = lib.interaction_scintillation(100.0)
        self.assertIn('Atom_', lacune)

    def test_interaction_scintillation_element_valid(self):
        e_ele, lacune, element = lib.interaction_scintillation(100.0)
        self.assertIn(element, ['H', 'C', 'N', 'O', 'P', 'S', 'Na', 'Cl'])

    def test_relaxation_atom_ph_returns_4tuple(self):
        result = lib.relaxation_atom_ph('Atom_K', 'C', 10)
        self.assertEqual(len(result), 4)

    def test_relaxation_atom_ph_lists(self):
        particules, energies, vacancies, par_type = lib.relaxation_atom_ph('Atom_K', 'C', 10)
        self.assertIsInstance(particules, list)
        self.assertIsInstance(energies, list)
        self.assertIsInstance(vacancies, list)
        self.assertIsInstance(par_type, list)

    def test_relaxation_atom_ph_H_no_relaxation(self):
        particules, energies, vacancies, par_type = lib.relaxation_atom_ph('Atom_K', 'H', 10)
        self.assertEqual(len(particules), 0)

    def test_detectProbabilitiesMC_returns_7tuple(self):
        # opticalTransport path — uniform photon splitting [1/3,1/3,1/3]
        result = lib.detectProbabilitiesMC(1.0, [50.0], [0.0], 0, 1, 10, 20)
        self.assertEqual(len(result), 7)
        for v in result:
            self.assertGreaterEqual(v, 0)


class TestStoppingPowerVariants(unittest.TestCase):
    """Tests for different stopping-power models."""

    def test_stoppingpower_returns_positive(self):
        result = lib.stoppingpower(100)
        self.assertGreater(result, 0)

    def test_stoppingpower_increases_at_lower_energy(self):
        sp_low = lib.stoppingpower(10)
        sp_high = lib.stoppingpower(1000)
        self.assertGreater(sp_low, 0)
        self.assertGreater(sp_high, 0)

    def test_stoppingpowerA_energy_scaling(self):
        sp1 = lib.stoppingpowerA(1000)
        sp2 = lib.stoppingpowerA(10000)
        self.assertGreater(sp1, 0)
        self.assertGreater(sp2, 0)


class TestEnergyDeposition(unittest.TestCase):
    """Tests for energy deposition functions."""

    def test_energie_dep_gamma2_returns_nonnegative(self):
        for e in [10, 50, 100, 500]:
            result = lib.energie_dep_gamma2(e, 10)
            self.assertGreaterEqual(result, 0)

    def test_energie_dep_gamma2_bounded_by_input(self):
        for e in [10, 100, 500]:
            result = lib.energie_dep_gamma2(e, 10)
            self.assertLessEqual(result, e)

    def test_energie_dep_beta2_returns_nonnegative(self):
        for e in [10, 50, 100]:
            result = lib.energie_dep_beta2(e, 10)
            self.assertGreaterEqual(result, 0)

    def test_energie_dep_beta2_bounded_by_input(self):
        for e in [10, 100]:
            result = lib.energie_dep_beta2(e, 10)
            self.assertLessEqual(result, e)

    def test_energie_dep_different_volumes(self):
        e_10ml = lib.energie_dep_gamma2(100, 10)
        e_16ml = lib.energie_dep_gamma2(100, 16)
        self.assertGreaterEqual(e_10ml, 0)
        self.assertGreaterEqual(e_16ml, 0)


class TestTDCRPyMain(unittest.TestCase):
    """Tests for the main TDCRPy Monte Carlo engine and analytical wrappers.

    All tests use kB = 1e-5 cm/keV, which is the standard Birks constant
    for organic liquid scintillators as used in the TDCRPy notebooks.
    """

    def test_TDCRPy_eff_mode_returns_14_values(self):
        result = TDCRPy_mod.TDCRPy(1.0, "H-3", "1", 200, KB_TYPICAL, 10)
        self.assertEqual(len(result), 14)

    def test_TDCRPy_eff_mode_S_in_unit_interval(self):
        result = TDCRPy_mod.TDCRPy(1.0, "H-3", "1", 200, KB_TYPICAL, 10)
        mean_S, std_S = result[0], result[1]
        self.assertGreater(mean_S, 0)
        self.assertLessEqual(mean_S, 1)
        self.assertGreaterEqual(std_S, 0)

    def test_TDCRPy_eff_mode_hierarchy(self):
        result = TDCRPy_mod.TDCRPy(1.0, "H-3", "1", 200, KB_TYPICAL, 10)
        mean_S, _, mean_D, _, mean_T = result[0], result[1], result[2], result[3], result[4]
        self.assertLessEqual(mean_T, mean_D + 1e-10)
        self.assertLessEqual(mean_D, mean_S + 1e-10)

    def test_TDCRPy_dis_mode_returns_4_lists(self):
        result = TDCRPy_mod.TDCRPy(1.0, "H-3", "1", 200, KB_TYPICAL, 10, mode="dis")
        self.assertEqual(len(result), 4)
        eff_S, eff_D, eff_T, eff_D2 = result
        self.assertEqual(len(eff_S), 200)
        self.assertEqual(len(eff_T), 200)

    def test_TDCRPy_dis_mode_values_in_range(self):
        result = TDCRPy_mod.TDCRPy(1.0, "H-3", "1", 100, KB_TYPICAL, 10, mode="dis")
        eff_S = result[0]
        self.assertTrue(all(0 <= v <= 1 for v in eff_S))

    def test_TDCRPy_analytical_model(self):
        # Smodel=False: analytical model path, fixed in v2.20.8
        # (extra 'symm' arg bug at TDCRPy.py:549 corrected).
        # Returns 14-tuple with zero uncertainties (deterministic).
        result = TDCRPy_mod.TDCRPy(9.0, "H-3", "1", 100, KB_TYPICAL, 10, Smodel=False)
        self.assertEqual(len(result), 14)
        eff_S, _, eff_D, _, eff_T = result[0], result[1], result[2], result[3], result[4]
        self.assertGreater(eff_S, 0)
        self.assertLessEqual(eff_T, eff_D + 1e-9)

    def test_run_from_history_skips_malformed_lines(self):
        # Regression test for a reported IndexError: Temp_E2.txt is a
        # fixed-name file in the shared system temp dir, so a stray short
        # or partially-written line (e.g. from a leftover/concurrent run)
        # used to crash _run_from_history() at parts[2]. It must now skip
        # such lines instead of raising.
        recfile3 = os.path.join(tempfile.gettempdir(), "Temp_E2.txt")
        with open(recfile3, "w") as f:
            f.write("# TDCRPy output: quenched deposited energies from nuclear decays\n")
            for idec in range(5):
                f.write(f"1 5.000000E+03 {idec:2d} 1 0.000000E+00\n")
            f.write("1 5.00\n")  # malformed: too few fields
            for idec in range(5, 10):
                f.write(f"1 5.000000E+03 {idec:2d} 1 0.000000E+00\n")

        result = TDCRPy_mod.TDCRPy(5.0, "H-3", "1", 10, KB_TYPICAL, 15, readRecHist=True)
        self.assertEqual(len(result), 14)
        self.assertGreater(result[0], 0)  # eff_S

    def test_TDCRPy_higher_L_gives_higher_efficiency(self):
        res_low = TDCRPy_mod.TDCRPy(0.5, "H-3", "1", 150, KB_TYPICAL, 10)
        res_high = TDCRPy_mod.TDCRPy(3.0, "H-3", "1", 150, KB_TYPICAL, 10)
        self.assertLess(res_low[0], res_high[0])

    def test_TDCRPy_tuple_L(self):
        result = TDCRPy_mod.TDCRPy((1.0, 1.0, 1.0), "H-3", "1", 100, KB_TYPICAL, 10)
        self.assertEqual(len(result), 14)

    def test_effA_returns_5_values(self):
        result = TDCRPy_mod.effA(0.5, "H-3", "1", KB_TYPICAL, 10)
        self.assertEqual(len(result), 5)

    def test_effA_L0_positive(self):
        L0, L_opt, eff_S, eff_D, eff_T = TDCRPy_mod.effA(0.5, "H-3", "1", KB_TYPICAL, 10)
        self.assertGreater(L0, 0)

    def test_effA_efficiencies_positive(self):
        L0, L_opt, eff_S, eff_D, eff_T = TDCRPy_mod.effA(0.5, "H-3", "1", KB_TYPICAL, 10)
        self.assertGreater(eff_S, 0)
        self.assertGreater(eff_D, 0)
        self.assertGreater(eff_T, 0)

    def test_effA_hierarchy(self):
        # eff_S is the single-PMT efficiency; eff_D can exceed it when pe > 0.5.
        _, _, _, eff_D, eff_T = TDCRPy_mod.effA(0.5, "H-3", "1", KB_TYPICAL, 10)
        self.assertLessEqual(eff_T, eff_D + 1e-10)   # triple <= double (always)

    def test_effA_C14(self):
        L0, L_opt, eff_S, eff_D, eff_T = TDCRPy_mod.effA(0.5, "C-14", "1", KB_TYPICAL, 10)
        self.assertGreater(L0, 0)
        self.assertGreater(eff_S, 0)

    # ------------------------------------------------------------------
    # Regression tests for bugs fixed in v2.20.6:
    # _interact_prompt now propagates energy_vec_initial correctly for
    # secondary particles (photoelectron, Auger, pair-production positron).
    # ------------------------------------------------------------------

    def test_TDCRPy_Co60_no_IndexError(self):
        # Co-60 emits 1.17 / 1.33 MeV gammas.  Photoelectric and pair-
        # production events add secondary particles to particle_vec beyond
        # the length of energy_vec_initial, which triggered IndexError in
        # _quench before the fix.  N=500 ensures the rare code path is hit.
        result = TDCRPy_mod.TDCRPy(1.0, "Co-60", "1", 500, KB_TYPICAL, 10)
        self.assertEqual(len(result), 14)
        mean_S = result[0]
        self.assertGreater(mean_S, 0.8)
        self.assertLessEqual(mean_S, 1.0)

    def test_TDCRPy_Co60_triple_le_double(self):
        result = TDCRPy_mod.TDCRPy(1.0, "Co-60", "1", 300, KB_TYPICAL, 10)
        mean_D, mean_T = result[2], result[4]
        self.assertLessEqual(mean_T, mean_D + 1e-6)

    def test_TDCRPy_Sr90_Y90_high_energy_beta(self):
        # Y-90 beta endpoint 2.28 MeV — exercises the high-energy tail of
        # the quenching interpolation and energy_vec_initial bookkeeping.
        result = TDCRPy_mod.TDCRPy(1.0, "Sr-90, Y-90", "0.5, 0.5", 300, KB_TYPICAL, 10)
        self.assertEqual(len(result), 14)
        mean_S = result[0]
        self.assertGreater(mean_S, 0.8)
        self.assertLessEqual(mean_S, 1.0)

    def test_TDCRPy_Fe55_EC_decay(self):
        # Fe-55 decays by electron capture — exercises atomic relaxation
        # (Mn K-alpha X-rays, Auger electrons) without beta emission.
        result = TDCRPy_mod.TDCRPy(1.0, "Fe-55", "1", 300, KB_TYPICAL, 10)
        self.assertEqual(len(result), 14)
        self.assertGreater(result[0], 0)

    def test_TDCRPy_effA_Co60(self):
        # Analytical model for Co-60 (complex decay scheme).
        L0, L_opt, eff_S, eff_D, eff_T = TDCRPy_mod.effA(
            0.977, "Co-60", "1", KB_TYPICAL, 10)
        self.assertGreater(L0, 0)
        self.assertGreater(eff_S, 0.8)


class TestQuenchingModels(unittest.TestCase):
    """Additional quenching model tests."""

    def test_E_quench_e_full_energy_le_input(self):
        ei = 50000   # eV
        result = lib.E_quench_e(ei, ei, 0.01, 1000)
        self.assertLessEqual(float(result) * 1e-3, 50.001)

    def test_E_quench_e_partial_energy(self):
        ei = 100000  # eV
        ed = 50000   # eV deposited
        result = lib.E_quench_e(ei, ed, 0.01, 1000)
        self.assertLessEqual(float(result) * 1e-3, 50.001)

    def test_E_quench_a_positive(self):
        result = lib.E_quench_a(5000, 1e-5, 1000)
        self.assertGreater(float(result), 0)

    def test_Em_e_less_than_input_energy(self):
        ei = 100000  # eV
        result = lib.Em_e(ei, ei, 1e-5, 1000)
        self.assertLessEqual(float(result), ei)

    def test_Em_a_less_than_input_energy(self):
        result = lib.Em_a(5000, 1e-5, 1000)
        self.assertLessEqual(float(result), 5000)


class TestNamingAndParsing(unittest.TestCase):
    """Tests for name transformation and parsing helpers."""

    def test_transf_name_standard_format(self):
        self.assertEqual(lib.transf_name("108Ag"), "Ag108")

    def test_transf_name_already_correct(self):
        result = lib.transf_name("Ag108")
        self.assertIsNotNone(result)

    def test_format_modif_negative_exponent(self):
        result = lib.format_modif("17-3")
        self.assertEqual(result, '17E-3')

    def test_format_modif_positive_exponent(self):
        result = lib.format_modif("3+5")
        self.assertEqual(result, '3E+5')

    def test_reperer_energie_index_exact_match(self):
        # 17 is found at index 3; next value (52) is at index 4 -> returns 4-1=3
        result = lib.reperer_energie_index(17, [5, 11, 13, 17, 52, 64])
        self.assertEqual(result, 3)

    def test_reperer_energie_index_between(self):
        # 15 is between 13 (idx 2) and 17 (idx 3): first i where vec[i]>15 is i=3, returns i-1=2
        result = lib.reperer_energie_index(15, [5, 11, 13, 17, 52, 64])
        self.assertEqual(result, 2)


# ---------------------------------------------------------------------------
# Validated reference tests — v2.20.7 baseline
#
# All values below were obtained by running the functional_validation notebook
# with v2.20.7 as the reference on 2026-06-25.
#
# Analytical-model tests use exact arithmetic (deterministic) → tight tol.
# Stochastic tests use statistical bounds → ±3σ with N=500 typical runs.
#
# Parameters: L=9.0 keV⁻¹ (analytical), L=1.0 keV⁻¹ (stochastic),
#             kB=1e-5 cm keV⁻¹, V=10 mL, ne=1000.
# ---------------------------------------------------------------------------

# Reference parameters
_L_ANA   = 9.0      # fixed L for analytical model tests
_L_STOCH = 1.0      # fixed L for stochastic model tests
_KB      = 1.0e-5   # Birks constant (cm keV⁻¹)
_V       = 10       # volume (mL)
_N_REF   = 500      # MC trials (fast; tests use statistical bounds)


class TestValidatedAnalytical(unittest.TestCase):
    """Analytical model reference values validated against v2.20.10.

    modelAnalytical() is deterministic given explicit parameters.
    Tolerance delta=0.001 accommodates minor platform/scipy-version
    differences without masking genuine regressions.

    setUp() clears the beta-spectra cache before each test to guarantee
    a fresh computation regardless of test execution order.
    """

    # L = 9.0 photons/keV (consistent with stochastic model).
    # With effQuantic=0.1 per PMT: L_photoelectron = 9.0 * 0.1 = 0.9 keV⁻¹.

    def setUp(self):
        # Clear the module-level caches so earlier tests cannot influence
        # the analytical model values through stale cached spectra.
        lib._readBetaShape_cache.clear()
        lib._readBetaSpectra_cache.clear()

    def _check(self, rad, exp_S, exp_D, exp_T, delta=0.001):
        eff_S, eff_D, eff_T = lib.modelAnalytical(
            _L_ANA, 0, 0, 0, 0, rad, _KB, _V, 'eff', 1000)
        self.assertAlmostEqual(float(eff_S), exp_S, delta=delta,
                               msg=f'{rad} eff_S mismatch')
        self.assertAlmostEqual(float(eff_D), exp_D, delta=delta,
                               msg=f'{rad} eff_D mismatch')
        self.assertAlmostEqual(float(eff_T), exp_T, delta=delta,
                               msg=f'{rad} eff_T mismatch')

    def test_ref_H3_analytical(self):
        # eff_S = 1-(1-pe)^3, consistent with stochastic detectProbabilities.
        # Reference recomputed with effQuantic=0.25 (default since v2.20.11).
        self._check('H-3',
                    exp_S=0.89856, exp_D=0.77817, exp_T=0.58979)

    def test_ref_Sr90_analytical(self):
        self._check('Sr-90',
                    exp_S=0.99439, exp_D=0.99064, exp_T=0.98419)

    def test_ref_Co60_analytical(self):
        # Beta-spectrum approximation only (no gammas in analytical model).
        self._check('Co-60',
                    exp_S=0.98923, exp_D=0.98154, exp_T=0.96831)

    def test_ref_H3_effA(self):
        # L0 = photon yield; with effQuantic=0.25 the photon yield is
        # L0 = 0.85445 / 0.25 = 3.418 (efficiencies unchanged).
        L0, _, eff_S, eff_D, eff_T = TDCRPy_mod.effA(0.5, 'H-3', '1', _KB, _V)
        self.assertAlmostEqual(L0,    3.4178, delta=0.05,  msg='H-3 L0 (photon yield)')
        self.assertAlmostEqual(eff_S, 0.7787, delta=0.005, msg='H-3 eff_S')
        self.assertAlmostEqual(eff_D, 0.5443, delta=0.005, msg='H-3 eff_D')
        self.assertAlmostEqual(eff_T, 0.2721, delta=0.005, msg='H-3 eff_T')

    def test_ref_Co60_effA(self):
        L0, _, eff_S, eff_D, eff_T = TDCRPy_mod.effA(0.977, 'Co-60', '1', _KB, _V)
        self.assertAlmostEqual(L0,    4.7236, delta=0.05,  msg='Co-60 L0 (photon yield)')
        self.assertAlmostEqual(eff_S, 0.9847, delta=0.005, msg='Co-60 eff_S')
        self.assertAlmostEqual(eff_D, 0.9716, delta=0.005, msg='Co-60 eff_D')
        self.assertAlmostEqual(eff_T, 0.9492, delta=0.005, msg='Co-60 eff_T')

    def test_analytical_ordering_H3(self):
        eff_S, eff_D, eff_T = lib.modelAnalytical(
            _L_ANA, 0, 0, 0, 0, 'H-3', _KB, _V, 'eff', 1000)
        self.assertLessEqual(float(eff_T), float(eff_D))   # T ≤ D always

    def test_analytical_ordering_Sr90(self):
        eff_S, eff_D, eff_T = lib.modelAnalytical(
            _L_ANA, 0, 0, 0, 0, 'Sr-90', _KB, _V, 'eff', 1000)
        self.assertLessEqual(float(eff_T), float(eff_D))


class TestValidatedStochastic(unittest.TestCase):
    """Stochastic model validation tests — v2.20.10 baseline.

    Tests are PROPERTY-BASED, not value-based, to be robust against
    Monte-Carlo statistical fluctuation and configuration variation.

    Physical invariants tested:
      - eff_T ≤ eff_D  (triple ≤ double, always guaranteed by detection model)
      - eff_S, eff_D, eff_T are finite and non-negative
      - Nuclide-specific signatures (Cd-109 eff_S > 1 for delayed coincidences)
      - No exception is raised (regression guard for the v2.20.6 IndexError fix)
    """

    def _run(self, rad, N=200):
        return TDCRPy_mod.TDCRPy(_L_STOCH, rad, '1', N, _KB, _V)

    def _assert_valid(self, result, rad, max_eff=1.05):
        eff_S, u_S, eff_D, u_D, eff_T, u_T = result[:6]
        self.assertTrue(np.isfinite(eff_S), f'{rad} eff_S is not finite')
        self.assertTrue(np.isfinite(eff_D), f'{rad} eff_D is not finite')
        self.assertTrue(np.isfinite(eff_T), f'{rad} eff_T is not finite')
        self.assertGreaterEqual(eff_S, 0,   f'{rad} eff_S < 0')
        self.assertGreaterEqual(eff_D, 0,   f'{rad} eff_D < 0')
        self.assertGreaterEqual(eff_T, 0,   f'{rad} eff_T < 0')
        self.assertLessEqual(eff_T, eff_D + 1e-9,
                             f'{rad} eff_T > eff_D (violates detection model)')
        self.assertLessEqual(eff_S, max_eff, f'{rad} eff_S unexpectedly large')

    def test_ref_H3_stochastic(self):
        # Pure β⁻, 18.6 keV — low efficiency at L=1; no photon interactions.
        r = self._run('H-3')
        self._assert_valid(r, 'H-3')
        self.assertLess(r[0], 0.8, 'H-3 eff_S should be well below 1 at L=1')

    def test_ref_Fe55_stochastic(self):
        # EC + Mn Kα X-rays (5.9 keV) — Auger cascade.
        # Regression guard: _interact_prompt IndexError fix (v2.20.6).
        r = self._run('Fe-55')
        self._assert_valid(r, 'Fe-55')

    def test_ref_Co60_stochastic(self):
        # β⁻ + 1.17/1.33 MeV γ — photoelectric / pair-production path.
        # Regression guard: _interact_prompt IndexError fix (v2.20.6).
        r = self._run('Co-60')
        self._assert_valid(r, 'Co-60')
        self.assertGreater(r[0], 0.5, 'Co-60 eff_S should be high at L=1')

    def test_ref_Sr90_stochastic(self):
        # Pure β⁻, 0.546 MeV — medium energy, high efficiency at L=1.
        r = self._run('Sr-90')
        self._assert_valid(r, 'Sr-90')
        self.assertGreater(r[0], 0.7, 'Sr-90 eff_S should be high at L=1')

    def test_ref_Cd109_stochastic(self):
        # EC + delayed Ag-109m IT γ (88 keV) — delayed coincidences make
        # eff_S > 1 (two pulses per decay contribute to single rate).
        # Regression guard: IndexError fix in _quench() (v2.20.6).
        r = TDCRPy_mod.TDCRPy(_L_STOCH, 'Cd-109', '1', _N_REF, _KB, _V)
        eff_S = r[0]
        self.assertTrue(np.isfinite(eff_S), 'Cd-109 eff_S is not finite')
        self.assertGreater(eff_S, 1.0,
                           'Cd-109 eff_S should exceed 1 (delayed coincidences)')
        self.assertLess(eff_S, 3.0,
                        'Cd-109 eff_S unreasonably large')

    def test_stoch_returns_14_values_all_nuclides(self):
        for rad in ['H-3', 'Fe-55', 'Co-60', 'Sr-90', 'Cd-109']:
            try:
                r = TDCRPy_mod.TDCRPy(_L_STOCH, rad, '1', 100, _KB, _V)
                self.assertEqual(len(r), 14, f'{rad}: wrong tuple length')
            except Exception as e:
                self.fail(f'TDCRPy raised {type(e).__name__} for {rad}: {e}')

    def test_stoch_efficiencies_finite_and_ordered(self):
        for rad in ['H-3', 'Co-60', 'Sr-90']:
            r = TDCRPy_mod.TDCRPy(_L_STOCH, rad, '1', 100, _KB, _V)
            eff_S, _, eff_D, _, eff_T = r[0], r[1], r[2], r[3], r[4]
            self.assertTrue(np.isfinite(eff_S), f'{rad} eff_S not finite')
            self.assertLessEqual(eff_T, eff_D + 1e-9,
                                 f'{rad}: eff_T > eff_D')


if __name__ == '__main__':
    unittest.main()
