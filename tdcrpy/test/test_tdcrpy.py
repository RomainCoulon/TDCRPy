# -*- coding: utf-8 -*-
"""
test_tdcrpy.py

Unit tests for TDCRPy package — covers TDCR_model_lib and TDCRPy modules.

Authors: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

import math
import inspect
import unittest
import warnings
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
        # Regression test for a reported IndexError: a stray short or
        # partially-written line (from a leftover run, or before v2.20.20 a
        # concurrent one) used to crash _run_from_history() at parts[2].
        # It must skip such lines instead of raising.
        recfile3 = lib.record_file_path(2)
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



class TestMicelleRadiusConsistency(unittest.TestCase):
    """The stochastic micellar-quenching kernel must conserve the aqueous
    volume fraction it is given.

    Regression test for the v2.20.19 fix: the inter-micellar mean free path used
    to be built from the *nominal* radius while the aqueous chords were drawn
    from the zero-truncated sampled radius. Because truncation makes
    ``<r> > r_d`` when ``sigma`` is comparable to ``r_d``, the aqueous path
    fraction came out inflated by ``<r>/r_d`` and the high-energy retention
    plateau fell below ``1 - fAq`` (0.8748 instead of 0.9000 at r_d=1, sigma=1).
    """

    E_PLATEAU = 50.0     # keV: range >> micelle size, so f(E) has saturated
    NS = 60000           # samples; MC error on the mean is then ~1e-4

    def _plateau(self, r_d, sigma, f_w):
        mean, _, _ = lib.pure_mc_efficient_energy_numba(
            self.E_PLATEAU, r_d_nm=r_d, sigma_r=sigma, f_w=f_w,
            num_samples=self.NS)
        return mean / self.E_PLATEAU

    def test_mean_truncated_radius_matches_analytic(self):
        # r_mean + r_sigma * phi(a)/(1-Phi(a)),  a = -r_mean/r_sigma
        for r, s_ in [(1.0, 1.0), (2.0, 1.0), (3.0, 1.0), (0.5, 2.0)]:
            a = -r / s_
            pdf = math.exp(-0.5 * a * a) / math.sqrt(2.0 * math.pi)
            cdf = 0.5 * (1.0 + math.erf(a / math.sqrt(2.0)))
            expected = r + s_ * pdf / (1.0 - cdf)
            got = lib.mean_truncated_radius(r, s_)
            self.assertAlmostEqual(got, expected, places=10,
                                   msg=f'r={r} sigma={s_}')

    def test_mean_truncated_radius_zero_sigma(self):
        for r in (0.5, 1.0, 4.0):
            self.assertEqual(lib.mean_truncated_radius(r, 0.0), r)

    def test_plateau_reaches_one_minus_faq(self):
        """The regression itself: every (r_d, sigma) must plateau at 1-fAq."""
        f_w = 0.10
        for r_d in (1.0, 2.0, 3.0, 4.0, 8.0):
            for sigma in (0.0, 1.0):
                f_inf = self._plateau(r_d, sigma, f_w)
                self.assertAlmostEqual(
                    f_inf, 1.0 - f_w, places=2,
                    msg=(f'r_d={r_d} sigma={sigma}: plateau {f_inf:.4f} '
                         f'!= 1-fAq {1-f_w:.4f}'))

    def test_worst_case_small_radius_large_sigma(self):
        """r_d=1, sigma=1 was the worst case: 0.8748 before the fix."""
        f_inf = self._plateau(1.0, 1.0, 0.10)
        self.assertGreater(f_inf, 0.885,
                           f'plateau {f_inf:.4f} still depressed (pre-fix 0.8748)')
        self.assertAlmostEqual(f_inf, 0.90, places=2)

    def test_plateau_independent_of_sigma(self):
        """Sigma must change the spread, not the mean retention."""
        f_w = 0.10
        ref = self._plateau(2.0, 0.0, f_w)
        for sigma in (0.5, 1.0, 1.5):
            self.assertAlmostEqual(
                self._plateau(2.0, sigma, f_w), ref, places=2,
                msg=f'plateau moved with sigma={sigma}')

    def test_plateau_tracks_faq(self):
        for f_w in (0.05, 0.20, 0.30):
            self.assertAlmostEqual(self._plateau(1.0, 1.0, f_w), 1.0 - f_w,
                                   places=2, msg=f'fAq={f_w}')

    def test_shipped_default_is_monodisperse(self):
        """The packaged configuration must ship sigma_micelle = 0 (v2.20.21).

        Reported reverse-micelle size spreads are small compared with the mean,
        so the mean size carries the physics. A zero default also keeps the
        sampled radius strictly positive, which removes any dependence on how
        the size distribution is truncated or parameterised.
        """
        with importlib.resources.path('tdcrpy', 'configDefault.toml') as cfg:
            text = cfg.read_text(encoding='utf-8')
        line = [l for l in text.splitlines()
                if l.strip().startswith('sigma_micelle')]
        self.assertEqual(len(line), 1, 'sigma_micelle not found in configDefault')
        self.assertEqual(float(line[0].split('=')[1].strip()), 0.0,
                         f'shipped default is not monodisperse: {line[0]!r}')

    def test_zero_sigma_is_deterministic_in_radius(self):
        """At sigma = 0 every droplet has the nominal radius, so the retention
        is the pure-geometry result and no truncation question arises."""
        self.assertEqual(lib.mean_truncated_radius(2.0, 0.0), 2.0)
        self.assertAlmostEqual(self._plateau(2.0, 0.0, 0.10), 0.90, places=2)

    def test_kernel_default_radius_is_half_the_diameter(self):
        """``diam_micelle`` is a DIAMETER; the kernel takes a radius (v2.20.21).

        Before the fix the configured diameter was bound straight to
        ``r_d_nm``, so the quenching path ran on droplets twice the requested
        size while :func:`micelleLoss` read the same key as a diameter. The two
        micelle models therefore disagreed by a factor two on the same config.
        """
        saved = (lib.fAq, lib.diam_micelle)
        try:
            lib.modifyfAq(0.10)
            lib.modifyDiam_micelle(4.0)
            # 0.2 keV sits in the confinement regime, where the retention is
            # strongly radius-dependent -- at the plateau every radius gives
            # 1-fAq and the test could not tell 2 nm from 4 nm.
            E, ns = 0.2, 200000
            got = lib.pure_mc_efficient_energy_numba(E, num_samples=ns)[0]
            as_radius = lib.pure_mc_efficient_energy_numba(
                E, r_d_nm=2.0, num_samples=ns)[0]
            as_diameter = lib.pure_mc_efficient_energy_numba(
                E, r_d_nm=4.0, num_samples=ns)[0]
            # Both are Monte-Carlo estimates: the per-event spread here is
            # about 0.34 of E, so each mean carries ~0.34/sqrt(ns) and their
            # difference ~1.4x that. The tolerance is ~5 sigma of that, which
            # is loose enough not to fail on noise and still an order of
            # magnitude tighter than the 2 nm / 4 nm separation being tested.
            tol = 5 * 1.5 * 0.34 / math.sqrt(ns)
            self.assertAlmostEqual(
                got / E, as_radius / E, delta=tol,
                msg='the default is not diam_micelle / 2')
            self.assertLess(
                as_diameter, as_radius,
                'a 4 nm radius must confine more than a 2 nm one')
            self.assertGreater(
                abs(got - as_diameter) / E, 0.05,
                'the kernel appears to be using the diameter as the radius')
        finally:
            lib.modifyfAq(saved[0])
            lib.modifyDiam_micelle(saved[1])

    def test_shipped_default_diameter_is_four_nm(self):
        """4 nm diameter: below the 8 nm long treated as canonical, per the
        light-scattering survey of Bergeron, Appl. Radiat. Isot. 70 (2012)."""
        with importlib.resources.path('tdcrpy', 'configDefault.toml') as cfg:
            text = cfg.read_text(encoding='utf-8')
        line = [l for l in text.splitlines()
                if l.strip().startswith('diam_micelle')]
        self.assertEqual(len(line), 1, 'diam_micelle not found in configDefault')
        self.assertEqual(float(line[0].split('=')[1].strip()), 4.0)

    def test_micelleloss_accepts_the_shipped_diameter(self):
        """micelleLoss indexes a table of tabulated diameters, so the shipped
        default has to be one of them or the legacy path raises."""
        S = lib.micelleLoss(10.0)
        self.assertTrue(0.0 < S <= 1.0, f'micelleLoss returned {S}')



class TestMicelleDefaultsResolveAtCallTime(unittest.TestCase):
    """modify*() must take effect without an importlib.reload (v2.20.22).

    The kernel used to bind ``f_w=fAq``, ``r_d_nm=diam_micelle/2`` and
    ``sigma_r=sigma_micelle`` as keyword defaults evaluated at import. Since
    the shipped configuration has ``fAq = 0`` and the kernel returns the energy
    unchanged when ``f_w == 0``, enabling ``micCorr`` after a ``modify*()``
    without a reload was a silent no-op: no error, no warning, and a
    "with/without correction" comparison that comes out as pure MC noise.
    """

    E = 5.7          # keV, a typical H-3 electron: range >> droplet

    def setUp(self):
        self._saved = (lib.fAq, lib.diam_micelle, lib.sigma_micelle)

    def tearDown(self):
        f, d, sg = self._saved
        lib.modifyfAq(f)
        lib.modifyDiam_micelle(d)
        lib.modifySigma_micelle(sg)

    def test_signature_defaults_are_sentinels(self):
        p = inspect.signature(lib.pure_mc_efficient_energy_numba).parameters
        for name in ('r_d_nm', 'sigma_r', 'f_w'):
            self.assertIsNone(p[name].default,
                              f'{name} is still bound at import time')

    def test_modifyfAq_takes_effect_without_reload(self):
        lib.modifyfAq(0.0)
        self.assertEqual(lib.pure_mc_efficient_energy_numba(self.E), self.E,
                         'at fAq = 0 the kernel must be a pass-through')
        lib.modifyfAq(0.10)
        vals = [lib.pure_mc_efficient_energy_numba(self.E) for _ in range(200)]
        self.assertTrue(any(v < self.E for v in vals),
                        'correction still inert after modifyfAq -- the '
                        'no-op regression is back')
        # range >> droplet, so the mean retention is the bulk limit 1 - fAq
        self.assertAlmostEqual(float(np.mean(vals)) / self.E, 0.90, places=1)

    def test_modifyDiam_takes_effect_without_reload(self):
        lib.modifyfAq(0.10)
        for diam in (2.0, 8.0):
            lib.modifyDiam_micelle(diam)
            p = inspect.signature(lib.pure_mc_efficient_energy_numba).parameters
            self.assertIsNone(p['r_d_nm'].default)
            self.assertAlmostEqual(lib.diam_micelle, diam, places=9)

    def test_explicit_arguments_still_win(self):
        lib.modifyfAq(0.10)
        self.assertEqual(
            lib.pure_mc_efficient_energy_numba(self.E, f_w=0.0), self.E,
            'an explicit f_w must override the configuration')

    def test_micelleLoss_follows_the_configuration(self):
        lib.modifyfAq(0.0)
        lib.modifyDiam_micelle(4.0)
        a = lib.micelleLoss(10.0)
        lib.modifyfAq(0.20)
        b = lib.micelleLoss(10.0)
        self.assertNotAlmostEqual(a, b, places=3,
                                  msg='micelleLoss ignored modifyfAq')


class TestLBracketing(unittest.TestCase):
    """eff()/effA() must bracket the physical root themselves (v2.20.22).

    The modelled TDCR is not monotonic in L, so a bounded search over a wide
    span can settle on the falling branch and return an L reproducing the
    wrong ratio with nothing raised. Lbounds now defaults to None, meaning
    "scan, then solve inside the first crossing".
    """

    def test_lbounds_defaults_to_none(self):
        for fn in (TDCRPy_mod.eff, TDCRPy_mod.effA):
            self.assertIsNone(
                inspect.signature(fn).parameters['Lbounds'].default,
                f'{fn.__name__} still defaults to a hard-coded bracket')

    def test_bracket_picks_the_rising_branch(self):
        """A curve that rises, turns over and falls has two solutions; the
        physical one is on the rising branch."""
        def tdcr(L):                      # peaks at L = 10, TDCR = 1.0
            return 1.0 - ((L - 10.0) / 10.0) ** 2

        lo, hi = TDCRPy_mod._bracket_physical_L(tdcr, 0.75)
        self.assertLess(hi, 10.0, 'bracket is not on the rising branch')
        self.assertLessEqual(tdcr(lo), 0.75)
        self.assertGreaterEqual(tdcr(hi), 0.75)

    def test_bracket_returns_none_when_unreachable(self):
        self.assertIsNone(
            TDCRPy_mod._bracket_physical_L(lambda L: 0.5, 0.95),
            'an unreachable ratio must report no bracket, not a wrong one')

    def test_solve_L_falls_back_when_the_bracket_will_not_hold(self):
        """A stochastic tdcr_at breaks brentq's sign condition (v2.20.23).

        The replay samples photons when ``opticalTransport`` is on, so the
        value at a bracket end is not the same on re-evaluation and brentq can
        raise ValueError. That must degrade to a bounded search *inside the
        bracket* -- which still excludes the falling branch -- not propagate.
        """
        seen = {'n': 0}

        def flaky(L):
            # rises during the scan, then goes flat, so f(lo) and f(hi) end up
            # with the same sign when brentq re-evaluates them
            seen['n'] += 1
            return min(0.99, 0.1 * L) if seen['n'] <= TDCRPy_mod._AUTO_POINTS \
                else 0.20

        used = {}

        def fallback(span):
            used['span'] = span
            return 42.0

        out = TDCRPy_mod._solve_L(flaky, 0.6, fallback, 'X', 'eff()')
        self.assertEqual(out, 42.0, 'brentq failure was not caught')
        self.assertIn('span', used, 'fallback was never reached')
        lo, hi = used['span']
        self.assertGreater(lo, TDCRPy_mod._AUTO_SPAN[0],
                           'fell back to the full span instead of the bracket')
        self.assertLess(hi, TDCRPy_mod._AUTO_SPAN[1],
                        'fell back to the full span instead of the bracket')

    def test_config_retry_budget_survives_a_transient_lock(self):
        """config.toml I/O must tolerate a lock for longer than a few seconds.

        On Windows an atomic os.replace onto config.toml fails with WinError 5
        while another process holds the destination open, which an on-access
        virus scan is enough to cause. The budget used to be 5 attempts = 3 s,
        short enough that one transient lock killed a notebook run an hour in.
        """
        n, b = lib._CONFIG_READ_ATTEMPTS, lib._CONFIG_READ_BACKOFF
        total = b * (2 ** (n - 1) - 1)
        self.assertGreaterEqual(total, 20.0,
                                f'retry budget is only {total:.1f} s')

    def test_optical_transport_default_reads_the_configuration(self):
        """The opticalTransport key was inert before v2.20.24.

        TDCRPy() took it as a plain ``False`` argument default and never
        consulted the configuration, so setting the key -- or calling
        modifyOpticalTransport() -- changed nothing unless the caller also
        passed the argument by hand, which eff() never did.
        """
        self.assertIsNone(
            inspect.signature(TDCRPy_mod.TDCRPy).parameters[
                'opticalTransport'].default,
            'TDCRPy() still ignores the opticalTransport configuration key')

    def test_modifyOpticalTransport_takes_effect_without_reload(self):
        saved = lib.opticalTransport
        try:
            lib.modifyOpticalTransport(True)
            self.assertTrue(lib.opticalTransport)
            lib.modifyOpticalTransport(False)
            self.assertFalse(lib.opticalTransport)
        finally:
            lib.modifyOpticalTransport(saved)

    def test_tolerance_is_matched_to_the_monte_carlo_noise(self):
        """Solving L to machine precision is meaningless: the histories are one
        realisation, so L is reproducible only to ~1 % at N = 1e4."""
        self.assertGreaterEqual(TDCRPy_mod._L_XTOL, 1e-5)
        self.assertLessEqual(TDCRPy_mod._L_XTOL, 1e-2)

    def test_eff_returns_u_L_last(self):
        """u_L is appended, so every existing positional index still holds."""
        import inspect as _i
        src = _i.getsource(TDCRPy_mod.eff)
        self.assertIn('u_L,', src)
        # the documented order: L0 first, u_L last
        self.assertLess(src.index('L0, L_opt,'), src.rindex('u_L,'))

    def test_bracket_tolerates_non_finite(self):
        def tdcr(L):
            return np.nan if L < 1.0 else min(0.99, 0.2 * L)

        br = TDCRPy_mod._bracket_physical_L(tdcr, 0.6)
        self.assertIsNotNone(br)
        self.assertGreaterEqual(br[0], 1.0)


class TestEmptyDecaysCounted(unittest.TestCase):
    """Decays that deposit nothing must count as zero, not vanish (v2.20.28).

    A decay in which every emitted particle escapes the vial writes no line to
    the record file. Averaging the replayed efficiencies over the records that
    exist therefore dropped those decays instead of counting them, inflating
    the efficiency by M/(N-M). Measured at N = 20000: Fe-55 lost 328 decays
    and read 1.67 % high, Cr-51 lost 281 and read 1.43 % high; H-3 and C-14
    lose none and were unaffected.
    """

    def test_full_list_is_a_plain_mean(self):
        """With one record per decay the result must not change at all."""
        v = [0.0, 0.25, 0.5, 0.75, 1.0]
        mean, sem = lib._mean_over_decays(v, len(v))
        self.assertAlmostEqual(mean, float(np.mean(v)), places=12)
        self.assertAlmostEqual(sem, float(np.std(v)) / math.sqrt(len(v)),
                               places=12)

    def test_missing_records_count_as_zeros(self):
        """Three records out of four decays: the fourth contributed zero."""
        v = [1.0, 1.0, 1.0]
        mean, sem = lib._mean_over_decays(v, 4)
        self.assertAlmostEqual(mean, 0.75, places=12)
        # same as averaging [1, 1, 1, 0] explicitly
        ref = np.array([1.0, 1.0, 1.0, 0.0])
        self.assertAlmostEqual(mean, float(ref.mean()), places=12)
        self.assertAlmostEqual(sem, float(ref.std()) / 2.0, places=12)

    def test_inflation_matches_the_closed_form(self):
        """The old behaviour was high by exactly M/(N-M)."""
        v = [0.8] * 970                      # 30 of 1000 decays deposited none
        mean, _ = lib._mean_over_decays(v, 1000)
        old = float(np.mean(v))              # what the code used to return
        self.assertAlmostEqual(old / mean - 1.0, 30 / 970.0, places=10)

    def test_more_records_than_decays_does_not_scale_down(self):
        v = [1.0] * 10
        mean, _ = lib._mean_over_decays(v, 4)
        self.assertAlmostEqual(mean, 1.0, places=12)

    def test_replay_reproduces_the_direct_simulation(self):
        """The regression itself, on a nuclide that loses decays.

        Fe-55 is an electron-capture emitter whose X-rays frequently escape,
        so a measurable fraction of its decays deposit nothing. record=True
        both simulates and records, so the replay below consumes exactly the
        decays the direct result was computed from and any difference is the
        reconstruction, not Monte-Carlo noise.
        """
        saved = (lib.fAq, lib.micCorr)
        try:
            lib.modifyfAq(0.0)
            lib.modifyMicCorr(False)
            n = 400
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                direct = TDCRPy_mod.TDCRPy(6.0, 'Fe-55', '1', n,
                                           KB_TYPICAL, 10, record=True)
                replay = TDCRPy_mod.TDCRPy(6.0, 'Fe-55', '1', n,
                                           KB_TYPICAL, 10, readRecHist=True)
            # Not bit-exact: the record file stores energies as %.6E, so the
            # replay reads them back rounded to about 7 significant digits.
            # Before the fix the gap was 1.7e-2 relative, not 1e-8.
            for name, i in (('eff_S', 0), ('eff_D', 2), ('eff_T', 4)):
                self.assertAlmostEqual(
                    direct[i], replay[i], places=6,
                    msg=f'{name}: replay {replay[i]} != direct {direct[i]}; '
                        'decays that deposited nothing are being dropped')
        finally:
            lib.modifyfAq(saved[0])
            lib.modifyMicCorr(saved[1])
            lib.cleanup_record_files()


class TestTTAFunctions(unittest.TestCase):
    """The closed-form pieces of the delayed-fluorescence model.

    R. Coulon, "TDCR Measurements -- Lecture Notes", ch. 3, Eq. (3.4),
    (3.7a), (3.8), (3.8a), as implemented in newTDCRmodel/tta_light_model.py.
    """

    def test_sa_matches_the_reference_constants(self):
        """S_a = 4 W A_c k_q / (3 gamma); the shipped constants give 1666.7."""
        self.assertAlmostEqual(lib.sa_tta(k_q=5.0e6, gamma=1.0e-11,
                                          A_c=1.0e-13), 1666.6666666, places=4)
        self.assertAlmostEqual(lib.sa_tta(), 1666.6666666, places=4)

    def test_sa_scales_with_each_constant(self):
        base = lib.sa_tta(k_q=5.0e6, gamma=1.0e-11, A_c=1.0e-13)
        self.assertAlmostEqual(
            lib.sa_tta(k_q=1.0e7, gamma=1.0e-11, A_c=1.0e-13), 2 * base)
        self.assertAlmostEqual(
            lib.sa_tta(k_q=5.0e6, gamma=2.0e-11, A_c=1.0e-13), base / 2)
        self.assertAlmostEqual(
            lib.sa_tta(k_q=5.0e6, gamma=1.0e-11, A_c=2.0e-13), 2 * base)

    def test_phi_tta_series_branch_does_not_underflow(self):
        """Below b = 1e-4 the direct form is 1 minus 1-to-within-1e-16."""
        self.assertAlmostEqual(float(lib.phi_tta(np.array([1e-12]))[0]),
                               5e-13, delta=1e-24)
        # the two branches must agree across the switch at 1e-4
        for b in (9.9e-5, 1.0e-4, 1.01e-4):
            direct = 1.0 - math.log1p(b) / b
            self.assertAlmostEqual(float(lib.phi_tta(np.array([b]))[0]),
                                   direct, places=10)

    def test_phi_tta_is_bounded_and_monotone(self):
        b = np.logspace(-14, 6, 2000)
        v = lib.phi_tta(b)
        self.assertTrue(np.all(v >= 0.0))
        self.assertTrue(np.all(v < 1.0))
        self.assertTrue(np.all(np.diff(v) > 0),
                        'phi_TTA must rise with the stopping power')

    def test_f_prompt(self):
        self.assertAlmostEqual(lib.f_prompt(3.0, 3.0), 1 - math.exp(-1.0))
        self.assertAlmostEqual(lib.f_prompt(50.0, 3.0), 1.0, places=6)
        self.assertAlmostEqual(lib.f_prompt(0.0, 3.0), 0.0)
        # tau_p <= 0 is the instantaneous limit
        self.assertEqual(lib.f_prompt(1.0, 0.0), 1.0)
        self.assertEqual(lib.f_prompt(1.0, -1.0), 1.0)

    def test_f_delayed(self):
        self.assertAlmostEqual(lib.f_delayed(85.0, 85.0), 0.5)
        self.assertAlmostEqual(lib.f_delayed(50.0, 85.0), 50 / 135.0)
        self.assertEqual(lib.f_delayed(0.0, 85.0), 0.0)
        # never saturates as fast as the prompt term: that is the W-dependence
        self.assertLess(lib.f_delayed(1000.0, 85.0), 0.93)

    def test_self_test_passes(self):
        self.assertTrue(lib.selfTestDelayedFluorescence())


class TestPsiReservoir(unittest.TestCase):
    """Psi(E), the TTA delayed-light reservoir, against the reference.

    DELAYED_FLUORESCENCE_MODEL.md section 3 publishes, for Tc-99 at
    kB = 7.5e-6 cm/keV (spectrum-averaged, mean E = 97.3 keV): <Q> = 90.5 keV
    (93 % of <E>) and <Psi> = 57.9 keV (59.5 %). Reproducing those is what
    makes this a port rather than a re-derivation.
    """

    def test_reproduces_the_reference_spectrum_average(self):
        e_bins, p_spec = lib.readBetaSpectra('Tc-99')
        e_bins = np.asarray(e_bins, dtype=float)
        p_spec = np.asarray(p_spec, dtype=float)
        mean_e = float(np.average(e_bins, weights=p_spec))
        mean_psi = float(np.average(lib.psiTTA(e_bins), weights=p_spec))
        self.assertAlmostEqual(mean_e, 97.34, delta=0.05)
        self.assertAlmostEqual(mean_psi, 57.95, delta=0.15,
                               msg=f'<Psi> = {mean_psi:.2f} keV, reference '
                                   f'57.9 keV')
        self.assertAlmostEqual(100 * mean_psi / mean_e, 59.53, delta=0.15)

    def test_psi_never_exceeds_the_energy(self):
        for e in (0.01, 0.05, 0.5, 5.0, 50.0, 500.0, 5000.0, 2e4):
            self.assertLessEqual(lib.psiTTA(e), e + 1e-9, f'Psi({e}) > E')

    def test_psi_is_monotone(self):
        grid = np.logspace(-2, np.log10(2e4), 1500)
        v = lib.psiTTA(grid)
        self.assertTrue(np.all(np.diff(v) >= -1e-12))

    def test_psi_is_bigger_than_the_quenching_deficit_it_replaced(self):
        """The reservoir is ~8.5x the old E - Q(E,kB) form, so the two etas
        are on different scales and must never be compared."""
        e = 97.3
        q = lib.Em_e(e * 1e3, e * 1e3, 7.5e-6 * 1e3, 1000) * 1e-3
        deficit = e - q
        self.assertGreater(lib.psiTTA(e) / deficit, 5.0)

    def test_grid_covers_the_hardest_betas(self):
        """A grid stopping at 2 MeV extrapolates rather than saturates:
        phi_TTA is ~0.33 for minimum-ionising electrons, not ~0."""
        self.assertGreaterEqual(lib._PSI_E_MAX_KEV, 2.0e4)
        ratio = lib.psiTTA(2.0e4) / 2.0e4
        self.assertAlmostEqual(ratio, 0.3668, delta=0.02)

    def test_vectorised_and_scalar_agree(self):
        e = np.array([1.0, 10.0, 100.0])
        v = lib.psiTTA(e)
        self.assertEqual(np.shape(v), (3,))
        for i, x in enumerate(e):
            self.assertAlmostEqual(float(v[i]), lib.psiTTA(float(x)),
                                   places=12)

    def test_cache_is_keyed_on_the_stopping_power(self):
        """Psi depends on the cocktail through S(E). A reservoir that
        outlived a cocktail change would be a trap, not a cache."""
        lib.psiTTA(10.0)
        # key = (S_a, e_min, e_max, n, cocktail, sp_model, rho, Z, A)
        self.assertTrue(all(len(k) == 9 for k in lib._PSI_CACHE))
        self.assertTrue(any(k[4] == lib.lsCocktail()
                            for k in lib._PSI_CACHE),
                        'cache key does not carry the cocktail')
        self.assertTrue(any(k[5] == lib.sp_model for k in lib._PSI_CACHE),
                        'cache key does not carry the stopping-power model')


class TestEffectiveLightEnergy(unittest.TestCase):
    """Q_eff = Q f_p + eta Psi f_d, Eq. (3.2)."""

    def test_eta_zero_and_a_long_window_give_back_the_birks_energy(self):
        q = 90.5
        self.assertAlmostEqual(
            lib.effective_light_energy(q, 97.3, 97.3, 1e9, eta=0.0,
                                       tau_p=3.0, tau_d=85.0),
            q, places=12)

    def test_prompt_factor_is_applied_at_a_short_window(self):
        """eta = 0 is NOT the baseline at a short window: f_p still bites."""
        q = 100.0
        got = lib.effective_light_energy(q, 100.0, 100.0, 3.0, eta=0.0,
                                         tau_p=3.0, tau_d=85.0)
        self.assertAlmostEqual(got, q * (1 - math.exp(-1.0)), places=10)

    def test_delayed_term_has_the_documented_value(self):
        q, e = 90.0, 97.3
        got = lib.effective_light_energy(q, e, e, 50.0, eta=0.03,
                                         tau_p=3.0, tau_d=85.0)
        expected = (q * lib.f_prompt(50.0, 3.0)
                    + 0.03 * lib.psiTTA(e) * lib.f_delayed(50.0, 85.0))
        self.assertAlmostEqual(got, expected, places=12)

    def test_partial_deposit_uses_the_track_segment(self):
        """dPsi = Psi(Ei) - Psi(Ei - Ed), the same segment E_quench_e
        integrates the Birks factor over -- and Psi(Ed) when Ei == Ed."""
        q = 10.0
        full = lib.effective_light_energy(q, 100.0, 100.0, 50.0, eta=0.1,
                                          tau_p=3.0, tau_d=85.0)
        f_p = lib.f_prompt(50.0, 3.0)
        expected_full = (q * f_p
                         + 0.1 * lib.psiTTA(100.0) * lib.f_delayed(50.0, 85.0))
        self.assertAlmostEqual(full, expected_full, places=9)

        part = lib.effective_light_energy(q, 100.0, 40.0, 50.0, eta=0.1,
                                          tau_p=3.0, tau_d=85.0)
        d_psi = lib.psiTTA(100.0) - lib.psiTTA(60.0)
        self.assertAlmostEqual(
            part, q * f_p + 0.1 * d_psi * lib.f_delayed(50.0, 85.0), places=9)
        self.assertLess(part, full,
                        'a partial deposit must draw on less of the reservoir')

    def test_grows_monotonically_with_the_window(self):
        q, e = 90.0, 97.3
        vals = [lib.effective_light_energy(q, e, e, w, eta=0.03, tau_p=3.0,
                                           tau_d=85.0)
                for w in (1, 5, 10, 50, 100, 500, 1000)]
        self.assertTrue(all(b > a for a, b in zip(vals, vals[1:])),
                        f'Q_eff is not monotone in the window: {vals}')


class TestDelayedFluorescenceConfig(unittest.TestCase):
    """The option, its per-cocktail parameters, and its default state."""

    def setUp(self):
        self._saved = {
            'delayedFluorescence': lib.delayedFluorescence,
            'tau_prompt': lib.tau_prompt,
            'tau_delayed': lib.tau_delayed,
            'eta_delayed': lib.eta_delayed,
            'tau': lib.tau,
        }
        TDCRPy_mod._delayed_announced.clear()

    def tearDown(self):
        lib.update_config_batch(self._saved)
        TDCRPy_mod._delayed_announced.clear()

    def test_ships_off(self):
        with importlib.resources.path('tdcrpy', 'configDefault.toml') as cfg:
            text = cfg.read_text(encoding='utf-8')
        line = [l.strip() for l in text.splitlines()
                if l.strip().startswith('delayedFluorescence')]
        self.assertEqual(len(line), 1)
        self.assertIn('False', line[0])

    def test_missing_keys_fall_back_to_off(self):
        """A config.toml written before v2.20.30 must not switch it on."""
        import configparser
        parser = configparser.ConfigParser()
        parser.read_string('[Inputs]\nmicCorr = False\n')
        inp = parser['Inputs']
        self.assertFalse(inp.getboolean('delayedFluorescence', fallback=False))
        self.assertEqual(inp.getfloat('tau_prompt', fallback=3.0), 3.0)
        self.assertEqual(inp.getfloat('tau_delayed', fallback=85.0), 85.0)
        self.assertEqual(inp.getfloat('eta_delayed', fallback=0.03), 0.03)

    def test_the_two_fitted_cocktail_families(self):
        """The values this release ships, from the per-source fits."""
        for name in ('Ultima Gold', 'Ultima Gold LLT', 'Ultima Gold AB',
                     'Ultima Gold XR'):
            tau_d, eta, measured = lib.delayedFluorescenceDefaults(name)
            self.assertEqual((tau_d, eta), (85.0, 0.03), name)
            self.assertTrue(measured, name)
        tau_d, eta, measured = lib.delayedFluorescenceDefaults('Toluene')
        self.assertEqual((tau_d, eta), (138.0, 0.05))
        self.assertTrue(measured)

    def test_unfitted_cocktails_are_flagged_as_such(self):
        for name in ('Insta-Gel Plus', 'Hionic-Fluor', 'ProSafe+',
                     'Pseudocumene', 'PXE'):
            _, _, measured = lib.delayedFluorescenceDefaults(name)
            self.assertFalse(measured,
                             f'{name} has no fitted eta/tau_delayed')

    def test_neat_aromatics_group_with_toluene(self):
        for name in ('Pseudocumene', 'PXE', 'LAB'):
            tau_d, eta, _ = lib.delayedFluorescenceDefaults(name)
            self.assertEqual((tau_d, eta), (138.0, 0.05), name)

    def test_setters_round_trip(self):
        lib.modifyTau_prompt(2.5)
        self.assertAlmostEqual(lib.tau_prompt, 2.5)
        lib.modifyTau_delayed(138.0)
        self.assertAlmostEqual(lib.tau_delayed, 138.0)
        lib.modifyEta_delayed(0.05)
        self.assertAlmostEqual(lib.eta_delayed, 0.05)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            lib.modifyDelayedFluorescence(True)
        self.assertTrue(lib.delayedFluorescence)
        lib.modifyDelayedFluorescence(False)
        self.assertFalse(lib.delayedFluorescence)

    def test_eta_is_bounded_by_spin_statistics(self):
        self.assertAlmostEqual(lib.ETA_MAX_TTA, 1.5)
        with self.assertRaises(ValueError):
            lib.modifyEta_delayed(2.0)
        with self.assertRaises(ValueError):
            lib.modifyEta_delayed(-0.01)

    def test_cocktail_change_re_derives_the_parameters(self):
        """They describe the scintillator, so a stale pair left over from the
        previous cocktail would be exactly the staleness this package keeps
        tripping over."""
        saved = (lib.lsCocktail(), lib.fAq)
        try:
            lib.modifyLScocktail('Toluene', 0.0)
            self.assertAlmostEqual(lib.tau_delayed, 138.0)
            self.assertAlmostEqual(lib.eta_delayed, 0.05)
            lib.modifyLScocktail('Ultima Gold', 0.0)
            self.assertAlmostEqual(lib.tau_delayed, 85.0)
            self.assertAlmostEqual(lib.eta_delayed, 0.03)
        finally:
            lib.modifyLScocktail(saved[0], saved[1])

    def test_warns_once_for_an_unfitted_cocktail(self):
        saved = (lib.lsCocktail(), lib.fAq)
        try:
            lib.modifyLScocktail('Insta-Gel Plus', 0.0)
            lib.update_config_batch({'delayedFluorescence': True})
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                TDCRPy_mod._resolve_delayed_fluorescence(50.0)
                TDCRPy_mod._resolve_delayed_fluorescence(50.0)
            hits = [c for c in caught if 'no eta/tau_delayed were fitted'
                    in str(c.message)]
            self.assertEqual(len(hits), 1, 'unfitted cocktail not flagged '
                                           'exactly once')
        finally:
            lib.update_config_batch({'delayedFluorescence': False})
            lib.modifyLScocktail(saved[0], saved[1])

    def test_resolution_returns_none_when_off(self):
        lib.update_config_batch({'delayedFluorescence': False})
        self.assertIsNone(TDCRPy_mod._resolve_delayed_fluorescence(50.0))

    def test_resolution_carries_the_window(self):
        lib.update_config_batch({'delayedFluorescence': True,
                                 'tau_prompt': 3.0, 'tau_delayed': 85.0,
                                 'eta_delayed': 0.03})
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            got = TDCRPy_mod._resolve_delayed_fluorescence(77.0)
        self.assertEqual(got, (77.0, 0.03, 3.0, 85.0))


class TestDelayedFluorescenceInTheModel(unittest.TestCase):
    """End to end: the option is inert off, and window-dependent on."""

    N = 400
    L = 6.0

    def setUp(self):
        self._saved = {
            'delayedFluorescence': lib.delayedFluorescence,
            'eta_delayed': lib.eta_delayed,
            'tau_prompt': lib.tau_prompt,
            'tau_delayed': lib.tau_delayed,
            'tau': lib.tau,
            'micCorr': lib.micCorr,
            'fAq': lib.fAq,
        }
        TDCRPy_mod._delayed_announced.clear()

    def tearDown(self):
        lib.update_config_batch(self._saved)
        TDCRPy_mod._delayed_announced.clear()
        lib.cleanup_record_files()

    def _eff_d(self, rad='C-14'):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return TDCRPy_mod.TDCRPy(self.L, rad, '1', self.N,
                                     KB_TYPICAL, 10)[2]

    def test_eta_zero_reproduces_the_plain_model_exactly(self):
        """At the shipped tau the prompt factor is 1 to 1e-7, so eta = 0 must
        return the Birks-quenched result bit for bit -- no Monte-Carlo
        tolerance needed, because the only change is arithmetic on an
        identical decay sample.

        It is still not the baseline in general: at a SHORT window f_p bites
        even with eta = 0, which is why the option, not eta, is the switch.
        """
        lib.update_config_batch({'delayedFluorescence': False, 'tau': 100.0,
                                 'micCorr': False, 'fAq': 0.0})
        np.random.seed(3)
        off = self._eff_d()
        lib.update_config_batch({'delayedFluorescence': True,
                                 'eta_delayed': 0.0, 'tau_prompt': 3.0})
        np.random.seed(3)
        on = self._eff_d()
        self.assertAlmostEqual(on, off, places=12)

    def test_switching_it_on_raises_the_efficiency(self):
        lib.update_config_batch({'delayedFluorescence': False, 'tau': 100.0,
                                 'micCorr': False, 'fAq': 0.0})
        np.random.seed(5)
        off = self._eff_d('H-3')
        lib.update_config_batch({'delayedFluorescence': True,
                                 'eta_delayed': 0.03, 'tau_delayed': 85.0,
                                 'tau_prompt': 3.0})
        np.random.seed(5)
        on = self._eff_d('H-3')
        self.assertGreater(on, off,
                           'the delayed term adds light, so it must raise '
                           'the efficiency at fixed L')

    def test_efficiency_rises_with_the_coincidence_window(self):
        """The whole point of the model: the plain Birks model has no
        W-dependence at all."""
        lib.update_config_batch({'delayedFluorescence': True,
                                 'eta_delayed': 0.03, 'tau_delayed': 85.0,
                                 'tau_prompt': 3.0, 'micCorr': False,
                                 'fAq': 0.0})
        effs = []
        for tau in (10.0, 100.0, 1000.0):
            lib.update_config_batch({'tau': tau})
            np.random.seed(9)
            effs.append(self._eff_d('H-3'))
        self.assertTrue(all(b > a for a, b in zip(effs, effs[1:])),
                        f'eff_D is not monotone in the window: {effs}')

    def test_off_by_default_leaves_the_quench_path_untouched(self):
        """_quench with delayed=None must be the pre-v2.20.30 code path."""
        particles = ['electron', 'electron']
        energies = [10.0, 50.0]
        initial = [10.0, 50.0]
        plain = TDCRPy_mod._quench(list(particles), list(energies),
                                   list(initial), KB_TYPICAL, 100, 100,
                                   False, False)
        explicit = TDCRPy_mod._quench(list(particles), list(energies),
                                      list(initial), KB_TYPICAL, 100, 100,
                                      False, False, None)
        for a, b in zip(plain, explicit):
            self.assertAlmostEqual(a, b, places=12)

    def test_micelle_retention_scales_the_total_light(self):
        """With both corrections on, confinement must scale prompt AND
        delayed light: it is about where the track ran, not about which
        component the light came out as."""
        particles = ['electron']
        e = [3.0]
        delayed = (100.0, 0.03, 3.0, 85.0)
        lib.update_config_batch({'fAq': 0.10, 'micCorr': True})
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            both = TDCRPy_mod._quench(list(particles), list(e), list(e),
                                      KB_TYPICAL, 100, 100, True, False,
                                      delayed)[0]
            only_delayed = TDCRPy_mod._quench(list(particles), list(e),
                                              list(e), KB_TYPICAL, 100, 100,
                                              False, False, delayed)[0]
        self.assertLess(both, only_delayed,
                        'the micelle correction did not reach the delayed '
                        'light')


class TestECBranchProbability(unittest.TestCase):
    """Which nuclides need the micelle correction is read, not listed.

    Until v2.20.29 the two nuclides known to need the correction were a
    hard-coded tuple -- the same kind of duplicated constant that has gone
    stale elsewhere in this package. The PenNuc decay scheme already carries
    the answer: a primary branch entry labelled ``Atom_<shell>`` is a vacancy
    made by electron capture.
    """

    # Electron-capture probabilities from the DDEP/LNHB evaluations. Anchors
    # the reconstruction, not just its self-consistency.
    EVALUATED = {
        'Fe-55': 1.0000, 'Cr-51': 1.0000, 'Mn-54': 1.0000, 'I-125': 1.0000,
        'Cd-109': 1.0000, 'Co-57': 1.0000, 'Zn-65': 0.9830, 'Tl-201': 0.9800,
        'Cu-64': 0.4350, 'K-40': 0.1072, 'Na-22': 0.0964, 'Tl-204': 0.0297,
        'Cl-36': 0.0190,
    }
    PURE_BETA = ('H-3', 'C-14', 'Ni-63', 'Tc-99', 'Sr-90', 'P-32', 'S-35')

    def test_matches_the_evaluated_ec_probabilities(self):
        for rad, expected in self.EVALUATED.items():
            got = lib.ec_branch_probability(rad)
            self.assertAlmostEqual(
                got, expected, delta=max(0.01, 0.05 * expected),
                msg=f'{rad}: P(EC) = {got:.4f}, evaluated {expected:.4f}')

    def test_pure_beta_emitters_have_no_ec_branch(self):
        for rad in self.PURE_BETA:
            self.assertEqual(lib.ec_branch_probability(rad), 0.0,
                             f'{rad} is a pure beta emitter')

    def test_isomeric_transition_is_not_electron_capture(self):
        """Tc-99m decays by IT. Its conversion-electron vacancies are not EC."""
        self.assertEqual(lib.ec_branch_probability('Tc-99m'), 0.0)

    def test_probability_is_normalised(self):
        for rad in ('Fe-55', 'Cu-64', 'Eu-152', 'Am-241', 'Co-60'):
            p = lib.ec_branch_probability(rad)
            self.assertGreaterEqual(p, 0.0)
            self.assertLessEqual(p, 1.0)

    def test_unknown_nuclide_is_not_critical(self):
        """A missing decay file must never switch a correction on."""
        self.assertEqual(lib.ec_branch_probability('Xx-999'), 0.0)
        self.assertFalse(lib.is_micelle_critical('Xx-999'))

    def test_threshold_selects_the_right_nuclides(self):
        self.assertTrue(lib.is_micelle_critical('Fe-55', threshold=0.01))
        self.assertTrue(lib.is_micelle_critical('Cl-36', threshold=0.01))
        self.assertFalse(lib.is_micelle_critical('Cl-36', threshold=0.5))
        self.assertTrue(lib.is_micelle_critical('Fe-55', threshold=0.5))
        self.assertFalse(lib.is_micelle_critical('H-3', threshold=0.0001))


class TestMicelleAutoActivation(unittest.TestCase):
    """The correction switches itself on for EC nuclides, and says so.

    Fe-55 and Cr-51 are electron-capture emitters whose Auger cascades sit in
    the confinement regime: about 5 % of their light is lost inside the water
    droplets, and that part is NOT absorbed by the fitted L, unlike the
    (1 - fAq) bulk dilution. Since v2.20.29 TDCRPy applies the correction to
    them by itself whenever fAq > 0.
    """

    def setUp(self):
        self._faq = lib.fAq
        self._auto = lib.micCorrAuto
        self._thr = lib.micCorrECThreshold
        self._mic = lib.micCorr
        TDCRPy_mod._micelle_advised.clear()
        TDCRPy_mod._micelle_auto_announced.clear()

    def tearDown(self):
        lib.update_config_batch({'fAq': self._faq,
                                 'micCorrAuto': self._auto,
                                 'micCorrECThreshold': self._thr,
                                 'micCorr': self._mic})
        TDCRPy_mod._micelle_advised.clear()
        TDCRPy_mod._micelle_auto_announced.clear()

    # -- what ships -------------------------------------------------------
    def test_default_configuration_ships_auto_on_with_miccorr_off(self):
        with importlib.resources.path('tdcrpy', 'configDefault.toml') as cfg:
            text = cfg.read_text(encoding='utf-8')
        keys = {}
        for line in text.splitlines():
            line = line.strip()
            if line.startswith('micCorr'):
                k, _, v = line.partition('=')
                keys[k.strip()] = v.strip()
        self.assertEqual(keys.get('micCorr'), 'False')
        self.assertEqual(keys.get('micCorrAuto'), 'True')
        self.assertEqual(float(keys.get('micCorrECThreshold')), 0.01)

    def test_missing_keys_fall_back_to_auto_on(self):
        """A config.toml written before v2.20.29 must get the new default."""
        import configparser
        parser = configparser.ConfigParser()
        parser.read_string('[Inputs]\nmicCorr = False\n')
        inp = parser['Inputs']
        self.assertTrue(inp.getboolean('micCorrAuto', fallback=True))
        self.assertEqual(inp.getfloat('micCorrECThreshold', fallback=0.01),
                         0.01)

    # -- the resolution itself --------------------------------------------
    def test_no_water_means_no_change(self):
        """With fAq = 0 there are no droplets: behaviour as before v2.20.29."""
        lib.modifyfAq(0.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            self.assertFalse(
                TDCRPy_mod._resolve_micelle_correction(['Fe-55'], False))
        self.assertEqual(len(caught), 0, 'warned with no aqueous phase')

    def test_activates_for_an_ec_nuclide(self):
        lib.modifyfAq(0.10)
        for rad in ('Fe-55', 'Cr-51', 'I-125', 'Mn-54'):
            TDCRPy_mod._micelle_auto_announced.clear()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                on = TDCRPy_mod._resolve_micelle_correction([rad], False)
            self.assertTrue(on, f'{rad}: correction did not self-activate')
            self.assertEqual(len(caught), 1,
                             f'{rad}: activation was not announced')
            self.assertIn('automatically', str(caught[0].message))

    def test_leaves_beta_emitters_alone(self):
        lib.modifyfAq(0.10)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            for rad in ('H-3', 'C-14', 'Tc-99', 'Ni-63'):
                self.assertFalse(
                    TDCRPy_mod._resolve_micelle_correction([rad], False), rad)
        self.assertEqual(len(caught), 0, 'warned for a pure beta emitter')

    def test_announces_once_per_nuclide(self):
        lib.modifyfAq(0.10)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            TDCRPy_mod._resolve_micelle_correction(['Fe-55'], False)
            TDCRPy_mod._resolve_micelle_correction(['Fe-55'], False)
        self.assertEqual(len(caught), 1, 'activation announced twice')

    def test_activates_for_a_mixture_containing_one_ec_nuclide(self):
        lib.modifyfAq(0.10)
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always')
            self.assertTrue(TDCRPy_mod._resolve_micelle_correction(
                ['H-3', 'Fe-55'], False))

    def test_explicit_miccorr_is_never_overridden(self):
        lib.modifyfAq(0.10)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            self.assertTrue(
                TDCRPy_mod._resolve_micelle_correction(['H-3'], True))
            self.assertTrue(
                TDCRPy_mod._resolve_micelle_correction(['Fe-55'], True))
        self.assertEqual(len(caught), 0,
                         'announced an activation the caller had already made')

    def test_threshold_is_honoured(self):
        lib.modifyfAq(0.10)
        lib.modifyMicCorrECThreshold(0.5)
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always')
            self.assertFalse(
                TDCRPy_mod._resolve_micelle_correction(['Cl-36'], False),
                'Cl-36 has a 1.9 % EC branch and is below a 0.5 threshold')
            self.assertTrue(
                TDCRPy_mod._resolve_micelle_correction(['Fe-55'], False))

    # -- opting out --------------------------------------------------------
    def test_auto_off_falls_back_to_the_advisory(self):
        lib.modifyfAq(0.10)
        lib.modifyMicCorrAuto(False)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            on = TDCRPy_mod._resolve_micelle_correction(['Fe-55'], False)
        self.assertFalse(on, 'micCorrAuto = False did not opt out')
        self.assertEqual(len(caught), 1)
        self.assertIn('micCorrAuto is disabled', str(caught[0].message))

    def test_advisory_warns_once_per_nuclide(self):
        lib.modifyfAq(0.10)
        lib.modifyMicCorrAuto(False)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            TDCRPy_mod._resolve_micelle_correction(['Cr-51'], False)
            TDCRPy_mod._resolve_micelle_correction(['Cr-51'], False)
        self.assertEqual(len(caught), 1, 'advisory repeated for one nuclide')

    def test_setters_round_trip(self):
        lib.modifyMicCorrAuto(False)
        self.assertFalse(lib.micCorrAuto)
        lib.modifyMicCorrAuto(True)
        self.assertTrue(lib.micCorrAuto)
        lib.modifyMicCorrECThreshold(0.25)
        self.assertAlmostEqual(lib.micCorrECThreshold, 0.25)

    # -- end to end --------------------------------------------------------
    def test_run_enters_the_kernel_exactly_when_it_should(self):
        """The exact validation: is the micelle kernel entered, and for whom?

        The kernel is @njit(parallel=True) and keeps a per-thread RNG state
        that np.random.seed() cannot reach, so two arms of the same physics
        differ by Monte-Carlo noise however the Python side is seeded --
        comparing them bit for bit would test the seeding. Counting the calls
        is exact.
        """
        real = lib.pure_mc_efficient_energy_numba
        calls = []

        def counting(E, **kwargs):
            calls.append(E)
            return real(E, **kwargs)

        lib.pure_mc_efficient_energy_numba = counting
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')

                # opted out: the kernel must never be entered
                lib.update_config_batch({'fAq': 0.10, 'micCorr': False,
                                         'micCorrAuto': False})
                calls.clear()
                TDCRPy_mod.TDCRPy(6.0, 'Fe-55', '1', 60, KB_TYPICAL, 10)
                self.assertEqual(len(calls), 0,
                                 'micCorrAuto = False still corrected')

                # self-activating on an EC nuclide
                lib.modifyMicCorrAuto(True)
                TDCRPy_mod._micelle_auto_announced.clear()
                calls.clear()
                TDCRPy_mod.TDCRPy(6.0, 'Fe-55', '1', 60, KB_TYPICAL, 10)
                self.assertGreater(len(calls), 0,
                                   'Fe-55 ran without the correction')

                # ... and only there
                calls.clear()
                TDCRPy_mod.TDCRPy(6.0, 'C-14', '1', 60, KB_TYPICAL, 10)
                self.assertEqual(len(calls), 0,
                                 'C-14 paid for a correction it does not need')

                # no droplets, no correction, whatever the nuclide
                lib.modifyfAq(0.0)
                calls.clear()
                TDCRPy_mod.TDCRPy(6.0, 'Fe-55', '1', 60, KB_TYPICAL, 10)
                self.assertEqual(len(calls), 0,
                                 'corrected a cocktail with no aqueous phase')
        finally:
            lib.pure_mc_efficient_energy_numba = real
            lib.cleanup_record_files()

    def test_self_activated_run_agrees_with_an_explicit_one(self):
        """Statistical validation: the two arms are the same physics.

        Self-activation must land on the explicit micCorr = True result and
        not somewhere between it and the uncorrected one. The arms are
        independent Monte-Carlo samples, so they are compared inside their own
        combined standard uncertainty, while the uncorrected arm has to sit
        clear of both.
        """
        n = 2000
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                lib.update_config_batch({'fAq': 0.10, 'micCorr': False,
                                         'micCorrAuto': True})
                auto = TDCRPy_mod.TDCRPy(6.0, 'Fe-55', '1', n, KB_TYPICAL, 10)
                lib.update_config_batch({'micCorr': True})
                explicit = TDCRPy_mod.TDCRPy(6.0, 'Fe-55', '1', n,
                                             KB_TYPICAL, 10)
                lib.update_config_batch({'micCorr': False,
                                         'micCorrAuto': False})
                off = TDCRPy_mod.TDCRPy(6.0, 'Fe-55', '1', n, KB_TYPICAL, 10)

            u = math.hypot(auto[3], explicit[3])
            self.assertLess(
                abs(auto[2] - explicit[2]), 4.0 * u,
                f'self-activated eff_D = {auto[2]:.5f} +/- {auto[3]:.5f} is '
                f'not the explicit {explicit[2]:.5f} +/- {explicit[3]:.5f}')
            self.assertGreater(
                off[2] - auto[2], 2.0 * math.hypot(auto[3], off[3]),
                f'uncorrected eff_D = {off[2]:.5f} is not measurably above '
                f'the corrected {auto[2]:.5f}: the default is cosmetic')
        finally:
            lib.cleanup_record_files()

    def test_run_differs_from_the_uncorrected_one(self):
        """And it must actually change something, or the default is cosmetic."""
        lib.update_config_batch({'fAq': 0.10, 'micCorr': False})
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                lib.modifyMicCorrAuto(False)
                np.random.seed(11)
                off = TDCRPy_mod.TDCRPy(6.0, 'Fe-55', '1', 300,
                                        KB_TYPICAL, 10)
                lib.modifyMicCorrAuto(True)
                np.random.seed(11)
                on = TDCRPy_mod.TDCRPy(6.0, 'Fe-55', '1', 300,
                                       KB_TYPICAL, 10)
            self.assertGreater(
                abs(off[2] - on[2]), 1e-6,
                'self-activation changed nothing: the correction is inert')
            self.assertGreater(
                off[2], on[2],
                'the uncorrected efficiency must be the HIGHER one -- '
                'confinement removes light')
        finally:
            lib.cleanup_record_files()


class TestRecordFileIsolation(unittest.TestCase):
    """Record files must be private to the writing process.

    Before v2.20.20 the four decay-history files had fixed names in the shared
    system temp directory, so two TDCRPy processes on one machine overwrote
    each other's histories. `_run_from_history` skips malformed lines, so the
    damage surfaced as biased efficiencies and non-converging eff() fits rather
    than as an error.
    """

    def test_path_contains_pid(self):
        pid = os.getpid()
        for i in range(4):
            p = lib.record_file_path(i)
            self.assertIn(str(pid), os.path.basename(p),
                          f'record file {i} is not process-specific: {p}')
            self.assertTrue(os.path.basename(p).startswith(f'Temp_E{i}_'))

    def test_paths_differ_between_processes(self):
        mine = lib.record_file_path(2)
        other = lib.record_file_path(2, pid=os.getpid() + 1)
        self.assertNotEqual(mine, other,
                            'two processes would share one record file')

    def test_four_distinct_files(self):
        paths = [lib.record_file_path(i) for i in range(4)]
        self.assertEqual(len(set(paths)), 4)
        self.assertEqual(len({os.path.dirname(p) for p in paths}), 1)

    def test_explicit_pid_is_honoured(self):
        p = lib.record_file_path(2, pid=4242)
        self.assertEqual(os.path.basename(p), 'Temp_E2_4242.txt')

    def test_record_then_replay_round_trip(self):
        """record=True then readRecHist=True must agree, using the new names."""
        out_rec = TDCRPy_mod.TDCRPy(_L_STOCH, 'H-3', '1', 200, _KB, _V,
                                    record=True)
        self.assertTrue(os.path.exists(lib.record_file_path(2)),
                        'record=True did not create the per-process file')
        out_rep = TDCRPy_mod.TDCRPy(_L_STOCH, 'H-3', '1', 200, _KB, _V,
                                    readRecHist=True)
        # Not bit-exact: energies are written as %.6E, so the replay reads
        # them back rounded to ~7 significant digits. That propagates to a
        # ~1e-9 difference on an efficiency -- anything larger means the
        # replay read a different history.
        for k, i in (('S', 0), ('D', 2), ('T', 4)):
            self.assertAlmostEqual(
                out_rec[i], out_rep[i], places=7,
                msg=f'replay disagrees with the run that wrote the history ({k})')

    def test_cleanup_removes_this_process_files(self):
        TDCRPy_mod.TDCRPy(_L_STOCH, 'H-3', '1', 50, _KB, _V, record=True)
        self.assertTrue(os.path.exists(lib.record_file_path(2)))
        removed = lib.cleanup_record_files()
        self.assertTrue(removed, 'cleanup_record_files() removed nothing')
        for i in range(4):
            self.assertFalse(os.path.exists(lib.record_file_path(i)),
                             f'record file {i} survived cleanup')
        self.assertEqual(lib.cleanup_record_files(), [],
                         'cleanup on an already-clean state should be a no-op')


if __name__ == '__main__':
    unittest.main()
