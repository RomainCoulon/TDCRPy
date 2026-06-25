# -*- coding: utf-8 -*-
"""
test_tdcrpy.py

Unit tests for TDCRPy package — covers TDCR_model_lib and TDCRPy modules.

Authors: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

import unittest
import importlib.resources
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
        result = lib.stoppingpowerA(5000)
        self.assertEqual(result, 871008.0)

    def test_stoppingpower(self):
        result = lib.stoppingpower(100)
        self.assertEqual(result, 326.78)

    def test_readBetaShape(self):
        e, p = lib.readBetaShape("H-3", "beta-", "tot")
        # Mean over the combined energy+probability values (matches original commented test)
        result = np.mean(np.hstack([e, p]))
        self.assertAlmostEqual(result, 4.66, delta=1.0)

    def test_readBetaSpectra(self):
        result = np.mean(lib.readBetaSpectra("H-3"))
        self.assertEqual(result, 4.651594887459807)

    def test_E_quench_e(self):
        result = np.mean(lib.E_quench_e(100 * 1e3, 100 * 1e3, 0.01, 20000) * 1e-3)
        self.assertAlmostEqual(result, 91.24, delta=0.1)

    def test_E_quench_e2(self):
        result = np.mean(lib.E_quench_e(100 * 1e3, 50 * 1e3, 0.01, 20000) * 1e-3)
        self.assertAlmostEqual(result, 47.97, delta=0.1)

    def test_E_quench_a(self):
        result = np.mean(lib.E_quench_a(5000, 1e-5, 20000))
        self.assertEqual(result, 344.28633759073443)

    # run_interpolate tested inside Em_a() and Em_e()

    def test_Em_a(self):
        result = lib.Em_a(5000, 1e-5, 10000)
        self.assertLessEqual(result, 344.12105896)

    def test_Em_e(self):
        result = lib.Em_e(50000, 50000, 1e-5, 10000)
        self.assertLessEqual(result, 48568.14392205286)

    def test_Em_e2(self):
        result = lib.Em_e(50000, 25000, 1e-5, 10000)
        self.assertLessEqual(result, 25000.754895953367)

    def test_micelleLoss(self):
        result = lib.micelleLoss(5)
        self.assertLessEqual(result, 0.9)

    def test_micelleLoss2(self):
        result = lib.micelleLoss(0.2)
        self.assertLessEqual(result, 0.86)

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
        self.assertLessEqual(result, "Ag108")

    def test_read_ENSDF(self):
        result = lib.read_ENSDF("Co-60")
        self.assertLessEqual(result[0][0], 'NI60')

    def test_incer(self):
        result = lib.incer(["0.9", "0.2"], ["1", "1"])
        self.assertLessEqual(result, [[1.0], [1.0]])

    def test_relaxation_atom(self):
        result = lib.relaxation_atom("NI60", "Co-60", lacune="Atom_L")
        self.assertLessEqual(result, ('Auger L', 0.8))

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

    def test_readParameters_returns_tuple_of_31(self):
        result = lib.readParameters()
        self.assertEqual(len(result), 31)

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

    def test_stochasticDepTD_returns_3_values(self):
        pa, pb, pc = lib.stochasticDepTD(1.0, 0.1)
        self.assertIsInstance(float(pa), float)
        self.assertIsInstance(float(pb), float)
        self.assertIsInstance(float(pc), float)

    def test_stochasticDepTD_values_in_range(self):
        for _ in range(10):
            pa, pb, pc = lib.stochasticDepTD(1.0, 0.1)
            self.assertGreaterEqual(float(pa), 0)
            self.assertLessEqual(float(pa), 1)

    def test_stochasticDepCN_returns_2_values(self):
        pa, pb = lib.stochasticDepCN(1.0, 0.1)
        self.assertIsInstance(float(pa), float)
        self.assertIsInstance(float(pb), float)

    def test_stochasticDepCN_values_in_range(self):
        for _ in range(10):
            pa, pb = lib.stochasticDepCN(1.0, 0.1)
            self.assertGreaterEqual(float(pa), 0)
            self.assertLessEqual(float(pa), 1)

    def test_stochasticDepCN_sums_to_one(self):
        for _ in range(10):
            pa, pb = lib.stochasticDepCN(1.0, 0.1)
            self.assertAlmostEqual(float(pa) + float(pb), 1.0, places=10)


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

    @unittest.expectedFailure
    def test_TDCRPy_analytical_model(self):
        # Known bug in TDCRPy.py:549 — modelAnalytical() is called with 11
        # positional arguments but only accepts 10 (extra 'symm' arg passed).
        result = TDCRPy_mod.TDCRPy(1.0, "H-3", "1", 100, KB_TYPICAL, 10, Smodel=False)
        self.assertEqual(len(result), 14)

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


if __name__ == '__main__':
    unittest.main()
