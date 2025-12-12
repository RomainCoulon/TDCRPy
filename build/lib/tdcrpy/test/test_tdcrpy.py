# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:04:46 2023

test_my_module.py

@author: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

import unittest
import importlib.resources
import numpy as np
import tdcrpy


with importlib.resources.path('tdcrpy', 'MCNP-MATRIX') as data_path:
    fp1 = data_path / 'matrice/fichier/matrice_10ml-photon_1_200k.txt'

class TestMyModule(unittest.TestCase):

    def test_normalise(self):
        result = tdcrpy.TDCR_model_lib.normalise([0.8,0.2,0.1])
        self.assertEqual(result, [0.7272727272727273, 0.18181818181818182, 0.09090909090909091])

    def test_sampling(self):
        result = tdcrpy.TDCR_model_lib.sampling([0.1,0.5,0.2,0.2])
        self.assertIn(result, [0,1,2,3])
        
    def test_readPenNuc2(self):
        result = tdcrpy.TDCR_model_lib.readPenNuc2("H-3")
        self.assertEqual(result[0], ["HE3"])
        self.assertEqual(result[1], [1.0])
        self.assertEqual(result[2], [18.591])
        
    def test_stoppingpowerA(self):
        result = tdcrpy.TDCR_model_lib.stoppingpowerA(5000)
        self.assertEqual(result, 871008.0)
        
    def test_stoppingpower(self):
        result = tdcrpy.TDCR_model_lib.stoppingpower(100)
        self.assertEqual(result, 326.78)
        
    # def test_readBetaShape(self):
    #     result = np.mean(tdcrpy.TDCR_model_lib.readBetaShape("H-3","beta-","tot"))
    #     self.assertEqual(result, 4.661922025723472)
    
    def test_readBetaSpectra(self):
        result = np.mean(tdcrpy.TDCR_model_lib.readBetaSpectra("H-3"))
        self.assertEqual(result, 4.651594887459807)
    
    def test_E_quench_e(self):
        result = np.mean(tdcrpy.TDCR_model_lib.E_quench_e(100*1e3,100*1e3,0.01,20000)*1e-3)
        self.assertEqual(result, 91.24694728992633)
        
    def test_E_quench_e2(self):
        result = np.mean(tdcrpy.TDCR_model_lib.E_quench_e(100*1e3,50*1e3,0.01,20000)*1e-3)
        self.assertEqual(result, 47.97271731307207)
        
    def test_E_quench_a(self):
        result = np.mean(tdcrpy.TDCR_model_lib.E_quench_a(5000,1e-5,20000))
        self.assertEqual(result, 344.28633759073443)
    
    # run_interpolate tested inside Em_a() and Em_e()
    
    def test_Em_a(self):
        result = tdcrpy.TDCR_model_lib.Em_a(5000, 1e-5, 10000)
        self.assertLessEqual(result, 344.12105896)
    
    def test_Em_e(self):
        result = tdcrpy.TDCR_model_lib.Em_e(50000, 50000, 1e-5, 10000)
        self.assertLessEqual(result, 48568.14392205286)
        
    def test_Em_e2(self):
        result = tdcrpy.TDCR_model_lib.Em_e(50000, 25000, 1e-5, 10000)
        self.assertLessEqual(result, 25000.754895953367)        
        
    def test_micelleLoss(self):
        result = tdcrpy.TDCR_model_lib.micelleLoss(5)
        self.assertLessEqual(result, 0.9)             

    def test_micelleLoss2(self):
        result = tdcrpy.TDCR_model_lib.micelleLoss(0.2)
        self.assertLessEqual(result, 0.86)         
        
    def test_read_matrice(self):
        result = tdcrpy.TDCR_model_lib.read_matrice(fp1,0)[50][50]
        self.assertEqual(result, 0.000718)
    
    # energie_dep_gamma() no more used.
    
    def test_energie_dep_gamma(self):
        result = tdcrpy.TDCR_model_lib.energie_dep_gamma2(100,10)
        self.assertLessEqual(result, 100)
    
    # energie_dep_beta() no more used.
    
    def test_energie_dep_beta(self):
        result = tdcrpy.TDCR_model_lib.energie_dep_beta2(100,10)
        self.assertLessEqual(result, 100)    
        
    # writeEffcurves() no more used.
    
    def test_transf_name(self):
        result = tdcrpy.TDCR_model_lib.transf_name("108Ag")
        self.assertLessEqual(result, "Ag108")

    def test_read_ENSDF(self):
        result = tdcrpy.TDCR_model_lib.read_ENSDF("Co-60")
        self.assertLessEqual(result[0][0], 'NI60')
    
    def test_incer(self):
        result = tdcrpy.TDCR_model_lib.incer(["0.9", "0.2"],["1", "1"])
        self.assertLessEqual(result, [[1.0], [1.0]])
        
    def test_relaxation_atom(self):
        result = tdcrpy.TDCR_model_lib.relaxation_atom("NI60","Co-60",lacune="Atom_L")
        self.assertLessEqual(result, ('Auger L', 0.8))
    
    def test_format_modif(self):
        result = tdcrpy.TDCR_model_lib.format_modif("17-3")
        self.assertLessEqual(result, '17E-3')
        
    def test_reperer_energie_index(self):
        result = tdcrpy.TDCR_model_lib.reperer_energie_index(17,[5,11,13,17,52,64])
        self.assertLessEqual(result, 3)
    
    # modelAnalytical(L,TD,TAB,TBC,TAC,rad,kB,V,mode,mode2,ne)
    
    # clear_terminal()
    
    # display_header()
    
if __name__ == '__main__':
    unittest.main()