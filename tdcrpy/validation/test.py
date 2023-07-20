# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 19:15:23 2023
Validation test of the TDCRPy code
@author: romain.coulon
"""
import tdcrpy as td


"""
Validation with standard solution for Co-60
"""

# Rad="Co-60"    # list of radionuclides (Na-24)
# pmf_1="1"
# kB =1.0e-5       # Birks constant in cm/keV
# RHO = 0.98
# nE = 1000
# TD = 0.977667386529166        # Measured TDCR value
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# L = 1
# N = 10

# resuts_1=td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, RHO, nE, "eff", "sym", Display=False)
# print(resuts_1)

# resuts_2=td.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, "sym", N=N)
# print(resuts_2)


"""
Validation with the analytical model
"""

Rad='H-3'
pmf_1="1"
kB =1.0e-5       # Birks constant in cm/keV
RHO = 0.98
nE = 1000
TD = 0.70        # Measured TDCR value
TAB = 0.992232838598821
TBC = 0.992343419459002
TAC = 0.99275350064608
L = 1
N = 5000

resuts_1=td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, RHO, nE, "eff", "sym", Display=False)
print(resuts_1)

resuts_2=td.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, "sym", N=N)
print(resuts_2)
