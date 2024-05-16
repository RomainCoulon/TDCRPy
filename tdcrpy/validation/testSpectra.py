# -*- coding: utf-8 -*-
"""
Created on Tue May  7 09:17:11 2024

@author: romain.coulon
"""

import tdcrpy as td
import numpy as np
import matplotlib.pyplot as plt
import cProfile

"""
TEST SAMPLING IN BETASHAPE SPECTRA
"""
rad="Co-60"
E0 = 96.16 # expected mean energy BetaShape 2.3.1 in keV
uE0 = 0.20 # standard uncertainty
print(td.TDCR_model_lib.readBetaShapeInfo(rad,"beta-",'tot'))
out=td.TDCR_model_lib.readBetaShape(rad,"beta-",'tot_myEstep') # regular bining of 0.1 keV
# out=td.TDCR_model_lib.readBetaShape("Co-60","beta-",'tot') # irregular bining
E1 = sum(np.asarray(out[0])*np.asarray(out[1]))
print(f"mean energy = {E0} +/ {uE0} keV")
print(f"E (analytical estimation) = {E1} keV")
print(f"error (analytical estimation) = {E1-E0} keV")
N = 10000
e_vec = []
for i in range(N):
    ind = td.TDCR_model_lib.sampling(out[1])
    e_vec.append(out[0][ind])    
EM = np.mean(e_vec)
print(f"E (Monte Carlo estimation) = {EM}")
print(f"error (Monte Carlo estimation) = {EM-E0}")

"""
TEST SAMPLING IN DEPOSITED ENERGY SPECTRA
"""
out=td.TDCR_model_lib.readBetaSpectra(rad) # regular bining of 0.1 keV
Ed = sum(np.asarray(out[0])*np.asarray(out[1]))
print(f"deposited E (analytical estimation) = {Ed} keV")
print(f"shift (analytical estimation) = {Ed-E1} keV")

"""
PROFILING
"""

def regularBin():
    out = td.TDCR_model_lib.readBetaShape(rad,"beta-",'tot_myEstep') # regular bining of 0.1 keV
    ind = td.TDCR_model_lib.sampling(out[1])
    pass

def irregularBin():
    out = td.TDCR_model_lib.readBetaShape("Co-60","beta-",'tot') # irregular bining
    ind = td.TDCR_model_lib.sampling(out[1])