# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 10:04:53 2023

@author: romain.coulon
"""


import TDCRPy as td
import scipy.optimize as opt


## INPUT OF THE MODEL
N=10                   # number of simulated decay (MC trials)
Rad="Co-60"    # list of radionuclides (Na-24)
pmf_1="1"        # relative abondance (pmf)
kB =1.0e-5       # Birks constant in cm/keV  
L=(1.13, 1.13, 1.13)           # Free paramete in keV-1


TDCR_measure = 0.977784        # Measured TDCR value
u_TDCR_measure = 0.000711      # standard uncertainty

TAB = 0.85
TBC = 0.85
TAC = 0.85

Display = False               # to display calculation results on the console
# RHO = 0.96         #density of absorber (Toluene) g/cm3
RHO = 0.98           #density of absorber (UG + H20) g/cm3
nE = 1000            #number of bin to discretize the energy vector for scintillation quenching calculation

mode = "res"
mode2 = "asym"

r=opt.minimize(td.TDCRPy, L, args=(TDCR_measure, TAB, TBC, TAC, Rad, pmf_1, N, kB, RHO, nE, mode, mode2), method='nelder-mead',options={'xatol': 1e-7, 'disp': True, 'maxiter':100})
L=r.x
print(r)
print(L)

mode = "eff"
out=td.TDCRPy(L,TDCR_measure, TAB, TBC, TAC, Rad, pmf_1, N, kB, RHO, nE, mode, mode2)
print(out)