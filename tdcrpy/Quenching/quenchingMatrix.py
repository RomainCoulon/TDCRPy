# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 13:45:13 2023

This programme aims to produce the quenching matrix from E_quench_e() and E_quench_a().

@author: romain.coulon
"""
import numpy as np
import tdcrpy
from tqdm import tqdm

def E_quench_e_M(kB): # in cm/MeV
    E = np.linspace(0, int(1e7), int(1e4)) # in eV
    with open("QuenchEnergyElectron_"+str(kB)+".txt", "w") as file2:
        for Ei in tqdm(E, desc="Processing", unit=" inter"):
            Eq = tdcrpy.TDCR_model_lib.E_quench_e(Ei,kB,int(1e4)) # in eV
            file2.write(str(Eq)+" ")
    return E.tolist(), Eq

E_quench_e_M(0.006)
E_quench_e_M(0.007)
E_quench_e_M(0.008)
E_quench_e_M(0.009)
E_quench_e_M(0.010)
E_quench_e_M(0.011)
E_quench_e_M(0.012)
E_quench_e_M(0.013)
E_quench_e_M(0.014)
E = E_quench_e_M(0.015)[0]

with open("inputVecteurElectron.txt", "w") as file:
    for Ei in E:
        file.write(str(Ei)+" ")

def E_quench_a_M(kB): # in cm/keV
    E = np.linspace(0, int(1e4), int(1e4)) # in keV
    with open("QuenchEnergyAlpha_"+str(kB*1e-3)+".txt", "w") as file2:
        for Ei in tqdm(E, desc="Processing", unit=" inter"):
            Eq = tdcrpy.TDCR_model_lib.E_quench_a(Ei,kB*1e-3,int(1e4)) # in eV
            file2.write(str(Eq)+" ")
    return E.tolist(), Eq

E_quench_a_M(0.006)
E_quench_a_M(0.007)
E_quench_a_M(0.008)
E_quench_a_M(0.009)
E_quench_a_M(0.010)
E_quench_a_M(0.011)
E_quench_a_M(0.012)
E_quench_a_M(0.013)
E_quench_a_M(0.014)
E = E_quench_a_M(0.015)[0]

with open("inputVecteurAlpha.txt", "w") as file:
    for Ei in E:
        file.write(str(Ei)+" ")
    