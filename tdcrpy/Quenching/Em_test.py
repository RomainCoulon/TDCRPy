# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 11:47:21 2023

@author: romain.coulon
"""

import matplotlib.pyplot as plt
import tdcrpy
import numpy as np
from tqdm import tqdm
import scipy.interpolate as  interp


kB_a = [6e-6, 7e-6, 8e-6, 9e-6, 1e-5, 1.1e-5] #, 1.2e-5, 1.3e-5, 1.4e-5, 1.5e-5] # cm/MeV

Ei_alpha_fid = open("inputVecteurAlpha.txt")
Ei_alpha = Ei_alpha_fid.readlines()
Ei_alpha = Ei_alpha[0].split(" ")
Ei_alpha = [float(x) for x in Ei_alpha[:-1]]

Em_alpha = []
for ikB in kB_a:
    fid = open("QuenchEnergyAlpha_"+str(ikB)+".txt")
    line = fid.readlines()
    line = line[0].split(" ")
    line = [float(x) for x in line[:-1]]
    Em_alpha.append(line)
    
Ei_electron_fid = open("inputVecteurElectron.txt")
Ei_electron = Ei_electron_fid.readlines()
Ei_electron = Ei_electron[0].split(" ")
Ei_electron = [float(x) for x in Ei_electron[:-1]]

kB_e = [0.006, 0.007, 0.008, 0.009, 0.010, 0.011, 0.012, 0.013, 0.014, 0.015] # cm/MeV
Em_electron = []
for ikB in kB_e:
    fid = open("QuenchEnergyElectron_"+str(ikB)+".txt")
    line = fid.readlines()
    line = line[0].split(" ")
    line = [float(x) for x in line[:-1]]
    Em_electron.append(line)

def Em(E,Ev,Emv,kB,p):
    """
    This fonction performs a cubic splin interpolation of pre-calculated quenching function using E_quench_e() and E_quench_a() functions.
    It aims to gain calculation time while inducing an acceptable calculation error.

    Parameters
    ----------
    E : float
        Input energy in eV for electron and in keV for alpha
    Ev : list
        Vector of deposited energies eV for electron and in keV for alpha
    Emv : list
        Vector of quenched energies eV for electron and in keV for alpha
    kB : float
        Birks constant in cm/MeV for electron and in cm/keV for alpha
    p : string
        Particule, electron or alpha

    Returns
    -------
    Float
        interpolated quenched energy in eV for electron and in keV for alpha

    """
    if p == "electron":
        # setup for electron
        Et = 5000 # threshold to interpolated values above 5 keV
        kB_vec = [0.006, 0.007, 0.008, 0.009, 0.010, 0.011, 0.012, 0.013, 0.014, 0.015] # in cm/MeV
    elif p == "alpha":
        # setup for alpha
        Et = 5 # threshold to interpolated values above 5 keV
        kB_vec = [6e-6, 7e-6, 8e-6, 9e-6, 1e-5, 1.1e-5, 1.2e-5, 1.3e-5, 1.4e-5, 1.5e-5] # in cm/keV
    
    if E <= Et:
        # run the accurate quenching model
        if p == "electron": r = tdcrpy.TDCR_model_lib.E_quench_e(E,kB,int(1e4))
        if p == "alpha": r = tdcrpy.TDCR_model_lib.E_quench_a(E,kB,int(1e4))
    else:
        # run interpolation
        if kB in kB_vec:
            # Exact value for the kB
            kBin = True
            ind_k = kB_vec.index(kB)
        else:
            # non exact value for the kB
            # find the index just above the true value
            for index, value in enumerate(kB_vec):
                ind_k = -1
                if value > kB:
                    ind_k = index
                    break        
            kBin = False
        
        for index, value in enumerate(Ev):
            # find the index just above the location of the exact value of the input energy 
            ind = -1
            if value > E:
                ind = index
                break

        m = 5 # set the window depht of the spline interpolation around the energy index
        if kBin:
            # case of exact kB value
            if ind<m and ind != -1:
                # troncated window on low values
                f = interp.UnivariateSpline(Ev[ind:ind+m], Emv[ind_k][ind:ind+m])
            elif ind>len(Ev)-m or ind==-1:
                # troncated window on high values
                f = interp.UnivariateSpline(Ev[ind-m:ind], Emv[ind_k][ind-m:ind])
            else:
                f = interp.UnivariateSpline(Ev[ind-m:ind+m], Emv[ind_k][ind-m:ind+m])
            r = f(E) # interpolated quenched energy
        else:
            # interpolation for the two indexes surounding the exact energy
            if ind<m and ind != -1:
                # troncated window on low values
                f1 = interp.UnivariateSpline(Ev[ind:ind+m], Emv[ind_k][ind:ind+m])
                f2 = interp.UnivariateSpline(Ev[ind:ind+m], Emv[ind_k-1][ind:ind+m])
            elif ind>len(Ev)-m or ind==-1:
                # troncated window on high values
                f1 = interp.UnivariateSpline(Ev[ind-m:ind], Emv[ind_k][ind-m:ind])
                f2 = interp.UnivariateSpline(Ev[ind-m:ind], Emv[ind_k-1][ind-m:ind])
            else:
                f1 = interp.UnivariateSpline(Ev[ind-m:ind+m], Emv[ind_k][ind-m:ind+m])
                f2 = interp.UnivariateSpline(Ev[ind-m:ind+m], Emv[ind_k-1][ind-m:ind+m])
            # linear interpolation for the estimation related to the exact kB value
            r = f2(E)+(f1(E) - f2(E))/(kB_vec[ind_k]-kB_vec[ind_k-1])*(kB-kB_vec[ind_k-1])
            
    return r



Err_e = []
Err_a = []
kBtest = 0.0104

Ein_ve = np.logspace(2, 7, 100) # in eV
for Ein in tqdm(Ein_ve, desc="Processing", unit=" inter"):
    Eout = Em(Ein,Ei_electron,Em_electron, kBtest, "electron")
    Eref = tdcrpy.TDCR_model_lib.E_quench_e(Ein,kBtest,int(1e4))
    Err_e.append(100*abs(Eout/Eref-1))

Ein_va = np.logspace(-1, 4, 100) # in keV
for Ein in tqdm(Ein_va, desc="Processing", unit=" inter"):    
    Eout = Em(Ein,Ei_alpha,Em_alpha, kBtest*1e-3, "alpha")
    Eref = tdcrpy.TDCR_model_lib.E_quench_a(Ein,kBtest*1e-3,int(1e4))
    # print("\n",Ein, Eout, Eref)
    Err_a.append(100*abs(Eout/Eref-1))

plt.figure('Error')
plt.clf()
plt.plot(np.asarray(Ein_va),Err_a,label=r'Err for alpha',linewidth=3, color='k')
plt.plot(np.asarray(Ein_ve*1e-3),Err_e,label=r'Err for electrons',linewidth=3, color='r')
plt.yscale("log")
plt.xscale("log")
plt.ylabel(r'Error (%)', fontsize=18)
plt.xlabel(r'$E_{i}$ /keV', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),ncol=2,fontsize=14)

plt.figure('quenching function for electrons')
plt.clf()
plt.title("Quenched energy of electrons",fontsize=16)
for i, ikB in enumerate(kB_e):
    plt.plot(np.asarray(Ei_electron),np.asarray(Em_electron[i]),label=r'$kB$ = '+str(round(ikB,3))+'cm/MeV' ,linewidth=3)
plt.ylabel(r'Em($E_{i}$) /keV', fontsize=18)
plt.xlabel(r'$E_{i}$ /keV', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.9),ncol=2,fontsize=14)

plt.figure('quenching function for alpha particle')
plt.clf()
plt.title("Quenched energy of alpha particles",fontsize=16)
for i, ikB in enumerate(kB_a[:6]):
    plt.plot(np.asarray(Ei_alpha),np.asarray(Em_alpha[i]),label=r'$kB$ = '+str(round(ikB*1e3,3))+'cm/MeV' ,linewidth=3)
plt.ylabel(r'Em($E_{i}$) /keV', fontsize=18)
plt.xlabel(r'$E_{i}$ /keV', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.9),ncol=2,fontsize=14)
    