# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 20:43:26 2023

@author: romain.coulon
"""
import matplotlib.pyplot as plt
import TDCR_model_lib as tl
import numpy as np
import scipy.signal as sg

Rad = "Co-60"

kB10 = 1E-5
kB11 = 1.1E-5
AB = 0
BC = 0
AC = 0
T = 0

def readEff(Rad, kB):
    file = open("EfficiencyCurves/"+''.join(Rad)+"/Eff_"+''.join(Rad)+'_[1]_'+str(kB)+".txt","r")
    L=[]
    pS=[]
    upS=[]
    for row in file:
        H=row.split(" ")
        L.append(float(H[0]))
        pS.append(float(H[1]))
        upS.append(float(H[2]))
    return L, pS, upS
    

L10, pS10, upS10 = readEff(Rad, kB10)
L11, pS11, upS11 = readEff(Rad, kB11)

pS10_smo = sg.savgol_filter(pS10, 53, 3) # window size used for filtering, order of fitted polynomial
pS11_smo = sg.savgol_filter(pS11, 53, 3) # window size used for filtering, order of fitted polynomial    

plt.figure("Efficiency curve")
plt.clf()
plt.title('  ')
plt.errorbar(L10, pS10, yerr = upS10, fmt=".b", label = "kB = 0.010 cm/MeV")
plt.plot(L10,pS10_smo,'-b', label = "Savitzky-Golay Filter")
plt.errorbar(L11, pS11, yerr = upS11, fmt=".r", label = "kB = 0.011 cm/MeV")
plt.plot(L11,pS11_smo,'-r', label = "Savitzky-Golay Filter")
# plt.xscale("log")
plt.xlabel(r"$L$ /keV$^{-1}$", fontsize = 14)
plt.ylabel(r"$\epsilon$", fontsize = 14)
plt.legend(fontsize = 12)
# plt.close()
    

    
