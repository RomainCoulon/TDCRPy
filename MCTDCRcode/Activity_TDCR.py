# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 20:43:26 2023

@author: romain.coulon
"""
import matplotlib.pyplot as plt
import TDCR_model_lib as tl
import numpy as np
import scipy.signal as sg

Rad = "H-3"

AB = 0
BC = 0
AC = 0
T = 0

# kB08 = 0.8E-5
# kB09 = 0.9E-5
kB10 = 1E-5
kB11 = 1.1E-5
kB12 = 1.2E-5

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
    
# L08, pS08, upS08 = readEff(Rad, kB08)
# L09, pS09, upS09 = readEff(Rad, kB09)
L10, pS10, upS10 = readEff(Rad, kB10)
L11, pS11, upS11 = readEff(Rad, kB11)
L12, pS12, upS12 = readEff(Rad, kB12)

# pS08_smo = sg.savgol_filter(pS08, 53, 3) # window size used for filtering, order of fitted polynomial
# pS09_smo = sg.savgol_filter(pS09, 53, 3) # window size used for filtering, order of fitted polynomial
pS10_smo = sg.savgol_filter(pS10, 21, 3) # window size used for filtering, order of fitted polynomial
pS11_smo = sg.savgol_filter(pS11, 21, 3) # window size used for filtering, order of fitted polynomial    
pS12_smo = sg.savgol_filter(pS12, 21, 3) # window size used for filtering, order of fitted polynomial

plt.figure("Efficiency curve")
plt.clf()
plt.title('  ')
# plt.errorbar(L08, pS08, yerr = upS08, fmt=".m", label = "kB = 0.008 cm/MeV")
# plt.plot(L08,pS08_smo,'-m', label = "Savitzky-Golay Filter")
# plt.errorbar(L09, pS09, yerr = upS09, fmt=".k", label = "kB = 0.009 cm/MeV")
# plt.plot(L09,pS09_smo,'-k', label = "Savitzky-Golay Filter")
plt.errorbar(L10, pS10, yerr = upS10, fmt=".b", label = "kB = 0.010 cm/MeV")
plt.plot(L10,pS10_smo,'-b', label = "Savitzky-Golay Filter")
plt.errorbar(L11, pS11, yerr = upS11, fmt=".r", label = "kB = 0.011 cm/MeV")
plt.plot(L11,pS11_smo,'-r', label = "Savitzky-Golay Filter")
plt.errorbar(L12, pS12, yerr = upS12, fmt=".y", label = "kB = 0.012 cm/MeV")
plt.plot(L12,pS12_smo,'-y', label = "Savitzky-Golay Filter")
# plt.xscale("log")
plt.xlabel(r"$L$ /keV$^{-1}$", fontsize = 14)
plt.ylabel(r"$\epsilon$", fontsize = 14)
plt.legend(fontsize = 12)
# plt.close()
    






    
