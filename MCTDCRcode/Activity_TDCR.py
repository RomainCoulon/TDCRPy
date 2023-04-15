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

kB = 1E-5

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
    

L, pS, upS = readEff(Rad, kB)

pS_smo = sg.savgol_filter(pS, 53, 3) # window size used for filtering, order of fitted polynomial
    
plt.figure("Efficiency curve")
plt.clf()
plt.title('  ')
plt.errorbar(L, pS, yerr = upS, fmt=".b", label = "S")
plt.plot(L,pS_smo,'-r', label = "Savitzky-Golay Filter")
# plt.xscale("log")
plt.xlabel(r"$L$ /keV$^{-1}$", fontsize = 14)
plt.ylabel(r"$\epsilon$", fontsize = 14)
plt.legend(fontsize = 12)
# plt.close()
    

    
