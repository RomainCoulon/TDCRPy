# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 20:43:26 2023

@author: romain.coulon
"""
import matplotlib.pyplot as plt
import TDCR_model_lib as tl
import numpy as np
import scipy.signal as sg
import scipy.optimize as opt

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

def eff0(L, Lv, pSv):
    pS = np.interp(L, Lv, pSv)
    pT  = pS**3
    pD  = 3*pS**2-2*pT
    return pT, pD

def res0(L, TD, Lv, pSv):
    return (eff0(L, Lv, pSv)[0]/eff0(L, Lv, pSv)[1]-T/D)**2

def eff(L, Lv, pSv):
    pSA = np.interp(L[0], Lv, pSv)
    pSB = np.interp(L[1], Lv, pSv)
    pSC = np.interp(L[2], Lv, pSv)
    pAB = pSA*pSB
    pBC = pSB*pSC
    pAC = pSA*pSC
    pT  = pSA*pSB*pSC
    pD  = pAB+pBC+pAC-2*pT
    return pT, pAB, pBC, pAC, pD

def res(L, TAB, TBC, TAC, Lv, pSv):
    a = eff(L, Lv, pSv)
    return (a[0]/a[1]-TAB)**2+(a[0]/a[2]-TBC)**2+(a[0]/a[3]-TAB)**2

def readProfil(Rad,kB):
    Lv, pSv, upSv = readEff(Rad, kB)
    pSv = sg.savgol_filter(pSv, 21, 3)
    return Lv, np.asarray(pSv), upSv

def effTDCR(TD, TAB, TBC, TAC, Rad, kB):
    Lv, pSv, _ = readProfil(Rad,kB)
    L = 0.5
    r0 = opt.minimize_scalar(res0, args=(TD, Lv, pSv), method='golden')
    L0 = (r0.x, r0.x, r0.x)
    r = opt.minimize(res, L0, args=(TAB, TBC, TAC, Lv, pSv), method='nelder-mead',options={'xatol': 1e-7, 'disp': True, 'maxiter':400})
    L=r.x[:3]
    return L0[0], L, eff0(L0[0], Lv, pSv)[1], eff(L, Lv, pSv)[4], eff0(L0[0], Lv, pSv)[0], eff(L, Lv, pSv)[0]

def plotEffProfil(Rad,kB):
    L = readProfil(Rad,kB)[0]
    effS = readProfil(Rad,kB)[1]
    effT = effS**3
    effD = 3*effS**2-2*effT
    tdcr = effT/effD
    
    plt.figure("Efficiency curve eff = f(L)")
    plt.clf()
    plt.title('  ')
    plt.plot(L, effS,'ob', label = r"$\epsilon_S$, $kB$ = "+str(kB)+" cm/keV")
    plt.plot(L, effD,'ok', label = r"$\epsilon_D$, $kB$ = "+str(kB)+" cm/keV")
    plt.plot(L, effT,'or', label = r"$\epsilon_T$, $kB$ = "+str(kB)+" cm/keV")
    plt.xscale("log")
    plt.xlabel(r"$L$ /keV$^{-1}$", fontsize = 14)
    plt.ylabel(r"$\epsilon$", fontsize = 14)
    plt.legend(fontsize = 12)
    # plt.close()
    
    plt.figure("Efficiency curve eff = f(TDCR)")
    plt.clf()
    plt.title('  ')
    plt.plot(tdcr, effS,'ob', label = r"$\epsilon_S$, $kB$ = "+str(kB)+" cm/keV")
    plt.plot(tdcr, effD,'ok', label = r"$\epsilon_D$, $kB$ = "+str(kB)+" cm/keV")
    plt.plot(tdcr, effT,'or', label = r"$\epsilon_T$, $kB$ = "+str(kB)+" cm/keV")
    # plt.xscale("log")
    plt.xlabel(r"$\epsilon_T/\epsilon_D$", fontsize = 14)
    plt.ylabel(r"$\epsilon$", fontsize = 14)
    plt.legend(fontsize = 12)
    # plt.close()

AB = 657.296
BC = 695.919
AC = 693.579
T = 448.868
D = AB+BC+AC-2*T
out = effTDCR(T/D, T/AB, T/BC, T/AC, "H-3", 1E-5) # H-3
# out = effTDCR(657.296, 695.919, 693.579, 448.868) # Co-60

digRound = 2
print("Efficiency of double coincidences = ", round(100*out[2], digRound),"%")
print("Efficiency of triple coincidences = ", round(100*out[4], digRound),"%")
print("Activity from double coincidences = ", round(D/out[2], digRound), "Bq")
print("Activity of triple coincidences = ", round(T/out[4], digRound), "Bq")
print("Efficiency of double coincidences (asym) = ", round(100*out[3], digRound),"%")
print("Efficiency of triple coincidences (asym) = ", round(100*out[5], digRound),"%")
print("Activity from double coincidences (asym) = ", round(D/out[3], digRound), "Bq")
print("Activity of triple coincidences (asym) = ", round(T/out[5], digRound), "Bq")

plotEffProfil("H-3",1E-5)