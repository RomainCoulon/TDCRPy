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

Rad = "H-3"

AB = 657.296
BC = 695.919
AC = 693.579
T = 448.868

D = AB+BC+AC-2*T

# kB08 = 0.8E-5
kB09 = 0.9E-5
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

# def readEffs(L, pS):
        

    
# L08, pS08, upS08 = readEff(Rad, kB08)
L09, pS09, upS09 = readEff(Rad, kB09)
L10, pS10, upS10 = readEff(Rad, kB10)
L11, pS11, upS11 = readEff(Rad, kB11)
L12, pS12, upS12 = readEff(Rad, kB12)

# pS08_smo = sg.savgol_filter(pS08, 53, 3) # window size used for filtering, order of fitted polynomial
pS09_smo = sg.savgol_filter(pS09, 21, 3) # window size used for filtering, order of fitted polynomial
pS10_smo = sg.savgol_filter(pS10, 21, 3) # window size used for filtering, order of fitted polynomial
pS11_smo = sg.savgol_filter(pS11, 21, 3) # window size used for filtering, order of fitted polynomial    
pS12_smo = sg.savgol_filter(pS12, 21, 3) # window size used for filtering, order of fitted polynomial


def eff0(L):
    i = np.where(np.asarray(L09)-L>0)[0][0]
    pS = pS10_smo[i]
    pT  = pS**3
    pD  = 3*pS**2-2*pT
    return pT, pD

def res0(L, D, T):
    print(L)
    return (eff0(L)[0]/eff0(L)[1]-T/D)**2

def eff(L):
    iA = np.where(np.asarray(L09)-L[0]>0)[0][0]
    pSA = pS10_smo[iA]
    iB = np.where(np.asarray(L09)-L[1]>0)[0][0]
    pSB = pS10_smo[iB]
    iC = np.where(np.asarray(L09)-L[2]>0)[0][0]
    pSC = pS10_smo[iC]
    pAB = pSA*pSB
    pBC = pSB*pSC
    pAC = pSA*pSC
    pT  = pSA*pSB*pSC
    pD  = pAB+pBC+pAC-2*pT
    return pT, pD

def res(L, D, T):
    print(L)
    return (eff(L)[0]/eff(L)[1]-T/D)**2

L = 0.5
r0 = opt.minimize_scalar(res0, args=(D, T), method='golden')
L = (r0.x, r0.x, r0.x)
r = opt.minimize(res, L, args=(D, T), method='nelder-mead',options={'xatol': 1e-6, 'disp': True, 'maxiter':100})
L=r.x[:3]

print("Efficiency of double coincidences = ", eff(L)[1])
print("Efficiency of triple coincidences = ", eff(L)[0])
print("Activity from double coincidences = ", D/eff(L)[1], "Bq")
print("Efficiency of triple coincidences = ", T/eff(L)[0], "Bq")

plt.figure("Efficiency curve")
plt.clf()
plt.title('  ')
# plt.errorbar(L08, pS08, yerr = upS08, fmt=".m", label = "kB = 0.008 cm/MeV")
# plt.plot(L08,pS08_smo,'-m', label = "Savitzky-Golay Filter")
plt.errorbar(L09, pS09, yerr = upS09, fmt=".k", label = "kB = 0.009 cm/MeV")
plt.plot(L09,pS09_smo,'-k', label = "Savitzky-Golay Filter")
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


def I2calc(td,tab,tbc,tac,rad,kB,*,f=0,xatol=1e-5):
    """This function aims to calculate the activity of an LS source by the TDCR primary standardization technique.
    This function is more especially used to calculate the key comparison indicator I2 for international comparison of standard solutions containing pure beta or alpha radionuclides.
    
    Reference:
        [Romain Coulon et al 2020 Metrologia 57 035009]
        https://iopscience.iop.org/article/10.1088/1681-7575/ab7e7b
            
    :param td: Triple to double coincidence count rate
    :type td: float
    :param tab: Triple coincidence count rate divided to the coincidence count rate between channels A and B
    :type tab: float
    :param tbc: Triple coincidence count rate divided to the coincidence count rate between channels B and C
    :type tbc: float 
    :param tac: Triple coincidence count rate divided to the coincidence count rate between channels A and C
    :type tac: float 
    :param rad: Name of the radionuclide (eg. C-14)
    :type rad: string     
    :param kB: Birks constant /(cm2.MeV^{-1})
    :type kB: float  
    :param f: (Optional) Probability of non-detection of the single photoelectron pulse (By default: set equal to 0)
    :type f: float
    :param xatol: (Optional) Value of the convergence paramter of the minimization procedure (By default: set equal to 1e-5)
    :type xatol: float    
    :param L0: Value of the global free parameter evaluated using the TDCR value /(phe.eV^{-1})
    :type L0: float
    :param L: Values of the free parameters of individual channels evaluated using the T/AB, T/BC, T/AC values /(phe.eV^{-1})
    :type L: array of floats    
    :param effD0: Efficiency of the double coincidence count rate evaluated using the TDCR value
    :type effD0: float
    :param effD: Efficiency of the double coincidence count rate evaluated using the T/AB, T/BC, T/AC values
    :type effD: float 
            
    :return: L0, L, effD0, effD
    :rtype: tuple
    """
    fin=f
    # Estimation of the quanched function Q    
    e,p,Q,em=QuencFunc(rad,kB)

    # Estimation of the global free parameter using the TDCR value
    r0=minimize_scalar(tdcrmodel0, args=(p,Q,em,td,0,fin), method='golden')
    L0=r0.x # Value of the free parameter
    # Estimation of the double coincidence efficiency
    effD0=tdcrmodel0(L0,p,Q,em,td,1,fin)

    # approximative prior knoleged of the free paramters to accelerate the calculation 
    pAd=1.0
    pBd=1.0
    pCd=1.0
    
    L=empty(3) # initialization of the individual free parameters / (phe.eV^{-1})
    L[0]=L0*pAd
    L[1]=L0*pBd
    L[2]=L0*pCd
    
    # Estimation of the global free parameter using the T/AB, T/BC, T/AC
    r=minimize(tdcrmodel, L, args=(p,Q,em,tab,tbc,tac,0,fin), method='nelder-mead',options={'xatol': xatol, 'disp': False, 'maxiter':100})
    L=r.x[:3]
    # Estimation of the double coincidence efficiency
    effD=tdcrmodel(L,p,Q,em,tab,tbc,tac,1,fin)
    return L0, L, effD0, effD 





    
