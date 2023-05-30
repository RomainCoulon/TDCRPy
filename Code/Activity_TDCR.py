# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 20:43:26 2023

@author: romain.coulon
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sg
import scipy.optimize as opt

def readEff(Rad, kB, SDT):
    file = open("G:\Python_modules\TDCRPy\TDCRPy\Code\\EfficiencyCurves/"+''.join(Rad)+"/Eff"+SDT+"_"+''.join(Rad)+'_[1]_'+str(kB)+".txt","r")
    L=[]
    p=[]
    up=[]
    for row in file:
        H=row.split(" ")
        L.append(float(H[0]))
        p.append(float(H[1]))
        up.append(float(H[2]))
    return L, p, up

def readProfil(Rad,kB,SDT):
    Lv, pSv, upSv = readEff(Rad, kB, SDT)
    pSv = sg.savgol_filter(pSv, 11, 3)
    return Lv, np.asarray(pSv), upSv

def res(L, TD, Lv, pDv, pTv):
    return (np.interp(L, Lv, pTv)/np.interp(L, Lv, pDv)-TD)**2

def effTDCR(TD, Rad, kB):
    kBv = [0.8e-5, 0.9e-5, 1.0e-5, 1.1e-5, 1.2e-5]
    Lv0 = readProfil(Rad,0.8e-5,"D")[0]
    pDv=[]; pTv=[]; Lv1=[]; effD1=[]; effT1=[]
    for i in kBv:
        pDv.append(readProfil(Rad,i,"D")[1])
        pTv.append(readProfil(Rad,i,"T")[1])
        Lv1.append(opt.minimize_scalar(res, args=(TD, Lv0, pDv[-1], pTv[-1]), method='golden').x)
        effD1.append(np.interp(Lv1[-1], Lv0, pDv[-1]))
        effT1.append(np.interp(Lv1[-1], Lv0, pTv[-1]))
    L = np.interp(kB, kBv, Lv1)
    effD = np.interp(kB, kBv, effD1)
    effT = np.interp(kB, kBv, effT1)
     
    # plt.figure("eff_D = f(kB)")
    # plt.clf()
    # plt.title('  ')
    # plt.plot(kBv, effD1,'-b')
    # plt.plot(kB, effD, 'or')
    # plt.xlabel(r"$kB$ /(cm keV$^{-1}$)", fontsize = 14)
    # plt.ylabel(r"$\epsilon_D$", fontsize = 14)
    # plt.legend(fontsize = 12)
    # # plt.close()
    
    # plt.figure("eff_T = f(kB)")
    # plt.clf()
    # plt.title('  ')
    # plt.plot(kBv, effT1,'-b')
    # plt.plot(kB, effT, 'or')
    # plt.xlabel(r"$kB$ /(cm keV$^{-1}$)", fontsize = 14)
    # plt.ylabel(r"$\epsilon_T$", fontsize = 14)
    # plt.legend(fontsize = 12)
    # # plt.close()
    
    # plt.figure("L = f(kB)")
    # plt.clf()
    # plt.title('  ')
    # plt.plot(kBv, Lv1,'-b')
    # plt.plot(kB, L, 'or')
    # plt.xlabel(r"$kB$ /(cm keV$^{-1}$)", fontsize = 14)
    # plt.ylabel(r"$L$ /keV$^{-1}$", fontsize = 14)
    # plt.legend(fontsize = 12)
    # # plt.close()
    
    return L, effD, effT

def plotEffProfil(Rad,kB):
    kBv = [0.8e-5, 0.9e-5, 1.0e-5, 1.1e-5, 1.2e-5]
    Lv0 = readProfil(Rad,0.8e-5,"D")[0]
    pDv=[]; pTv=[]; pSv=[]
    for i in kBv:
        pDv.append(readProfil(Rad,i,"D")[1])
        pTv.append(readProfil(Rad,i,"T")[1])
        pSv.append(readProfil(Rad,i,"S")[1])
    effS=[]; effD=[]; effT=[]
    print("S",np.asarray(pSv)[:,0])
    for i in range(len(Lv0)):
        effS.append(np.interp(kB, kBv, np.asarray(pSv)[:,i]))
        effD.append(np.interp(kB, kBv, np.asarray(pDv)[:,i]))
        effT.append(np.interp(kB, kBv, np.asarray(pTv)[:,i]))
    tdcr = np.asarray(effT)/np.asarray(effD)
    
    plt.figure("Efficiency curve eff = f(L)")
    plt.clf()
    plt.title('  ')
    plt.plot(Lv0, effS,'-b', label = r"$\epsilon_{S}$, $kB$ = "+str(kB)+" cm/keV")
    plt.plot(Lv0, effD,'-k', label = r"$\epsilon_{D}$, $kB$ = "+str(kB)+" cm/keV")
    plt.plot(Lv0, effT,'-r', label = r"$\epsilon_{T}$, $kB$ = "+str(kB)+" cm/keV")
    plt.xscale("log")
    plt.xlabel(r"$L$ /keV$^{-1}$", fontsize = 14)
    plt.ylabel(r"$\epsilon$", fontsize = 14)
    plt.legend(fontsize = 12)
    # plt.close()
    
    plt.figure("Efficiency curve eff = f(TDCR)")
    plt.clf()
    plt.title('  ')
    plt.plot(tdcr, effS,'-b', label = r"$\epsilon_{S}$, $kB$ = "+str(kB)+" cm/keV")
    plt.plot(tdcr, effD,'-k', label = r"$\epsilon_{D}$, $kB$ = "+str(kB)+" cm/keV")
    plt.plot(tdcr, effT,'-r', label = r"$\epsilon_{T}$, $kB$ = "+str(kB)+" cm/keV")
    # plt.xscale("log")
    plt.xlabel(r"$\epsilon_T/\epsilon_D$", fontsize = 14)
    plt.ylabel(r"$\epsilon$", fontsize = 14)
    # bbox_inches = "tight", format = "png", dpi = 500
    plt.legend(fontsize = 12)
    # plt.close()
    
    
def plotEffProfilkB(Rad,SDT): # plot the fitted efficiecny curves for a range of kB 
    kBv = [0.8e-5, 0.9e-5, 1.0e-5, 1.1e-5, 1.2e-5]   # kB vector
    plt.figure("Efficiency curve eff = f(TDCR)")
    plt.clf()
    plt.title('  ')
    for i in kBv:                                    # loop in the kB vector
        pDv = readProfil(Rad,i,"D")[1]               # probability vector of events
        pTv = readProfil(Rad,i,"T")[1]               # probability vector of triple coincidence events
        pSv = readProfil(Rad,i,"S")[1]               # probability vector of single events
        tdcr = np.asarray(pTv)/np.asarray(pDv)       # tdcr vector
        if SDT == "S":
            plt.plot(tdcr, pSv, label = r"$kB$ = "+str(i)+" cm/keV")
            plt.ylabel(r"$\epsilon_S$", fontsize = 14)
        if SDT == "D":
            plt.plot(tdcr, pDv, label = r"$kB$ = "+str(i)+" cm/keV")
            plt.ylabel(r"$\epsilon_D$", fontsize = 14)
        if SDT == "T":
            plt.plot(tdcr, pDv, label = r"$kB$ = "+str(i)+" cm/keV")
            plt.ylabel(r"$\epsilon$", fontsize = 14)
        # plt.xscale("log")
        plt.xlabel(r"$\epsilon_T/\epsilon_D$", fontsize = 14)
        # bbox_inches = "tight", format = "png", dpi = 500
        plt.legend(fontsize = 12)
        plt.savefig("EfficiencyCurves/"+Rad+"/tdcr_"+Rad+".png")
        # plt.close()


## PLOT EFFICIENCY CURVES FOR SEVERAL KB VALUES
# plotEffProfilkB("H-3","D")

## PLOT EFFICIENCY CURVE A GIVEN KB VALUE
# plotEffProfil("Fe-55", 1.05E-5)

## CALCULATION OF THE EFFICIENCY

# AB = 657.296
# BC = 695.919
# AC = 693.579
# T = 448.868
# D = AB+BC+AC-2*T
# out = effTDCR(T/D, "H-3", 1.05E-5) # H-3
# plotEffProfil("H-3", 1.05E-5)

# AB = 864.95
# BC = 1006.58
# AC = 968.45
T = 700
D = 1000
XX = (D + 2*T)/3
print(XX)
print(T/D)
# D = AB+BC+AC-2*T
TDCR = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
for t in TDCR:
    out = effTDCR(t, "H-3", 1.2E-5) # H-3
    eff = out[1]
    u_eff = 0
    A = D/eff
    uA = 0
    
    # DT = 0.977784348 # Co-60 NIST
    # uDT = 0.000711023 # Co-60 NIST
    # D = 84679.21179
    # nMc = 1e3
    # eff_v = []
    # kB1 = 0.8e-5
    # kB2 = 1.2e-5
    # for i in range(int(nMc)):
    #     dti = np.random.normal(DT,uDT)
    #     kBi = np.random.uniform(kB1,kB2)
    #     eff_v.append(effTDCR(dti, "Co-60", kBi)[1]) # Co-60 NIST
    # eff = np.median(eff_v)
    # u_eff = np.std(eff_v)
    # A = D/eff
    # uA = D*u_eff/eff**2
    # result 87099      67             (err = 0.0002) kB = 0.01
    # result 87045      88             kB = [0.008 - 0.012]
    # I2 86640.31523	14.64334414
    # A_NIST = 86920    200
    
    # out = effTDCR(657.296, 695.919, 693.579, 448.868) # Co-60
    digRound = 2
    print("Efficiency of double coincidences = ", round(100*eff, digRound),"+/-", round(100*u_eff, digRound),"%")
    # print("Efficiency of triple coincidences = ", round(100*out[2], digRound),"%")
    # print("Activity from double coincidences = ", round(A, digRound), "+/-", round(uA, digRound) ,"Bq")
    # print("Activity of triple coincidences = ", round(T/out[2], digRound), "Bq")
    # print("Efficiency of double coincidences (asym) = ", round(100*out[3], digRound),"%")
    # print("Efficiency of triple coincidences (asym) = ", round(100*out[5], digRound),"%")
    # print("Activity from double coincidences (asym) = ", round(D/out[3], digRound), "Bq")
    # print("Activity of triple coincidences (asym) = ", round(T/out[5], digRound), "Bq")
    
