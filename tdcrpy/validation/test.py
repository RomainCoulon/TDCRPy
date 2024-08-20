# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 19:15:23 2023
Validation test of the TDCRPy code
@author: romain.coulon
"""
import tdcrpy as td
import numpy as np
import matplotlib.pyplot as plt
import cProfile
import pstats
import random
from scipy.optimize import curve_fit
from tqdm import tqdm
import sys
sys.path.insert(1, 'G:\\Python_modules\\BIPM_RI_PyModules')
import TDCRcalculation as cl


L = 1
# L = (1.1, 1.05, 1.15)
TD = 0.977667386529166
# TD = (0.977667386529166, 0.992232838598821, 0.992343419459002, 0.99275350064608)
# TD = (0.977667386529166, 0.995232838598821, 0.990343419459002, 0.99275350064608)
Rad="Co-60"
pmf_1="1"
N = 1000
kB =1.0e-5
V = 10
mode = "eff"

out = td.TDCRPy.TDCRPy(L, Rad, pmf_1, N, kB, V, Display = False, record = False, readRecHist = False)
print(out)

# out = td.TDCRPy.TDCRPy(L, Rad, pmf_1, N, kB, V, Display = False, record = True, readRecHist = False)
# print("result", out)
# out = td.TDCRPy.TDCRPy(L, Rad, pmf_1, N, kB, V, Display = False, record = False, readRecHist = True)
# print("result", out)

# out = td.TDCRPy.eff(TD, Rad, pmf_1, kB, V, N=1000, L=1, maxiter=20, xatol=1e-7)
# print(out)






"""
Quenching
"""

# E=15 # keV
# nE=1000
# kB_vec = np.linspace(0.6e-2, 1.5e-2, 20)

# Em = []
# for kBi in kB_vec:
#     out = td.TDCR_model_lib.Em_e(E*1e3, E*1e3, kBi, nE)*1e-3
#     print(out)
#     Em.append(out)
    
# plt.figure(r"Em($E$) vs kB")
# plt.clf()
# plt.plot(kB_vec, Em, "-k", label=r'$E_i$ = 15 keV')
# # plt.colorbar()
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.legend(fontsize=16)
# plt.xlabel(r"$kB$ / (cm/MeV)", fontsize=18)
# plt.ylabel(r"Em($E$) / keV", fontsize=18)
# # plt.xscale("log")
# # plt.yscale("log")
# plt.show()


"""
Sample in bivar distri
"""


# E = td.TDCR_model_lib.energie_dep_beta2(1000,v=10)
# E = td.TDCR_model_lib.energie_dep_gamma2(10,v=10)
# print(E)


"""
Sampling in pmf
"""

# A = ["A","B","C","D"]
# pmf = [0 , 0.2, 0, 0.8]

# ind = td.TDCR_model_lib.sampling(pmf)
# print(A[ind])


"""
Sampling in pdf
"""
# e = [0, 1, 2, 3, 4, 5]
# pdf = [0.2,0.1,0.3,0.4,0.2]
# print(sum(pdf))
# pdf/=sum(np.asarray(pdf))
# N = 10000
# ev=[]
# for i in tqdm(range(N), desc="Processing", unit=" bins"):
#     ind = td.TDCR_model_lib.sampling(pdf) # sample in pdf
#     # ev.append(energie_dep_beta2(e[ind],V))
#     ev.append(e[ind])
# counts, bins = np.histogram(ev, bins=e, density=True)
# p2=counts/sum(counts)

# em0 = sum(np.asarray(e[:-1])*np.asarray(pdf))
# em1 = sum(bins[:-1]*p2)
# print(f"\nmean emitted E = {em0} keV {len(e)} {len(pdf)}")
# print(f"mean deposited E = {em1} keV {len(bins)} {len(p2)}\n")
# print(p2[-1])

# # bin_centers = (bins[:-1] + bins[1:]) / 2
# plt.figure("0")
# plt.clf()
# # plt.bar(bin_centers, p2, width=(bins[1] - bins[0]), color='g', alpha=0.6, label="deposited")
# plt.plot(bins[:-1], p2, '-g', alpha=0.6, label="deposited")
# plt.plot(e[:-1], pdf,'-r', alpha=0.6, label="betaShape")
# plt.legend()
# plt.xlabel("$E$ /keV")
# plt.ylabel(r"$p$ /keV$^{-1}$")

"""
Spectrum sampling and suming
"""

# # 3H em = 5,68 keV
# e_b,p_b = td.TDCR_model_lib.readBetaShape("H-3","beta-","tot",contH=False)   # read the data of BetaShape
# e_b = np.asarray(e_b)
# p_b = np.asarray(p_b)
# print(sum(p_b))
# p_b/=sum(p_b)

# # # integration
# # em =0
# # for i, pi in enumerate(p_b):
# #     em += pi*e_b[i]
# em = sum(e_b[:-1]*p_b)

# print("integration ", em)

# # sampling
# e_b,p_b = td.TDCR_model_lib.readBetaShape("H-3","beta-","tot",contH=False)   # read the data of BetaShape
# e_b = np.asarray(e_b)
# p_b = np.asarray(p_b)
# p_b/=sum(p_b)
# ei = []
# N = 1000000
# for i in range(N):
#     ind = td.TDCR_model_lib.sampling(p_b)
#     ei.append(e_b[ind])
# em = np.mean(ei)
# u_em = np.std(ei)/np.sqrt(N)

# print("sampling ", em, u_em)

# counts, bins = np.histogram(ei, bins=e_b, density=True)
# p2=counts/sum(counts)

# # bin_centers = (bins[:-1] + bins[1:]) / 2
# plt.figure("0")
# plt.clf()
# # plt.bar(bin_centers, p2, width=(bins[1] - bins[0]), color='g', alpha=0.6, label="deposited")
# plt.plot(bins[:-1], p2, '-g', alpha=0.6, label="deposited")
# plt.plot(e_b[:-1], p_b,'-r', alpha=0.6, label="betaShape")
# plt.legend()
# plt.xlabel("$E$ /keV")
# plt.ylabel(r"$p$ /keV$^{-1}$")


"""
Micelle
"""
# l = 1
# TD =0.9782008792186575
# Rad = "Co-60"
# pmf_1 = "1"
# N=1
# kB = 1e-5
# V=16
# mode="eff"
# mode2="sym"
# TAB = 0.9926063247051856
# TBC = 0.992839996182854
# TAC = 0.9924331860555556

# # out = td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V, mode, mode2, barp= False, Display = True, Smodel=True)
# out = td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V, mode, mode2, barp= False, Display = False, Smodel=False)
# # out = td.TDCRPy.TDCRPy(l, tD[qs],tAB[qs],tBC[qs],tAC[qs], rad, "1", 1, kB*1e-3, 16, "res", "sym", Smodel=Stochastic)
# # out = td.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, V, mode2)
# print(out)


"""
Build procedure
"""
# Lstart = 0.1
# Lend = 5
# Lstep = 0.1
# TD =0.89465
# rad = "C-14"
# pmf_1 = "1"
# N=1000
# kB = 1e-5
# V=10
# mode="res"
# mode2="sym"
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608

# nE=200
# Nstep=4
# span=Lstep
# Smodel0=True
# for j in range(Nstep):
#     if j==0:
#         L = np.arange(Lstart,Lend,span)
#         res =[]
#         for l in tqdm(L, desc="Processing", unit=" trail"):
#             # out = td.TDCR_model_lib.modelAnalytical(l, TD,TAB,TBC,TAC, rad, kB*1e-3, V, "res", "sym", nE)
#             out = td.TDCRPy.TDCRPy(l, TD,TAB,TBC,TAC, rad, "1", N, kB*1e-3, V, "res", "sym", Smodel=Smodel0)
#             res.append(out)
#         L0=L[np.argmin(res)]
#     else:
#         span/=10
#         L1 = np.arange(L0-10*span,L0+10*span,span)
#         res =[]
#         for l in tqdm(L1, desc="Processing", unit=" trail"):
#             # out = td.TDCR_model_lib.modelAnalytical(l, TD,TAB,TBC,TAC, rad, kB*1e-3, V, "res", "sym", nE)
#             out = td.TDCRPy.TDCRPy(l, TD,TAB,TBC,TAC, rad, "1", N, kB*1e-3, V, "res", "sym", Smodel=Smodel0)
#             res.append(out)
#         L0=L1[np.argmin(res)]
# out = td.TDCRPy.TDCRPy(L0, TD,TAB,TBC,TAC, rad, "1", N, kB*1e-3, V, "eff", "sym", Smodel=Smodel0); eff=out[2]
# eff1 = td.TDCRPy.TDCRPy(L0-span, TD,TAB,TBC,TAC, rad, "1", N, kB*1e-3, V, "eff", "sym", Smodel=Smodel0)[2]
# eff2 = td.TDCRPy.TDCRPy(L0+span, TD,TAB,TBC,TAC, rad, "1", N, kB*1e-3, V, "eff", "sym", Smodel=Smodel0)[2]
# ureff = abs(eff2-eff1)/(2*eff)

# print("\nTDCRPy",L0,eff,ureff)

# out=cl.I2calc(TD, TAB, TBC, TAC, rad, kB)
# print("I2calc",out[0]*1e3,out[2])



"""
energie_dep_gamma2(e_inci,v)
"""
# V=np.arange(8,21,0.1)
# N=10000
# E=15 # keV
# e_vec = []
# ue_vec = []
# for v in V:
#     x=[]
#     for i in range(N):
#         out=td.TDCR_model_lib.energie_dep_gamma2(15,v)
#         x.append(out)
#     e_vec.append(np.mean(x))
#     ue_vec.append(np.std(x)/np.sqrt(N))
#     print(v, e_vec[-1], ue_vec[-1])

# plt.figure(r"mean($E_d$) vs volume")
# plt.clf()
# plt.errorbar(V, e_vec, yerr=ue_vec, fmt="-k", label=r'$E_i$ = 15 keV')
# # plt.colorbar()
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.legend(fontsize=16)
# plt.xlabel(r"$v$ /mL", fontsize=18)
# plt.ylabel(r"mean($E_d$)/(keV)", fontsize=18)
# # plt.xscale("log")
# # plt.yscale("log")
# plt.show()


"""
energie_dep_beta2(e_inci,v)
"""
# V=np.arange(0,21,1)
# N=10000
# E=1500 # keV
# for v in V:
#     x=[]
#     for i in range(N):
#         out=td.TDCR_model_lib.energie_dep_beta2(E,v)
#         x.append(out)
#     print(v, np.mean(x), np.std(x)/np.sqrt(N))

"""
simple Eff
"""
# L = 1
# TD =0.89465



# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608

# out = td.TDCRoptimize.eff(TD,TAB,TBC,TAC,"Fe-55","1",1e-5,10,"sym",N=10000)

# out = td.TDCRoptimize.readEfficiency(TD, "Fe-55", "1", 1.0e-5, 10)


# print("TDCRPy sym,", out)

"""
Write efficency curves
"""

# def writeEffcurves(tdcr, uTdcr, effD, ueffD, rad, p, kB, V):
#     file = open("../EfficiencyCurves/"+''.join(rad)+"/Eff_"+''.join(rad)+'_'+''.join(str(p))+'_'+str(kB)+'_'+str(V)+".txt","w")
#     for i, ti in enumerate(tdcr):
#         file.write(str(ti)+" "+str(uTdcr[i])+" "+str(effD[i])+" "+str(ueffD[i])+"\n")
#     file.close()

# # L = np.arange(0.1,2,0.1)
# # L = np.arange(0.001,10,0.1)
# # L = np.logspace(-4, 1, num=10)
# L = np.logspace(-1, 1, num=100)
# # L = np.arange(0.1,0.5,0.1)
# TD = 0.854
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# Rad="Fe-55"
# pmf_1="1"
# N = 10000
# # N = 300
# kB = [0.6e-5, 0.7e-5, 0.8e-5, 0.9e-5, 1.0e-5, 1.1e-5, 1.2e-5, 1.3e-5, 1.4e-5, 1.5e-5]  
# # kB = [1.0e-5]
# V = 10
# mode = "eff"
# mode2 = "sym"

# for kBi in kB:
#     tdcr_v = []
#     utcrv_v = []
#     eff_v = []
#     ueff_v = []
#     for l in L:
#         print(f"free parameter = {round(l,3)} keV-1")
#         out = td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, Rad, pmf_1, N, kBi, V, mode, mode2, barp=True)
#         tdcr_v.append(out[4]/out[2])
#         utcrv_v.append(tdcr_v[-1]*np.sqrt((out[3]/out[2])**2+(out[5]/out[4])**2))
#         eff_v.append(out[2])
#         ueff_v.append(out[3])
    
#     writeEffcurves(tdcr_v,utcrv_v,eff_v,ueff_v,Rad, pmf_1, kBi, V)
    
# plt.figure("Efficiency curve")
# plt.clf()
# plt.errorbar(tdcr_v, eff_v, xerr=utcrv_v, yerr=ueff_v)
# plt.xlabel('TDCR')
# plt.ylabel('Efficiency')
# plt.show()


"""
Read efficency curves
"""


# rad = "H-3"
# L = 1


# td.TDCRPy.TDCRPy(1, 0.7, 0.7, 0.7, 0.7, "Fe-55", "1", 10, 1e-5, 10, "sym", "eff", Display=True)


"""
Plot stopping power
"""


# plt.clf()

# # E = np.logspace(0, 6, 100) # for electron
# E = np.logspace(-3, 4, 100) # for alpha
# # E = np.linspace(2900, 3000, 100) # for alpha

# w = []
# for Ei in E:
#     # Em.append(td.TDCR_model_lib.Em_e(Ei, Ei, kBi*1e3, 10000)) #  for electron
#     w.append(td.TDCR_model_lib.stoppingpowerA(Ei)) # for alpha
#     # w.append(stoppingpowerA(Ei)) # for alpha

# def second_order_poly(x, a, b, c):
#     return a * x**2 + b * x + c

# popt, pcov = curve_fit(second_order_poly, E, w)
# print(popt)



# plt.plot(E, w)
# # plt.plot(E, popt[0]*E**2+popt[1]*E+popt[2],label='fit')
# # plt.colorbar()
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.legend(fontsize=16)
# plt.xlabel(r"$E$ /keV", fontsize=18)
# plt.ylabel(r"w/(keV.cm2/g)", fontsize=18)
# plt.xscale("log")
# # plt.yscale("log")
# plt.show()






# plt.clf()

# kB = np.linspace(0.6e-5, 1.5e-5, 3)
# print(kB)
# for kBi in kB:
#     E = np.logspace(0, 6.5, 100) # for electron
#     # E = np.logspace(0, 4, 100) # for alpha
#     Em = []
#     for Ei in E:
#         Em.append(td.TDCR_model_lib.Em_e(Ei, Ei, kBi*1e3, 1000, Et=1500)) #  for electron
#         # Em.append(td.TDCR_model_lib.Em_a(Ei, kBi, 10000)) # for alpha
    
    
#     # plt.plot(E, Em/E, label = str(round(kBi*1e3,3))+" cm/MeV")
#     plt.plot(E*1e-3, Em/E, label = str(round(kBi*1e3,3))+" cm/MeV")
#     # plt.colorbar()
#     plt.xticks(fontsize=16)
#     plt.yticks(fontsize=16)
#     plt.legend(fontsize=16)
#     plt.xlabel(r"$E$ /keV", fontsize=18)
#     plt.ylabel(r"Em($E$)/E", fontsize=18)
#     plt.xscale("log")
#     # plt.yscale("log")
#     plt.show()

# A = td.TDCR_model_lib.stoppingpowerA(0.008)
# B = td.TDCR_model_lib.E_quench_a(1,1e-5,1000)
# print(A, B)

"""
Read response matrixes
"""

# from scipy.ndimage import gaussian_filter
# import cv2


# epsilon = 1e-5
# limite_sup = 1e-2
# sht = 0

# A = td.TDCR_model_lib.Matrice10_e_1
# A = td.TDCR_model_lib.Matrice10_e_2
# A = td.TDCR_model_lib.Matrice10_e_3
# # A = td.TDCR_model_lib.Matrice10_p_1
# # A = td.TDCR_model_lib.Matrice10_p_2
# # A = td.TDCR_model_lib.Matrice10_p_3
# C = np.flipud(A[1:])
# C = np.clip(C, a_min=epsilon, a_max=limite_sup)
# C = np.log(C)
# C = cv2.GaussianBlur(C, (5, 5), 20)
# print("step",A[0][1]-A[0][0], min(A[0]), max(A[0]))
# extent = [A[0,0], A[0,-1], A[0,0], A[0,-1]]
# x = np.arange(A[0,0], A[0,-1], A[0,-1]/9)
# y = np.arange(A[0,0], A[0,-1], A[0,-1]/9)


# plt.clf()
# plt.imshow(C, extent=extent, cmap='viridis', interpolation='nearest')
# plt.colorbar()
# plt.xticks(np.arange(A[0,0]-sht, A[0,-1]-sht, A[0,-1]/10))
# plt.yticks(np.arange(A[0,0]-sht, A[0,-1]-sht, A[0,-1]/10))
# plt.xlabel(r"$E_i$ /keV", fontsize=14)
# plt.ylabel(r"$E_d$ /keV", fontsize=14)
# plt.show()

"""
Tests decay data uncertainty propagation
"""
# L = 1
# TD = 0.977667386529166
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# Rad="Fe-55"
# pmf_1="1"
# N = 10
# kB =1.0e-5
# V = 10
# # mode = "dis"
# mode ="eff"
# mode2 = "sym"

# out = td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V, mode, mode2, Display=True, barp=False, uncData=False)
# if mode == "dis": td.TDCR_model_lib.display_distrib(out[0], out[1], out[2])
# else: print(out)

"""
Efficiency curves !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# L = np.logspace(-3,2,100) # H-3
# L = np.logspace(-3,2,100) # C-14
# L = np.logspace(-3,2,10) # Sr-89
# L = np.logspace(-3,2,10) # Fe-55
# L = np.logspace(-1,1,100) # Fe-55
# TD = 0.977667386529166
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# # rad_v = ["Fe-55"]
# rad_v = ["F-18"]
# pmf_1="1"
# V = 10
# mode = "eff"
# mode2 = "sym"
# N = 10000
# #kB = [0.6e-5, 1.0e-5, 1.5e-5]
# # kBtxt = [6, 10, 15]
# kB = [1.0e-5]
# kBtxt = [10]


# plt.figure("visu")
# for Rad in rad_v:
#     for i, kBi in enumerate(kB):
#         plt.clf()
#         f = open("result_A_"+Rad+"_"+str(kBtxt[i])+".txt", "w") # analystical
#         effD = []; effT = []; ueffD = []
#         for x, l in enumerate(L):
#             print('progress',100*x/len(L),' %')
#             out1 = td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, Rad, pmf_1, N, kBi, V, mode, mode2)
#             out2 = td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, Rad, pmf_1, N, kBi, V, mode, mode2)
#             out3 = td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, Rad, pmf_1, N, kBi, V, mode, mode2)
#             out4 = td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, Rad, pmf_1, N, kBi, V, mode, mode2)
#             out5 = td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, Rad, pmf_1, N, kBi, V, mode, mode2)
#             meanS = np.mean([out1[0], out2[0], out3[0], out4[0], out5[0]])
#             stdS = np.std([out1[0], out2[0], out3[0], out4[0], out5[0]])
#             meanD = np.mean([out1[2], out2[2], out3[2], out4[2], out5[2]])
#             stdD = np.std([out1[2], out2[2], out3[2], out4[2], out5[2]])
#             meanT = np.mean([out1[4], out2[4], out3[4], out4[4], out5[4]])
#             stdT = np.std([out1[4], out2[4], out3[4], out4[4], out5[4]])            
#             f.write(str(l)+" "+str(meanS)+" "+str(stdS)+" "+str(meanD)+" "+str(stdD)+" "+str(meanT)+" "+str(stdT)+" \n")
#             effD.append(meanD)
#             ueffD.append(stdD)
#             effT.append(meanT)
#         f.close()
#         effD=np.asarray(effD)
#         ueffD=np.asarray(ueffD)
#         effT=np.asarray(effT)
#         plt.errorbar(effT/effD, effD, yerr=ueffD, label=str(kBi))
        # plt.xscale("lin")

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# td.TDCRPy.TDCRPy(1, 0.97, 0.97, 0.97, 0.97, "Fe-55", "1", 10000, 1.0e-5, 10, "eff", "sym", Display=False,barp=True)

# # l = 1
# kBi = 1.0e-5
# f = open("result_AN_.txt", "w") # analystical
# for Rad in rad_v:
#     R = []
#     print(Rad)
#     for j in range(5):
#         out = td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, Rad, pmf_1, N, kBi, V, mode, mode2, barp=False)
#         R.append(out[2])
#     f.write(Rad+" "+str(np.median(R))+" "+str(np.std(R))+" \n")
# f.close()


# A=td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, "Tc-99", "1", 10000, kBi, V, mode, mode2, barp=True,Display= False) # MC 0.9771859524255478
# print(A) # TDCR17 0.9724
# e, p = td.TDCR_model_lib.readBetaShape("Tc-99", 'beta-', 'tot')    # no Ed 0.9690042889831171 - Ed 0.9690042792604032
# e2, p2 = td.TDCR_model_lib.readBetaSpectra("Tc-99")                # no Ed 0.9712517070635934

# A=td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, "S-35", "1", 10000, kBi, V, mode, mode2, barp=False, Display= False) # MC 0.9357931563925717 0.9404170688549414 0.9404404083209962
# print(A) # TDCR17 0.9443
# e, p = td.TDCR_model_lib.readBetaShape("S-35", 'beta-', 'tot')     # no Ed  0.9394876324351543 - Ed 0.9394854546918996 E_quench 0.9360969563486823
# e2, p2 = td.TDCR_model_lib.readBetaSpectra("S-35")                # no Ed 0.9432625905405294

#A=td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, "Ni-63", "1", 10, kBi, V, mode, mode2, barp=False, Display= True) # MC 0.8314873609197219 0.8291652424879523 0.8328369624432058
# print(A) # TDCR17 
# e, p = td.TDCR_model_lib.readBetaShape("Ni-63", 'beta-', 'tot')     # no Ed  0.820637703 
# e2, p2 = td.TDCR_model_lib.readBetaSpectra("Ni-63")                # no Ed 

# A=td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, "Pu-241", "1", 10000, kBi, V, mode, mode2, barp=False, Display= False) # MC 0.5281756667545862 0.520045761411465 0.5308791207244081
# print(A) # TDCR17 
# e, p = td.TDCR_model_lib.readBetaShape("Pu-241", 'beta-', 'tot')     # no Ed  0.514531626, 0.5200989957785515
# e2, p2 = td.TDCR_model_lib.readBetaSpectra("Pu-241")                # no Ed 


# rad = "S-35"
# e, p = td.TDCR_model_lib.readBetaShape(rad , 'beta-', 'tot')    
# e2, p2 = td.TDCR_model_lib.readBetaSpectra(rad) 
# e=np.asarray(e); e2=np.asarray(e2); p=np.asarray(p); p2=np.asarray(p2)



# Nb = 100000
# Ev = []
# Ed = []
# Em = []
# for i in range(Nb):
    
    
#     index_beta_energy = td.TDCR_model_lib.sampling(p)                           # sampling energy of beta
#     Ev.append(e[index_beta_energy])
    
#     # Ev.append(random.choices(e, weights=p, k=1)[0])

    
    
#     Ed.append(td.TDCR_model_lib.energie_dep_beta2(Ev[-1],v=V))
#     Em.append(td.TDCR_model_lib.Em_e(Ev[-1]*1e3, Ed[-1]*1e3, kBi*1e3, 1000)*1e-3)

# em2 = []
# for i, ei in enumerate(e2):
#     em2.append(td.TDCR_model_lib.Em_e(ei*1e3,ei*1e3,kBi*1e3,1000)*1e-3)


# print("Mean of the distribution (emission) ",sum(e*p),np.mean(Ev),np.std(Ev)/np.sqrt(Nb))
# print("Mean of the distribution (interaction) ",sum(e2*p2),np.mean(Ed),np.std(Ed)/np.sqrt(Nb))
# print("Mean of the distribution (quenched) ",sum(em2*p2),np.mean(Em),np.std(Em)/np.sqrt(Nb))
# # print(sum(e*p),sum(e2*p2))


# plt.figure("beta")
# plt.clf()
# # plt.plot(e,p,"+k",label='BetaShape')
# # plt.plot(e2,p2,"+r",label='Beta spectrum')
# plt.hist(Ev,bins=e,label='Sample (emission)',density=True)
# plt.hist(Ed,bins=e,label='Sample (deposition)',density=True, alpha = 0.5)
# plt.hist(Em,bins=e,label='Sample (quenched)',density=True, alpha=0.5)
# plt.legend()







# A=[]; B=[]
# for j in range(10000):
#     A.append(td.TDCR_model_lib.energie_dep_beta2(5000, 10))
#     B.append(td.TDCR_model_lib.energie_dep_beta(5000))

# print(np.mean(A), np.mean(B))

"""
Efficiency curves MC
"""

# L = np.arange(0.5,2,0.1)
# TD = 0.977667386529166
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# Rad="Sr-89"
# pmf_1="1"
# N = 10000
# kB =1.0e-5
# V = 10
# mode = "eff"
# mode2 = "sym"

# out = td.TDCRoptimize.effCurves(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V)
# print("L")
# for i in L:
#     print(i)
# print("tdcr")
# for i in out[2]:
#     print(i)
# print("u(tdcr)")
# for i in out[3]:
#     print(i)
# print("Eff")
# for i in out[0]:
#     print(i)
# print("u(Eff)")
# for i in out[1]:
#     print(i)


"""
Test display
"""
# td.TDCRPy.TDCRPy(1, 0, 0, 0, 0, "Fe-55", "1", 1, 1e-5, 10, "eff", "sym", Display=True, barp=False)


"""
Compare with I2calc()
"""
# import tdcrpy
# from TDCRcalculation import I2calc
# # data HMI std 42

# AB =  1663.601
# BC =  1668.343
# AC =  1664.642
# D =  1974.418
# T =  1511.084

# ## version I2calc
# F1, FFF, effa, effb=I2calc(T/D,T/AB,T/BC,T/AC,"H-3",1e-2)
# FOM = F1*1e3
# FOMA = FFF[0]*1e3
# FOMB = FFF[1]*1e3
# FOMC = FFF[2]*1e3
# EFF1 = effa
# EFF2 = effb
# I2 = D/EFF2
# print("I2Calc,", FOM, EFF1, EFF2, I2)


# F1, _, _, _, effa, _, _, _ = tdcrpy.TDCRoptimize.eff(T/D,T/AB,T/BC,T/AC,"H-3","1",1e-5,10,"sym",N=10)
# FOM = F1
# EFF1 = effa
# I2 = D/EFF1
# print("TDCRPy sym,", FOM, EFF1, I2)

# _, FFF, _, _, effb, _, _, _ = tdcrpy.TDCRoptimize.eff(T/D,T/AB,T/BC,T/AC,"H-3","1",1e-5,10,"asym",N=10)
# FOMA = FFF[0]
# FOMB = FFF[1]
# FOMC = FFF[2]
# EFF2 = effb
# I2 = D/EFF2
# print("TDCRPy asym,", FOMA, FOMB, FOMC, EFF2, I2)

"""
Profiling
"""
# TD = 0.977667386529166        # Measured TDCR value
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# L = 1.0
# e, p = td.TDCR_model_lib.readBetaShape("H-3", 'beta-', 'tot')

# def main():
    
#     """
#     test TDCRPy
#     """
#     # td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, "Sr-89", "1", 1, 1.0e-5, 10, "eff", "sym")
#     # td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, "Co-60", "1", 1, 1.0e-5, 10, "eff", "sym")
#     td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, "Am-241", "1", 1, 1.0e-5, 10, "eff", "sym")
#     # td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, "H-3", "1", 1, 1.0e-5, 10, "eff", "sym")
    
    
    
    
#     """
#     test sampling
#     """
#     # for j in range(1000): td.TDCR_model_lib.sampling(p)
#     # for j in range(1000): random.choices(range(len(p)), weights=p)[0]
    
# if __name__ == "__main__":
#     main()

# cProfile.run("main()", filename="profile_results.txt")
# stats = pstats.Stats("profile_results.txt")
# stats.sort_stats(pstats.SortKey.TIME)
# stats.print_stats()


"""
for the quick start
"""

# import tdcrpy

# L = (1.5, 1.2, 1.4)
# TD = 0.977667386529166
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# Rad="Co-60"
# pmf_1="1"
# N = 100
# kB =1.0e-5
# V = 10
# mode = "eff"
# mode2 = "asym"


# result = tdcrpy.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V, mode, mode2)

# print("R = ",result)

# print("efficiency S = ", round(result[0],4), "+/-", round(result[1],4))
# print("efficiency D = ", round(result[2],4), "+/-", round(result[3],4))
# print("efficiency T = ", round(result[4],4), "+/-", round(result[5],4))


# import tdcrpy

# TD = 0.977667386529166
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# Rad="Co-60"
# pmf_1="1"
# N = 250
# kB =1.0e-5
# V = 10
# mode2 = "asym"

# result = tdcrpy.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, V, mode2, N=N)

# print("Global free parameter = \t", round(result[0],4), " keV-1")
# print("Free parameter (PMT A) = \t", round(result[1][0],4) , " keV-1")
# print("Free parameter (PMT B) = \t", round(result[1][1],4) , " keV-1")
# print("Free parameter (PMT C) = \t", round(result[1][2],4) , " keV-1")
# print("efficiency S = \t", round(result[2],4), "+/-", round(result[3],4))
# print("efficiency D = \t", round(result[4],4), "+/-", round(result[5],4))
# print("efficiency T = \t", round(result[6],4), "+/-", round(result[7],4))




"""
Validation with standard solution for Co-60 (comparison 2023)
"""

# file_path = "result_Co60_sym.txt"  # Replace "example.txt" with the path of your desired file.
# file = open(file_path, "w")

# rad="Co-60"    # list of radionuclides (Na-24)
# pmf_1="1"
# kB =1.0e-5       # Birks constant in cm/keV
# RHO = 0.98
# V = 10
# N = 1000
# TD = 0.977667386529166        # Measured TDCR value
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# # TAB = 0.992232838598821
# # TBC = 0.992343419459002
# # TAC = 0.99275350064608
# L = 1.5
# # N = [10, 20, 50, 100, 200, 500, 800, 1000, 2000, 3000, 5000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 60000, 70000, 80000, 90000, 100000]

# effKCDB =   0.97426
# u_effKCDB = 0.00032

# resuts_2=td.TDCRoptimize.eff(TD, TAB, TBC, TAC, rad, pmf_1, kB, V, "asym", N=10000)
# print('asym 10000 kB = 0.010')
# print(resuts_2)


# span = 5
# Nstep = 3
# Lstart=0.5 
# Lend=3 
# for jst in range(Nstep):
#     if jst==0:
#         L = np.arange(Lstart,Lend,span)
#         res0 =[]
#         for l in tqdm(L, desc="Processing", unit=" trail"):
#             # out = td.TDCR_model_lib.modelAnalytical(l, TD[j],TAB[j],TBC[j],TAC[j], rad, kB*1e-3, 10, "res", "sym", nE)
#             out = td.TDCRPy.TDCRPy(l, TD, TAB, TBC, TAC, rad, "1", 1000, kB, 10, "res", "sym")
#             res0.append(out)
#         L0=L[np.argmin(res0)]
#     else:
#         span/=10
#         L1 = np.arange(L0-10*span,L0+10*span,span)
#         res =[]
#         for l in tqdm(L1, desc="Processing", unit=" trail"):
#             # out = td.TDCR_model_lib.modelAnalytical(l, TD[j],TAB[j],TBC[j],TAC[j], rad, kB*1e-3, 10, "res", "sym", nE)
#             out = td.TDCRPy.TDCRPy(l, TD,TAB,TBC,TAC, rad, "1", 1000, kB*1e-3, 10, "res", "sym")
#             res.append(out)
#         L0=L1[np.argmin(res)]
# out = td.TDCRPy.TDCRPy(L0, TD,TAB,TBC,TAC, rad, "1", 10000, kB*1e-3, 10, "eff", "sym"); eff=out[2]
# eff1 = td.TDCRPy.TDCRPy(L0-span, TD,TAB,TBC,TAC, rad, "1", 10000, kB*1e-3, 10, "eff", "sym")[2]
# eff2 = td.TDCRPy.TDCRPy(L0+span, TD,TAB,TBC,TAC, rad, "1", 10000, kB*1e-3, 10, "eff", "sym")[2]
# ureff = (eff2-eff1)/(2*eff)
# #### procedure TDCRPy

# plt.plot(L,res0)   # filtered 2
# plt.plot(L0,min(res0),"ok")

# print("efficiency",out[2])
# print("true eff", 0.9743)

# Symmetrical model
# for Ni in N:
#     print("symetrical model")
#     # td.TDCR_model_lib.tic()
#     resuts_2=td.TDCRoptimize.eff(TD, TAB, TBC, TAC, rad, pmf_1, kB, V, "sym", N=Ni)
#     # a = td.TDCR_model_lib.toc()
#     print("/n",Ni,resuts_2,"/n/n")
#     file.write(str(resuts_2[0])+" "+str(resuts_2[2])+" "+str(resuts_2[3])+"\n")
# file.close()

# file_path = "result_Co60_asym.txt"  # Replace "example.txt" with the path of your desired file.
# file = open(file_path, "w")

# # Asymmetrical model
# for Ni in N:
#     print("asymetrical model")
#     # td.TDCR_model_lib.tic()
#     resuts_2=td.TDCRoptimize.eff(TD, TAB, TBC, TAC, rad, pmf_1, kB, V, "asym", N=Ni)
#     # a = td.TDCR_model_lib.toc()
#     print("\n",Ni,resuts_2,"\n\n")
#     file.write(str(resuts_2[0])+" "+str(resuts_2[2])+" "+str(resuts_2[3])+"\n")
# file.close()

"""
Validation with Fe-55
"""

# Rad="Co-60"    # list of radionuclides (Na-24)
# pmf_1="1"
# kB =1.0e-5       # Birks constant in cm/keV
# RHO = 0.98
# nE = 1000
# D = 1000
# T = 700        # Measured TDCR value
# TAB = 800
# TBC = 800
# TAC = 800
# L = 1
# N = 5
# TD = T/D

# # resuts_1=td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, "eff", "sym")
# # print(resuts_1)

# resuts_2=td.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, "sym", N=N)
# print("\n",resuts_2)



"""
Validation with the analytical model
"""

# Rad="Sr-89"
# pmf_1="1"
# kB =1.0e-5       # Birks constant in cm/keV
# V = 10
# TD = 0.70        # Measured TDCR value
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# L = 1
# N = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]

# print("Rad",Rad)
# for iN in N:
#     td.TDCR_model_lib.tic()    
#     resuts_1=td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, iN, kB, V, "eff", "sym", Display=False)
#     print(iN,resuts_1[0],resuts_1[1],resuts_1[2],resuts_1[3],resuts_1[4],resuts_1[5])
#     td.TDCR_model_lib.toc()
# resuts_2=td.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, "sym", N=N)
# print(resuts_2)

# """
# Validation with the analytical model
# """

# kB =1.0e-5       # Birks constant in cm/keV
# nE=[10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000]
# e = 1000 # energy incidente
# for inE in nE:
#     # results = td.TDCR_model_lib.E_quench_e(e*1e3, kB*1e3, inE)*1e-3
#     results = td.TDCR_model_lib.E_quench_a(e,kB,inE)
#     print(results)

#===================================================================================================
#============================tracer stopping power of electron and alpha ===========================

# e1 = np.linspace(10,1.99e4,10000)
# e2 = np.linspace(1.99e4,1e8,20000)
# sp1 = []
# for i in e1:
#     sp1.append(td.TDCR_model_lib.stoppingpower(i))
# sp2=[]
# for i in e2:
#     sp2.append(td.TDCR_model_lib.stoppingpower(i))

# ea = np.linspace(1,8e3,4000)
# spa = []
# for i in ea:
#     spa.append(td.TDCR_model_lib.stoppingpowerA(i)*1e-3)
# x = ea*1e3
# plt.plot(e1,sp1,color='blue',label="pouvoir d'arrête électron")
# plt.plot(e2,sp2,color='blue')
# plt.plot(x,spa,label="pouvoir d'arrête alpha")
# plt.xscale('log')
# plt.yscale('log')
# plt.title("pouvoir d'arrête d'electron et alpha")
# plt.ylabel("pouvoir d'arrête/MeV.cm-1")
# plt.xlabel('énergie cinétique/eV')
# plt.legend(fontsize=10,loc='best')
# plt.savefig('Quenching/stoppingpowerE_A.png')

#=====================================================================================

## Display beta spectra

# def readTDCR17spectra(rad):
#     file = open("decayData/spectra/spectrumTDCR17_"+rad+".dat")
#     data = file.read()
#     file.close()
#     data = data.split("\n")
#     e = []; r = []
#     for i in data:
#         if "keV" not in i:
#             a = i.split("  ")
#             if len(a)>1:
#                 e.append(float(a[1])*1e-3)
#                 r.append(float(a[2]))
#     return e, r

# rplot = "S-35"
# out = td.TDCR_model_lib.readBetaShape(rplot, "beta-", "tot"); print(sum(out[1]))
# out_t = readTDCR17spectra(rplot); print(sum(out_t[1])) # TDCR17 spectra
# plt.figure("beta spectrum")
# plt.clf()
# plt.plot(out[0],out[1]/(out[0][1]-out[0][0]),label='Beta spectrum - BetaShape',ls='-',color='red',lw=3)
# plt.plot(out_t[0],np.asarray(out_t[1])/(out_t[0][1]-out_t[0][0]),label='Beta spectrum - TDCR17',ls='-',color='blue',lw=3)
# plt.xscale('linear')
# plt.legend(fontsize=12,loc='best')
# plt.xlabel(r'$E$ /keV', fontsize=12)
# plt.ylabel(r'd$N$/d$E$ /keV$^{-1}$',fontsize=12)
# plt.savefig("decayData/spectra/BetaSpectrum_"+rplot+".png")

#======================================================================================

#========================= Tracer les courbes avec kB différents ======================


# s1 = []
# s2 = []
# s3 = []
# x = np.linspace(5,8000,4000) 

# for i in x:
#     s1.append(td.TDCR_model_lib.E_quench_a(i,kB=7e-6)/i)
#     s2.append(td.TDCR_model_lib.E_quench_a(i,kB=1e-5)/i)
#     s3.append(td.TDCR_model_lib.E_quench_a(i,kB=1.4e-5)/i)

# plt.plot(x,s1,label='kB=0.007cm/MeV',color='green',lw=2)
# plt.plot(x,s2,label='kB=0.01cm/MeV',ls=':',color='red',lw=3)
# plt.plot(x,s3,label='kB=0.014cm/MeV')
# plt.xscale('log')
# #plt.yscale('log')
# plt.legend(fontsize=12,loc='best')
# plt.xlabel('E de particule/keV')
# plt.ylabel("énergie d'extinction/E")
# plt.savefig("Quenching/beta 100-10K E_Q sur E.png")


# td.TDCR_model_lib.readPenNuc2("H-3")
