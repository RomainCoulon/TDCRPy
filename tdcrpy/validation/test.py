# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 19:15:23 2023
Validation test of the TDCRPy code
@author: romain.coulon
"""
import tdcrpy as td
import numpy as np
import matplotlib.pyplot as plt

# """
# Validation with standard solution for Co-60
# """

# Rad="Co-60"    # list of radionuclides (Na-24)
# pmf_1="1"
# kB =1.0e-5       # Birks constant in cm/keV
# RHO = 0.98
# V = 10
# nE = 1000
# TD = 0.977667386529166        # Measured TDCR value
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# L = 1
# N = 10

# resuts_1 = []
# for i in range(3):
#     resuts_1.append(td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V, "eff", "sym"))
# print(resuts_1)

# # resuts_2=td.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, "sym", N=N)
# # print(resuts_2)


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
# TD = 0.70        # Measured TDCR value
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# L = 1
# N = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]

# print("Rad",Rad)
# for iN in N:
#     resuts_1=td.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, iN, kB, "eff", "sym", Display=False)
#     print(resuts_1[0],resuts_1[1],resuts_1[2],resuts_1[3],resuts_1[4],resuts_1[5])

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
