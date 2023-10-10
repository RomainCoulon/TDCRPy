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

import sys
sys.path.insert(1, 'G:\Python_modules\BIPM_RI_PyModules')


"""
Efficiency curves
"""
L = np.arange(0.5,2,0.1)
TD = 0.977667386529166
TAB = 0.992232838598821
TBC = 0.992343419459002
TAC = 0.99275350064608
Rad="Sr-89"
pmf_1="1"
N = 10000
kB =1.0e-5
V = 10
mode = "eff"
mode2 = "sym"

out = td.TDCRoptimize.effCurves(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V)
print("L")
for i in L:
    print(i)
print("tdcr")
for i in out[2]:
    print(i)
print("u(tdcr)")
for i in out[3]:
    print(i)
print("Eff")
for i in out[0]:
    print(i)
print("u(Eff)")
for i in out[1]:
    print(i)


"""
Test display
"""
# td.TDCRPy.TDCRPy(1, 0, 0, 0, 0, "Cd-109", "1", 1, 1e-5, 10, "eff", "sym", Display=True, barp=False)


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

# file_path = "result.txt"  # Replace "example.txt" with the path of your desired file.
# file = open(file_path, "w")

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
# L = 1.5
# N = [10, 20, 50, 100, 200, 500, 800, 1000, 2000, 3000, 5000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 60000, 70000, 80000, 90000, 100000]

# # Symetrical model
# for Ni in N:
#     print("symetrical model")
#     # td.TDCR_model_lib.tic()
#     resuts_2=td.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, V, "sym", N=Ni)
#     # a = td.TDCR_model_lib.toc()
#     print("/n",Ni,resuts_2,"/n/n")
#     file.write(str(resuts_2[0])+" "+str(resuts_2[2])+" "+str(resuts_2[3])+"\n")

# Asymetrical model
# for Ni in N:
#     print("asymetrical model")
#     # td.TDCR_model_lib.tic()
#     resuts_2=td.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, V, "asym", N=Ni)
#     # a = td.TDCR_model_lib.toc()
#     print("\n",Ni,resuts_2,"\n\n")
#     file.write(str(resuts_2[0])+" "+str(resuts_2[2])+" "+str(resuts_2[3])+"\n")


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
