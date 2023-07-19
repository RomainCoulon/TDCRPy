# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 10:04:23 2023

@author: romain.coulon
"""
import TDCR_model_lib as tl
import TDCRPy as td
import matplotlib.pyplot as plt

## INPUT OF THE MODEL
# N=1                   # number of simulated decay (MC trials)
N= 10000
Rad="Co-60, H-3"            # list of radionuclides (Na-24)
# Rad = ["Cs-137"]
pmf_1="0.5, 0.5"                # relative abondance (pmf)
kB =[1.0e-5]
# kB = [0.8e-5, 0.9e-5, 1.0e-5, 1.1e-5, 1.2e-5]    # Birks constant in cm/keV
L=[1.13]
# L = np.logspace(-3,1,50) # Free paramete in keV-1

Record = False                  # to record the efficiency curves
# Record = True                  # to record the efficiency curves
# Display = False               # to display calculation results on the console
Display = True                # to display calculation results on the console
# RHO = 0.96         #density of absorber (Toluene) g/cm3
RHO = 0.98           #density of absorber (UG + H20) g/cm3
nE = 1000            #number of bin to discretize the energy vector for scintillation quenching calculation
mode = "eff"
#meta_vec = ["Am-242m","Pa-234m","Pm-148m","Pr-144m","Xe-133m","Te-127m","Ag-110m","Ag-108m","Tc-99m","Nb-95m","Y-90m","Mn-52m"]

TDCR_measure = 0.977784        # Measured TDCR value
u_TDCR_measure = 0.000711      # standard uncertainty
TAB=0.8
TBC=0.8
TAC=0.8

for kB_i in kB: # Loop on the kB
    mean_efficiency_S = []  # efficiency of single counte rate
    std_efficiency_S = []   # std
    mean_efficiency_T = []  # efficiency of triple coincidence counte rate
    std_efficiency_T = []   # std
    mean_efficiency_D = []  # efficiency of double coincidence counte rate
    std_efficiency_D = []   # std
    TDCR_calcul = []
    for L_i in L: # loop on the free parameter values
        efficiency_S = []        # results calculated efficiency for single events
        efficiency_D = []        # results calculated efficiency for double coincidence
        efficiency_T = []        # results calculated efficiency for triple coincidence
        
        out=td.TDCRPy(L_i,TDCR_measure,TAB, TBC, TAC, Rad,pmf_1,10,kB_i,RHO,nE,mode="eff",mode2="sym",Display=True)
        efficiency_S.append(out[0])
        std_efficiency_S.append(out[1])
        efficiency_D.append(out[2])
        std_efficiency_D.append(out[3])
        efficiency_T.append(out[4])
        std_efficiency_T.append(out[5])
        
        if len(L) > 1: print("\n\t\t Progress = ", round(100*(L.tolist().index(L_i)+1)/len(L), 1), " %")
        print("\t\tradionuclide(s): ", Rad)
        print("\t\t TDCR calculation _ kB = ", kB_i, "cm/keV")
        print("\t\t Free parameter = ", L_i, " keV-1")
        print("\t\t Efficiency of Triple coincident events = ", round(100*efficiency_T[-1],3), "+/-", round(100*std_efficiency_T[-1],3), " %")
        print("\t\t Efficiency of Double coincident events = ", round(100*efficiency_D[-1],3), "+/-", round(100*std_efficiency_D[-1],3), " %")
        print("\t\t Efficiency of Single events = ", round(100*efficiency_S[-1],3), "+/-", round(100*std_efficiency_S[-1],3), " %")

    if len(mean_efficiency_S)>1:
        plt.figure("Efficiency curve I")
        plt.clf()
        plt.title(''.join(Rad))
        plt.errorbar(L, mean_efficiency_S, yerr = std_efficiency_S, fmt=".b", label = "S")
        plt.errorbar(L, mean_efficiency_D, yerr = std_efficiency_D, fmt=".k", label = "D")
        plt.errorbar(L, mean_efficiency_T, yerr = std_efficiency_T, fmt=".r", label = "T")
        plt.xscale("log")
        plt.xlabel(r"$L$ /keV$^{-1}$", fontsize = 14)
        plt.ylabel(r"$\epsilon$", fontsize = 14)
        plt.legend(fontsize = 12)
        if Record: plt.savefig("EfficiencyCurves/"+''.join(Rad)+"/fom_"+''.join(Rad)+"_"+str(kB_i)+".png")
        plt.close()
                
        plt.figure("Efficiency curve II")
        plt.clf()
        plt.title(''.join(Rad))
        plt.errorbar(TDCR_calcul, mean_efficiency_S, yerr=std_efficiency_S, fmt=".b", label = "S")
        plt.errorbar(TDCR_calcul, mean_efficiency_D, yerr=std_efficiency_D, fmt=".k", label = "D")
        plt.errorbar(TDCR_calcul, mean_efficiency_T, yerr=std_efficiency_T, fmt=".r", label = "T")
        plt.xlabel(r"$\epsilon_T/\epsilon_D$", fontsize = 14)
        plt.ylabel(r"$\epsilon_D$", fontsize = 14)
        plt.legend(fontsize = 12)
        if Record: plt.savefig("EfficiencyCurves/"+''.join(Rad)+"/tdcr_"+''.join(Rad)+"_"+str(kB_i)+".png")
        plt.close()
    
    if Record:
        tl.writeEffcurves(L, mean_efficiency_S, std_efficiency_S, Rad, pmf_1, kB_i, "S")
        tl.writeEffcurves(L, mean_efficiency_D, std_efficiency_D, Rad, pmf_1, kB_i, "D")
        tl.writeEffcurves(L, mean_efficiency_T, std_efficiency_T, Rad, pmf_1, kB_i, "T")