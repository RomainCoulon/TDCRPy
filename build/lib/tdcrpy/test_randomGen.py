# -*- coding: utf-8 -*-
"""
Created on Wed May 14 16:17:16 2025

@author: romain.coulon
"""
import numpy as np

# m = 0.1
# N = 1000000

# ps0 = 1-np.exp(-m)

# psi=np.random.poisson(m,N)
# ps1 = sum(psi>0)/N
# ups1 = np.sqrt(sum(psi>0))/N

# print(ps0,ps1,ups1)
# print(ps0-ps1,ups1)
# print(abs(ps0-ps1)<2*ups1)

import importlib
import tdcrpy
tdcrpy.TDCR_model_lib.modifyEffQ("0.1, 0.1, 0.1")
# tdcrpy.TDCR_model_lib.modifyOptModel("stochastic-dependence")
tdcrpy.TDCR_model_lib.modifyOptModel("poisson")
L = [1.0, 1.0, 1.0]
e_q = [100]
diffP = 1
importlib.reload(tdcrpy.TDCR_model_lib)

Q = tdcrpy.TDCR_model_lib.readEffQ0()
Q = Q.split(",")
Q = [float(i) for i in Q]
QL = [float(Qi)*L[i] for i, Qi in enumerate(Q)]

e_q2 = [0]; t1 = 0; evenement = 1; extDT = 50; measTime = 60000

S,D,T,_,_,_,_ = tdcrpy.TDCR_model_lib.detectProbabilities(QL, e_q, e_q2, t1, evenement, extDT, measTime)
SmcI=[];DmcI=[];TmcI=[]
nIter=100000
for i in range(nIter):
    Smc,Dmc,Tmc,_,_,_,_ = tdcrpy.TDCR_model_lib.detectProbabilitiesMC(L, e_q, e_q2, t1, evenement, extDT, measTime, dispParam=True)
    SmcI.append(Smc); DmcI.append(Dmc); TmcI.append(Tmc)

print('\n')
tdcrpy.TDCR_model_lib.readParameters(disp=True)

print("\nEffQ = ", Q, "\tEffQ*L = ", QL, "\n")

print("\nEFF, EFFmc, +/-")
print("single eff = ",round(S,4),round(np.mean(SmcI),4),round(np.std(SmcI)/np.sqrt(nIter),4))
print("double eff = ",round(D,4),round(np.mean(DmcI),4),round(np.std(DmcI)/np.sqrt(nIter),4))
print("triple eff = ",round(T,4),round(np.mean(TmcI),4),round(np.std(TmcI)/np.sqrt(nIter),4))
print('\nDEVIATION < 2 sigma')
print("single eff = ",abs(round(S,4)-round(np.mean(SmcI),4))<2*round(np.std(SmcI)/np.sqrt(nIter),4))
print("double eff = ",abs(round(D,4)-round(np.mean(DmcI),4))<2*round(np.std(DmcI)/np.sqrt(nIter),4))
print("triple eff = ",abs(round(T,4)-round(np.mean(TmcI),4))<2*round(np.std(TmcI)/np.sqrt(nIter),4))
print('\nPRECISION')
print("single eff = ",round(100*np.std(SmcI)/(np.sqrt(nIter)*round(np.mean(SmcI),4)),4)," %")
print("double eff = ",round(100*np.std(DmcI)/(np.sqrt(nIter)*round(np.mean(DmcI),4)),4)," %")
print("triple eff = ",round(100*np.std(TmcI)/(np.sqrt(nIter)*round(np.mean(TmcI),4)),4)," %")