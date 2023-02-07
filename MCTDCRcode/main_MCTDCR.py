# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:01:49 2023

Surogate Monte-Carlo code to calculate efficiency in TDCR measurements

@author: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

import TDCR_model_lib as tl
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st


## INPUT OF THE MODEL

N=1000 # number of simulated decay (MC trials)


Rad=["H-3", "Co-60"]    # list of radionuclides
pmf_1=[0.10, 1]         # relative abondance (pmf)

kB = 0.01                  # Birks constant in cm/MeV
L = 0.0056                  # the free parameter /MeV-1

TDCR_measure = 0.977784        # Measured TDCR value
u_TDCR_measure = 0.000711

## RUN THE MC CALCULATION
efficiency_D = []
efficiency_T = []
tdcr = []

for i in range(N): # Main Loop
    particle_vec=[]
    energy_vec=[]
    
    ## Sampling of the radionuclide
    index_rad = tl.sampling(pmf_1)
    rad_i = Rad[index_rad]
    print("\n Sampled radionuclide = ", rad_i)
    
    # Read PenNuc File
    out_PenNuc = tl.readPenNuc(rad_i)
    particle = out_PenNuc[0]
    p_branch = out_PenNuc[1]
    e_branch = out_PenNuc[2]
    LevelDaughter = out_PenNuc[3]
    levelNumber = out_PenNuc[4]
    prob = out_PenNuc[5]
    levelEnergy = out_PenNuc[6]
    transitionType = out_PenNuc[7]
    e_trans = out_PenNuc[8] 
    next_level = out_PenNuc[9] 
    Q_value = out_PenNuc[10]
    
    ## sampling of the decay branch
    multiplicity_branch = sum(np.asarray(p_branch))
    index_branch = tl.sampling(p_branch)
    particle_branch = particle[index_branch]
    energy_branch =  e_branch[index_branch]
    probability_branch = p_branch[index_branch]
    levelOftheDaughter = LevelDaughter[index_branch]
    print("\t Sampled decay branch:")
    print("\t\t Particle = ", particle_branch)
    print("\t\t Energy of the particle = ", energy_branch, " keV")
    print("\t\t Level of the daughter nucleus = ", levelOftheDaughter)

    # Scoring
    e_sum = energy_branch # Energy summary
    particle_vec.append(particle_branch)
    energy_vec.append(energy_branch)
    
    print("\t Subsequent isomeric transition")
    while levelOftheDaughter > 0:
      i_level = levelNumber.index(levelOftheDaughter)
      ## sampling of the transition in energy levels of the daughter nucleus
      index_t = tl.sampling(prob[i_level+1])
      print("\t\t ***")
      print("\t\t Energy of the level = ", levelEnergy[i_level], " keV")
      print("\t\t Transition type = ", transitionType[i_level+1][index_t])
      print("\t\t Energy of the transition = ", e_trans[i_level+1][index_t], "keV")
      print("\t\t next level = ", next_level[i_level+1][index_t])
      levelOftheDaughter = next_level[i_level+1][index_t]
      
      # Scoring
      if transitionType[i_level+1][index_t] == "GA":
          particle_vec.append("gamma")
      else:
          particle_vec.append("electron")
      energy_vec.append(e_trans[i_level+1][index_t])
      e_sum += e_trans[i_level+1][index_t] # Energy summary
    
    
    print("\t Summary of the nuclear decay")
    print("\t\t particles : ", particle_vec)
    print("\t\t energy : ", energy_vec, "keV")
    print("\t\t remaing energy (atomic transitions) : ", Q_value-e_sum, " keV")

    ## Atomic transitions to develop...!!!!!!!!!!
    
    
    ## Calculation of deposited energy
    
    for i, p in enumerate(particle_vec):
        if p == "beta":
            e_beta, p_beta, n_bin = tl.readBetaSpectrum(rad_i)
            index_beta_energy = tl.sampling(p_beta)
            particle_vec[i] = "electron"
            energy_vec[i] = e_beta[index_beta_energy]*1e-3
        if p == "gamma":
            particle_vec[i] = "electron" # false Compton scattering... to develop...!!!!!!!!!!!!!!
               
    print("\t Summary of the final charged particles")
    print("\t\t particles : ", particle_vec)
    print("\t\t energy : ", energy_vec, "keV")
 
    ## Calculation of the scintillation quenching
    
    for i, p in enumerate(particle_vec):
        if p == "alpha":
            energy_vec[i] = energy_vec[i]/(1+kB*tl.stoppingpowerA(energy_vec[i]*1e3)) # false to correct
        if p == "electron":
            energy_vec[i] = energy_vec[i]/(1+kB*tl.stoppingpowerE(energy_vec[i]*1e3))
    
    print("\t\t quenched energy : ", energy_vec, "keV")
    
    ## Calculation of the TDCR ratio
    efficiency_T.append((1-np.exp(-L*sum(energy_vec)/3))**3)
    efficiency_D.append(3*(1-np.exp(-L*sum(energy_vec)/3))**2-2*efficiency_T[-1])
    tdcr.append(efficiency_T[-1]/efficiency_D[-1])

mean_efficiency_T = np.mean(efficiency_T)
std_efficiency_T = np.std(efficiency_T)
mean_efficiency_D = np.mean(efficiency_D)
std_efficiency_D = np.std(efficiency_D)
mean_TDCR = np.mean(tdcr)
std_TDCR = np.std(tdcr)
#std_TDCR_calcul = TDCR_calcul*np.sqrt((std_efficiency_D/mean_efficiency_D)**2+(std_efficiency_T/mean_efficiency_T)**2)/np.sqrt(N) # no correlation false !!!

residual = (mean_TDCR-TDCR_measure)**2 # Calculation of the residual value to minimize
u_residual = np.sqrt(4*(mean_TDCR+TDCR_measure)**2*(std_TDCR**2+u_TDCR_measure**2))
Z_test = abs(mean_TDCR-TDCR_measure)/np.sqrt(std_TDCR**2+u_TDCR_measure**2)

print("\t TDCR calculation")
print("\t\t Efficiency of Triple coincident events : ", round(100*mean_efficiency_T,3), "+/-", round(100*std_efficiency_T,3), " %")
print("\t\t Efficiency of Double coincident events : ", round(100*mean_efficiency_D,3), "+/-", round(100*std_efficiency_D,3), " %")
print("\t\t Calculated TDCR value : ", round(mean_TDCR,5), "+/-", round(std_TDCR,5))
print("\t\t Measured TDCR value : ", round(TDCR_measure,5), "+/-", round(u_TDCR_measure,5))
print("\t\t Residual : ", residual, '+/-', u_residual)
print("Z_test = ", Z_test)
    

x = np.arange(np.mean(efficiency_T),1.001,0.001)
plt.figure("efficiency distribution")
plt.clf()
plt.hist(np.asarray(efficiency_D),bins=x,label="Efficiency of double coincidences")[0]
plt.hist(np.asarray(efficiency_T),bins=x,label="Efficiency of triple coincidences")[0]
plt.xlabel("Efficiency", fontsize = 14)
plt.ylabel(r"Number of counts", fontsize = 14)
plt.legend(fontsize = 12)

x = np.arange(np.mean(tdcr),1.001,0.001)
plt.figure("TDCR distribution")
plt.clf()
plt.hist(np.asarray(tdcr),bins=x,label="Calculated TDCR")[0]
plt.plot(x,st.norm.pdf(x, TDCR_measure, u_TDCR_measure),label="measured TDCR")[0]
plt.xlabel("Efficiency", fontsize = 14)
plt.ylabel(r"Number of counts", fontsize = 14)
plt.legend(fontsize = 12)



