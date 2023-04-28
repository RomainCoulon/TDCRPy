# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:01:49 2023

A Monte-Carlo code to calculate detection efficiency in TDCR measurements

@author: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

## IMPORT PYTHON MODULES
import TDCR_model_lib as tl
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st


## INPUT OF THE MODEL
N=10000                   # number of simulated decay (MC trials)
Rad=["C-14"]            # list of radionuclides (Na-24)
pmf_1=[1]                # relative abondance (pmf)
kB = [0.8e-5, 0.9e-5, 1.0e-5, 1.1e-5, 1.2e-5]    # Birks constant in cm/keV
# L=[1e-1]
# L = np.logspace(-3, 0, 25) # Free paramete in keV-1 (for Cs-137)
L = np.logspace(-3,1,200) # Free paramete in keV-1 (for Co-60)
# L = np.logspace(-3,-1,10) # Free paramete in keV-1 (for Am-241)
# L = np.logspace(-3,1,30) # Free paramete in keV-1 (for Sr-90)
# L = np.logspace(-2,1,50) # Free paramete in keV-1 (for H-3)
TDCR_measure = 0.977784        # Measured TDCR value
u_TDCR_measure = 0.000711      # standard uncertainty
Record = True                  # to record the efficiency curves
Display = False                # to display calculation results on the console
# RHO = 0.96         #density of absorber (Toluene) g/cm3
RHO = 0.98           #density of absorber (UG + H20) g/cm3
nE = 1000            #number of bin to discretize the energy vector for scintillation quenching calculation

if np.size(pmf_1) > 1:
    if sum(pmf_1) !=1: print("warning p not equal to 1")
elif pmf_1[0] != 1: print("warning")


"""
Read PenNuc File
"""
out_PenNuc = []
particle = []
p_branch = []
e_branch = []
LevelDaughter = []
levelNumber = []
prob = []
levelEnergy = []
transitionType = []
e_trans = []
next_level = []
Q_value = []
for rad_i in Rad:
    out_PenNuc = tl.readPenNuc(rad_i)
    particle.append(out_PenNuc[0])          # Particle(s) from the Mother
    p_branch.append(out_PenNuc[1])          # Probablity of the branch
    e_branch.append(out_PenNuc[2])          # Energy of the branch
    LevelDaughter.append(out_PenNuc[3])     # Level of the Daughter nucleus just after the particle emission
    levelNumber.append(out_PenNuc[4])       # The vector of level of the daughter to get information of all possible isomeric transitions
    prob.append(out_PenNuc[5])              # Probibility density for each of the daughter level
    levelEnergy.append(out_PenNuc[6])       # Energy of each level
    transitionType.append(out_PenNuc[7])    # type of each possible transitions (internal transitions or gamma emission)
    e_trans.append(out_PenNuc[8])           # Energy of the transition
    next_level.append(out_PenNuc[9])        # Next level on the daughter nucleus
    Q_value.append(out_PenNuc[10])          # Energy of the reaction


"""
Read BetaShape
"""
e_beta = []
p_beta = []
for i, rad_i in enumerate(Rad):
    e_beta_0 = []
    p_beta_0 = []
    for j in range(len(particle[i])):
        if particle[i][j] == "beta":
            e_beta_0.append(tl.readBetaShape(rad_i, "beta-", "trans"+str(j))[0])
            p_beta_0.append(tl.readBetaShape(rad_i, "beta-", "trans"+str(j))[1])
        elif particle[i][j] == "beta+":
            e_beta_0.append(tl.readBetaShape(rad_i, "beta+", "trans"+str(j))[0])
            p_beta_0.append(tl.readBetaShape(rad_i, "beta+", "trans"+str(j))[1])
        else:
            e_beta.append("")
            p_beta.append("")   
        e_beta.append(e_beta_0)
        p_beta.append(p_beta_0)



for kB_i in kB:
    mean_efficiency_S = []  # efficiency of single counte rate
    std_efficiency_S = []   # std
    mean_efficiency_T = []  # efficiency of triple coincidence counte rate
    std_efficiency_T = []   # std
    mean_efficiency_D = []  # efficiency of double coincidence counte rate
    std_efficiency_D = []   # std
    TDCR_calcul = []
    for L_i in L: # loop on the free parameter values
        # tl.tic()
        ## RUN THE MC CALCULATION
        efficiency_S = []        # results calculated efficiency for single events
        efficiency_D = []        # results calculated efficiency for double coincidence
        efficiency_T = []        # results calculated efficiency for triple coincidence
    
        for i in range(N): # Main Loop
        
           #tl.tic()
           particle_vec=[]
           energy_vec=[]
        
           ## Sampling of the radionuclide
           index_rad = tl.sampling(pmf_1)
           rad_i = Rad[index_rad]
           if Display: print("\n Sampled radionuclide = ", rad_i, "- L = ", L_i, ' keV-1 - kB = ', kB_i, ' cm/keV')
        
           """
           I. DESINTEGRATION NUCLEAIRE
           """
                  
           ## sampling of the decay branch
           multiplicity_branch = sum(np.asarray(p_branch[index_rad]))
           index_branch = tl.sampling(p_branch[index_rad])
           particle_branch = particle[index_rad][index_branch]            # sampled particle emitted by the mother
           energy_branch =  e_branch[index_rad][index_branch]             # energy of the particle emitted by the mother
           probability_branch = p_branch[index_rad][index_branch]         # probability of the sampled branch
           levelOftheDaughter = LevelDaughter[index_rad][index_branch]    # Level of the daughter just after the particle emission from the mother
           if Display: print("\t Sampled decay branch:")
           if Display: print("\t\t Particle = ", particle_branch)
           if Display: print("\t\t Energy of the particle = ", energy_branch, " keV")
           if Display: print("\t\t Level of the daughter nucleus = ", levelOftheDaughter)
           
           # Scoring
           e_sum = energy_branch                               # Update the Energy summary
           particle_vec.append(particle_branch)                # Update of the particle vector
           energy_vec.append(energy_branch)                    # Update of the energy of the particle
    
           if Display: print("\t Subsequent isomeric transition")          # finish with the mother / now with the daughter
           while levelOftheDaughter > 0:                       # Go on the loop while the daughter nucleus is a its fundamental level (energy 0)
             i_level = levelNumber[index_rad].index(levelOftheDaughter)   # Find the position in the daughter level vector 
             ## sampling of the transition in energy levels of the daughter nucleus
             index_t = tl.sampling(prob[index_rad][i_level+1])            # sampling of the transition
             if Display: print("\t\t Energy of the level = ", levelEnergy[index_rad][i_level], " keV")
             if Display: print("\t\t Transition type = ", transitionType[index_rad][i_level+1][index_t])
             if Display: print("\t\t Energy of the transition = ", e_trans[index_rad][i_level+1][index_t], "keV")
             if Display: print("\t\t next level = ", next_level[index_rad][i_level+1][index_t])
                
             # Scoring
             if transitionType[index_rad][i_level+1][index_t] == "GA": # if it is a gamma that has been emitted
               particle_vec.append("gamma")               # Update of the particle vector
               energy_vec.append(e_trans[index_rad][i_level+1][index_t])    # Update the energy vector
             else:                                          # if not, it is a internal conversion, so an electron
               particle_vec.append("electron")               # !!!!!!!!! it is OK for our model? Does the electron leave with the kinetic enegy of the transition 
               energy_vec.append(e_trans[index_rad][i_level+1][index_t])    # Update the energy vector
               if transitionType[index_rad][i_level+1][index_t] == "EK":
                  particle_vec.append("Atom_K") # record that an electron is missing on the K shell of the dughter nucleus
                  energy_vec.append(0)
               if transitionType[index_rad][i_level+1][index_t] == "EL1":
                  particle_vec.append("Atom_L1") # record that an electron is missing on the L1 shell of the dughter nucleus
                  energy_vec.append(0)
               if transitionType[index_rad][i_level+1][index_t] == "EL2":
                  particle_vec.append("Atom_L2") # record that an electron is missing on the L2 shell of the dughter nucleus
                  energy_vec.append(0)
               if transitionType[index_rad][i_level+1][index_t] == "EL3":
                  particle_vec.append("Atom_L3") # record that an electron is missing on the L3 shell of the dughter nucleus
                  energy_vec.append(0)
               if transitionType[index_rad][i_level+1][index_t] == "EM":
                  particle_vec.append("Atom_M") # record that an electron is missing on the M shell of the dughter nucleus
                  energy_vec.append(0)
               if transitionType[index_rad][i_level+1][index_t] == "EN":
                  particle_vec.append("Atom_N") # record that an electron is missing on the N shell of the dughter nucleus
                  energy_vec.append(0)
             e_sum += e_trans[index_rad][i_level+1][index_t]              # Energy summary
    
             levelOftheDaughter = next_level[index_rad][i_level+1][index_t]   # set the next level
        
           # Finish with the daughter Nucleus
           if Display: print("\t Summary of the nuclear decay")
           if Display: print("\t\t particles : ", particle_vec)
           if Display: print("\t\t energy : ", energy_vec, "keV")
           if Display: print("\t\t remaing energy (atomic transitions) : ", Q_value[index_rad]-e_sum, " keV")
    
    
           """
           II. LA RELAXATION ATOMIQUE
           """
    
           ## Look at the EADL https://www.nndc.bnl.gov/nndc/proceedings/2010csewgusndp/Tuesday/USNDP/eadl.pdf
           ## Also BrIccEmisDB ...
    
    
           """
           III. INTERACTION RAYONNEMENT/MATIERE + SPECTRES D'EMISSION
           """
           for i, p in enumerate(particle_vec):
             if p == "beta":
                 # e_beta, p_beta, n_bin = tl.readBetaSpectrum(rad_i) # deprecated
                 # e_beta, p_beta = tl.readBetaShape(rad_i, "beta-", "tot")
                 index_beta_energy = tl.sampling(p_beta[index_rad][-(1+index_branch)])
                 particle_vec[i] = "electron"
                 energy_vec[i] = e_beta[index_rad][-(1+index_branch)][index_beta_energy]
                 # Sampling Matrice comme gamma
             
             if p == "beta+":
                 # e_beta, p_beta = tl.readBetaShape(rad_i, "beta+", "tot")
                 index_beta_energy = tl.sampling(p_beta[index_rad][-(1+index_branch)])
                 particle_vec[i] = "positron"
                 energy_vec[i] = e_beta[index_rad][-(1+index_branch)][index_beta_energy]
                 # Sampling Matrice comme gamma
    
             if p == "gamma" or p == "x":
                 energy_vec[i] = tl.energie_dep_gamma(energy_vec[i])
                 particle_vec[i] = "electron" # false Compton scattering... to develop...!!!!!!!!!!!!!!
             
             if p[:4] == "Atom": # Electron capture
                 energy_vec[i] = 0
    
           if Display: print("\t Summary of the final charged particles")
           if Display: print("\t\t particles : ", particle_vec)
           if Display: print("\t\t energy : ", energy_vec, "keV")
    
    
           # tl.tic()
           """
           IV. LA SCINTILLATION
           """
           ## Now we have the (particle, energy) vectors that we would like
    
           ## Calculation of the scintillation quenching with the Birks Model
           for i, p in enumerate(particle_vec):
                e_discrete = np.linspace(0,energy_vec[i],nE) # vector for the quenched  energy calculation keV
                delta_e = e_discrete[2]-e_discrete[1]  #keV
                if p == "alpha":
                    energy_vec[i] = np.cumsum(delta_e/(1+kB_i*tl.stoppingpowerA(e_discrete)))[-1]
                    # energy_vec[i] = 0
                    # for j in e_discrete:
                    #     energy_vec[i] += delta_e/(1+kB_i*tl.stoppingpowerA(j)) # input (keV) / output (keV)
                if p == "electron":
                    # energy_vec = np.cumsum(delta_e/(1+kB_i*1e3*tl.stoppingpower(e_discrete*1e3)))
                    energy_vec[i] = 0
                    for j in e_discrete:
                        energy_vec[i] += delta_e/(1+kB_i*1e3*tl.stoppingpower(j*1e3)) # stoppingpower :input in (eV) / output (MeV)
                    
           if Display: print("\t\t quenched energy : ", energy_vec, "keV")
    
           # tl.toc()
    
           """
           V. LE MESURE TDCR
           """
    
           ## Calculation of the TDCR ratio 
           ## We fill our 3 results vectors
           #energy_vec=[3,3,3,90,90,90,90,90,90,90]  # test against Broda ARI 58 (2003) 585-594
           #energy_vec=[90,90,90,3,3,3,3,3,3,3]  # test against Broda ARI 58 (2003) 585-594
           #energy_vec=[25,25,25,3,3,3,3,3,3,3]  # test against Broda ARI 58 (2003) 585-594
           #energy_vec=[3,3,3,3,3,3,3,3,3,3]  # test against Broda ARI 58 (2003) 585-594
           p_nosingle = np.exp(-L_i*np.sum(np.asarray(energy_vec))/3) # probability to have 0 electrons in a PMT
           p_single = 1-p_nosingle                                    # probability to have at least 1 electrons in a PMT
           efficiency_S.append(p_single)
           efficiency_T.append(p_single**3)
           efficiency_D.append(3*(p_single)**2-2*efficiency_T[-1])
           
           
           #tl.toc()
    
    
        # We calculate the final estimator
        mean_efficiency_T.append(np.mean(efficiency_T)) # average
        std_efficiency_T.append(np.std(efficiency_T)/np.sqrt(N))   # standard deviation
        mean_efficiency_D.append(np.mean(efficiency_D))
        std_efficiency_D.append(np.std(efficiency_D)/np.sqrt(N))
        mean_efficiency_S.append(np.mean(efficiency_S))
        std_efficiency_S.append(np.std(efficiency_S)/np.sqrt(N))
        TDCR_calcul.append(mean_efficiency_T[-1]/mean_efficiency_D[-1])
        
        print("\t TDCR calculation")
        if len(L) > 0: print("\t\t Progress = ", round(100*(L.tolist().index(L_i)+1)/len(L), 1), " %")
        # tl.toc()
        print("\t\t Free parameter : ", L_i, " keV-1")
        print("\t\t Efficiency of Triple coincident events : ", round(100*mean_efficiency_T[-1],3), "+/-", round(100*std_efficiency_T[-1],3), " %")
        print("\t\t Efficiency of Double coincident events : ", round(100*mean_efficiency_D[-1],3), "+/-", round(100*std_efficiency_D[-1],3), " %")
        print("\t\t Efficiency of Single events : ", round(100*mean_efficiency_S[-1],3), "+/-", round(100*std_efficiency_S[-1],3), " %")
    #    x = np.arange(np.mean(efficiency_T),1.001,0.001)
    #    plt.figure("efficiency distribution")
    #    plt.clf()
    #    plt.hist(np.asarray(efficiency_D),bins=x,label="Efficiency of double coincidences")[0]
    #    plt.hist(np.asarray(efficiency_T),bins=x,label="Efficiency of triple coincidences")[0]
    #    plt.xlabel("Efficiency", fontsize = 14)
    #    plt.ylabel(r"Number of counts", fontsize = 14)
    #    plt.legend(fontsize = 12)
    #    plt.savefig('Effdistribution.png')
    
    #    x = np.arange(np.mean(tdcr),1.001,0.001)
    #    plt.figure("TDCR distribution")
    #    plt.clf()
    #    plt.hist(np.asarray(tdcr),bins=x,label="Calculated TDCR")[0]
    #    plt.plot(x,st.norm.pdf(x, TDCR_measure, u_TDCR_measure),label="measured TDCR")[0]
    #    plt.xlabel("Efficiency", fontsize = 14)
    #    plt.ylabel(r"Number of counts", fontsize = 14)
    #    plt.legend(fontsize = 12)
    #    plt.savefig('TDCRdistribution.png')
    
    # TDCR_calcul_vec = np.asarray(efficiency_T)/np.asarray(efficiency_D)
    
    # Eff0_S_reg = tl.regress(L, mean_efficiency_S)
    # Eff0_D_reg = tl.regress(L, mean_efficiency_D)
    # Eff0_T_reg = tl.regress(L, mean_efficiency_T)
    
    plt.figure("Efficiency curve I")
    plt.clf()
    plt.title(''.join(Rad))
    plt.errorbar(L, mean_efficiency_S, yerr = std_efficiency_S, fmt=".b", label = "S")
    plt.errorbar(L, mean_efficiency_D, yerr = std_efficiency_D, fmt=".k", label = "D")
    plt.errorbar(L, mean_efficiency_T, yerr = std_efficiency_T, fmt=".r", label = "T")
    # plt.plot(Eff0_S_reg[:,0],Eff0_S_reg[:,1],'-b')
    # plt.plot(Eff0_D_reg[:,0],Eff0_D_reg[:,1],'-k')
    # plt.plot(Eff0_T_reg[:,0],Eff0_T_reg[:,1],'-r')
    plt.xscale("log")
    #plt.plot(np.mean(TDCR_calcul_vec)*np.ones(N),efficiency_T,".b")[0]
    #plt.plot([TDCR_measure, TDCR_measure], [min(mean_efficiency_D), max(mean_efficiency_D)], '-r', label="Measurement")
    plt.xlabel(r"$L$ /keV$^{-1}$", fontsize = 14)
    plt.ylabel(r"$\epsilon$", fontsize = 14)
    plt.legend(fontsize = 12)
    if Record: plt.savefig("EfficiencyCurves/"+''.join(Rad)+"/fom_"+''.join(Rad)+"_"+str(kB_i)+".png")
    plt.close()
    
    # Eff_S_reg = tl.regress(TDCR_calcul, mean_efficiency_S)
    # Eff_D_reg = tl.regress(TDCR_calcul, mean_efficiency_D)
    # Eff_T_reg = tl.regress(TDCR_calcul, mean_efficiency_T)
    
    plt.figure("Efficiency curve II")
    plt.clf()
    plt.title(''.join(Rad))
    plt.errorbar(TDCR_calcul, mean_efficiency_S, yerr=std_efficiency_S, fmt=".b", label = "S")
    plt.errorbar(TDCR_calcul, mean_efficiency_D, yerr=std_efficiency_D, fmt=".k", label = "D")
    plt.errorbar(TDCR_calcul, mean_efficiency_T, yerr=std_efficiency_T, fmt=".r", label = "T")
    # plt.plot(Eff_S_reg[:,0],Eff_S_reg[:,1],'-b')
    # plt.plot(Eff_D_reg[:,0],Eff_D_reg[:,1],'-k')
    # plt.plot(Eff_T_reg[:,0], Eff_T_reg[:,1],'-r')
    plt.xlabel(r"$\epsilon_T/\epsilon_D$", fontsize = 14)
    plt.ylabel(r"$\epsilon_D$", fontsize = 14)
    plt.legend(fontsize = 12)
    if Record: plt.savefig("EfficiencyCurves/"+''.join(Rad)+"/tdcr_"+''.join(Rad)+"_"+str(kB_i)+".png")
    plt.close()
    
    if Record:
        tl.writeEffcurves(L, mean_efficiency_S, std_efficiency_S, Rad, pmf_1, kB_i, "S")
        tl.writeEffcurves(L, mean_efficiency_D, std_efficiency_D, Rad, pmf_1, kB_i, "D")
        tl.writeEffcurves(L, mean_efficiency_T, std_efficiency_T, Rad, pmf_1, kB_i, "T")
    