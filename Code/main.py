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
# N=1                   # number of simulated decay (MC trials)
N= 10000
<<<<<<< HEAD
Rad=["Sr-89"]            # list of radionuclides (Na-24)
=======
Rad=["Be-7"]            # list of radionuclides (Na-24)
>>>>>>> 419e1c69cfbda298e7934bae07a846069d4dcfaf
# Rad = ["Cs-137"]
pmf_1=[1]                # relative abondance (pmf)
# kB =[1.0e-5]
kB = [0.8e-5, 0.9e-5, 1.0e-5, 1.1e-5, 1.2e-5]    # Birks constant in cm/keV
<<<<<<< HEAD
# L=[1e-1]
L = np.logspace(-3,1,50) # Free paramete in keV-1

=======
L=[1e-1]
# L = np.logspace(-3, 0, 25) # Free paramete in keV-1 (for Cs-137)
#L = np.logspace(-3,1,50) # Free paramete in keV-1 (for Co-60)
# L = np.logspace(-3,1,50) # Free paramete in keV-1 (for Am-241)
# L = np.logspace(-3,1,50) # Free paramete in keV-1 (for Sr-90)
# L = np.logspace(-2,1,50) # Free paramete in keV-1 (for H-3, Fe-55)
>>>>>>> 419e1c69cfbda298e7934bae07a846069d4dcfaf

TDCR_measure = 0.977784        # Measured TDCR value
u_TDCR_measure = 0.000711      # standard uncertainty
Record = True                  # to record the efficiency curves
Display = False               # to display calculation results on the console
# Display = False                # to display calculation results on the console
# RHO = 0.96         #density of absorber (Toluene) g/cm3
RHO = 0.98           #density of absorber (UG + H20) g/cm3
nE = 10            #number of bin to discretize the energy vector for scintillation quenching calculation

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
DaughterVec = []     # Daughters
Pdaughter = []   # Probabiblity related to daughters
for rad_i in Rad: 
    out_PenNuc = tl.readPenNuc(rad_i)
    particle_r=[]; p_branch_r=[]; e_branch_r=[]; LevelDaughter_r=[]; levelNumber_r=[]
    prob_r=[]; levelEnergy_r=[]; transitionType_r=[]; e_trans_r=[]; next_level_r=[]
    Q_value_r=[]; DaughterVec_r=[]; Pdaughter_r=[]
    for daughter_i in range(len(out_PenNuc)):
    # for daughter_i in range(1):
        particle_r.append(out_PenNuc[daughter_i][0])          # Particle(s) from the Mother
        p_branch_r.append(out_PenNuc[daughter_i][1])          # Probablity of the branch
        e_branch_r.append(out_PenNuc[daughter_i][2])          # Energy of the branch
        LevelDaughter_r.append(out_PenNuc[daughter_i][3])     # Level of the Daughter nucleus just after the particle emission
        levelNumber_r.append(out_PenNuc[daughter_i][4])       # The vector of level of the daughter to get information of all possible isomeric transitions
        prob_r.append(out_PenNuc[daughter_i][5])              # Probibility density for each of the daughter level
        levelEnergy_r.append(out_PenNuc[daughter_i][6])       # Energy of each level
        transitionType_r.append(out_PenNuc[daughter_i][7])    # type of each possible transitions (internal transitions or gamma emission)
        e_trans_r.append(out_PenNuc[daughter_i][8])           # Energy of the transition
        next_level_r.append(out_PenNuc[daughter_i][9])        # Next level on the daughter nucleus
        
        Q_value_r.append(out_PenNuc[daughter_i][10])          # Energy of the reaction
        DaughterVec_r.append(out_PenNuc[daughter_i][11])      # Daughters
        Pdaughter_r.append(out_PenNuc[daughter_i][12])        # Probabiblity related to daughters
    particle.append(particle_r)          # Particle(s) from the Mother
    p_branch.append(p_branch_r)          # Probablity of the branch
    e_branch.append(e_branch_r)          # Energy of the branch
    LevelDaughter.append(LevelDaughter_r)     # Level of the Daughter nucleus just after the particle emission
    levelNumber.append(levelNumber_r)       # The vector of level of the daughter to get information of all possible isomeric transitions
    prob.append(prob_r)              # Probibility density for each of the daughter level
    levelEnergy.append(levelEnergy_r)       # Energy of each level
    transitionType.append(transitionType_r)    # type of each possible transitions (internal transitions or gamma emission)
    e_trans.append(e_trans_r)           # Energy of the transition
    next_level.append(next_level_r)        # Next level on the daughter nucleus
    
    Q_value.append(Q_value_r)          # Energy of the reaction
    DaughterVec.append(DaughterVec_r)      # Daughters
    Pdaughter.append(Pdaughter_r)        # Probabiblity related to daughters
               
"""
Read BetaShape
"""
e_beta = []
p_beta = []
for i, rad_i in enumerate(Rad): # radionuclide loop
    e_beta_0 = []
    p_beta_0 = []
    out_PenNuc = tl.readPenNuc(rad_i)
    for u in range(len(out_PenNuc)): # daughter loop
        e_beta_1 = []
        p_beta_1 = []
        betam = 0   # counter of beta- transition 
        betap = 0   # counter of beta+ transition
        for j in range(len(particle[i][u])): # transition loop
         if particle[i][u][j] == "beta":
             e_beta_1.append(tl.readBetaShape(rad_i, "beta-", "trans"+str(betam))[0])
             p_beta_1.append(tl.readBetaShape(rad_i, "beta-", "trans"+str(betam))[1])
             betam += 1
         elif particle[i][u][j] == "beta+":
             e_beta_1.append(tl.readBetaShape(rad_i, "beta+", "trans"+str(betap))[0])
             p_beta_1.append(tl.readBetaShape(rad_i, "beta+", "trans"+str(betap))[1])
             betap += 1
         else:
             e_beta_1.append("")
             p_beta_1.append("")   
        e_beta_0.append(e_beta_1)
        p_beta_0.append(p_beta_1)
    e_beta.append(e_beta_0)
    p_beta.append(p_beta_0)



for kB_i in kB: # Loop on the kB
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
    
        for i in range(N): # Main Loop - Monte Carlo trials
        
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
           ## sampling of the daughter
           iDaughter=tl.sampling(np.asarray(Pdaughter[index_rad])/sum(np.asarray(Pdaughter[index_rad])))
           if Display: print("\t Sampled daughter:")
           if Display: print("\t\t Daughter = ", DaughterVec[index_rad][iDaughter])           
           
           ## sampling of the decay branch
           multiplicity_branch = sum(np.asarray(p_branch[index_rad][iDaughter]))
           if p_branch[index_rad][iDaughter] != []:
               index_branch = tl.sampling(p_branch[index_rad][iDaughter])
               particle_branch = particle[index_rad][iDaughter][index_branch]            # sampled particle emitted by the mother
               energy_branch =  e_branch[index_rad][iDaughter][index_branch]             # energy of the particle emitted by the mother
               probability_branch = p_branch[index_rad][iDaughter][index_branch]         # probability of the sampled branch
               levelOftheDaughter = LevelDaughter[index_rad][iDaughter][index_branch]    # Level of the daughter just after the particle emission from the mother
               if particle_branch[:4] == "Atom":
                  particle_branch = [particle_branch]
                  particle_branch.append(DaughterVec[index_rad][iDaughter])
               if Display: print("\t Sampled decay branch:")
               if Display: print("\t\t Particle = ", particle_branch)
               if Display: print("\t\t Energy of the particle = ", energy_branch, " keV")
               if Display: print("\t\t Level of the daughter nucleus = ", levelOftheDaughter)
           
               # Scoring
               e_sum = energy_branch                               # Update the Energy summary
               particle_vec.append(particle_branch)                # Update of the particle vector
               energy_vec.append(energy_branch)                    # Update of the energy of the particle
           else:
               if Display: print("\t Sampled decay branch:")
               if Display: print("\t\t Particle = isomeric transition, no particle")
               if Display: print("\t\t Level number =", levelNumber[index_rad][iDaughter][0])
               levelOftheDaughter = Q_value[index_rad][iDaughter]
               e_sum = 0
               
               
    
           if Display: print("\t Subsequent isomeric transition")          # finish with the mother / now with the daughter
           while levelOftheDaughter > 0:                       # Go on the loop while the daughter nucleus is a its fundamental level (energy 0)
             if p_branch[index_rad][iDaughter] != []:
                 i_level = levelNumber[index_rad][iDaughter].index(levelOftheDaughter)   # Find the position in the daughter level vector
                 ## sampling of the transition in energy levels of the daughter nucleus
                 index_t = tl.sampling(prob[index_rad][iDaughter][i_level+1])  
                 if Display: print("\t\t Energy of the level = ", levelEnergy[index_rad][iDaughter][i_level], " keV")
                 if Display: print("\t\t Transition type = ", transitionType[index_rad][iDaughter][i_level+1][index_t])
                 if Display: print("\t\t Energy of the transition = ", e_trans[index_rad][iDaughter][i_level+1][index_t], "keV")
                 if Display: print("\t\t next level = ", next_level[index_rad][iDaughter][i_level+1][index_t])
                 
                 # Scoring
                 if transitionType[index_rad][iDaughter][i_level+1][index_t] == "GA": # if it is a gamma that has been emitted
                   particle_vec.append("gamma")               # Update of the particle vector
                   energy_vec.append(e_trans[index_rad][iDaughter][i_level+1][index_t])    # Update the energy vector
                 else:                                          # if not, it is a internal conversion, so an electron
                   particle_vec.append("electron")               # !!!!!!!!! it is OK for our model? Does the electron leave with the kinetic enegy of the transition 
                   energy_vec.append(e_trans[index_rad][iDaughter][i_level+1][index_t])    # Update the energy vector
                   if transitionType[index_rad][iDaughter][i_level+1][index_t] == "EK":
                      particle_vec.append(["Atom_K",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the K shell of the dughter nucleus
                      energy_vec.append(0)
                   if transitionType[index_rad][iDaughter][i_level+1][index_t] == "EL1":
                      particle_vec.append(["Atom_L1",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the L1 shell of the dughter nucleus
                      energy_vec.append(0)
                   if transitionType[index_rad][iDaughter][i_level+1][index_t] == "EL2":
                      particle_vec.append(["Atom_L2",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the L2 shell of the dughter nucleus
                      energy_vec.append(0)
                   if transitionType[index_rad][iDaughter][i_level+1][index_t] == "EL3":
                      particle_vec.append(["Atom_L3",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the L3 shell of the dughter nucleus
                      energy_vec.append(0)
                   if transitionType[index_rad][iDaughter][i_level+1][index_t] == "EM":
                      particle_vec.append(["Atom_M",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the M shell of the dughter nucleus
                      energy_vec.append(0)
                   if transitionType[index_rad][iDaughter][i_level+1][index_t] == "EN":
                      particle_vec.append(["Atom_N",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the N shell of the dughter nucleus
                      energy_vec.append(0)
                 e_sum += e_trans[index_rad][iDaughter][i_level+1][index_t]              # Energy summary
        
                 levelOftheDaughter = next_level[index_rad][iDaughter][i_level+1][index_t]   # set the next level
                 
             else:
                 i_level = levelNumber[index_rad][iDaughter][0]
                 index_t = tl.sampling(prob[index_rad][iDaughter][i_level]) 
                 # index_t = 0            # sampling of the transition
                 
                 # Scoring
                 if transitionType[index_rad][iDaughter][i_level][index_t] == "GA": # if it is a gamma that has been emitted
                   particle_vec.append("gamma")               # Update of the particle vector
                   energy_vec.append(e_trans[index_rad][iDaughter][i_level][index_t])    # Update the energy vector
                 else:                                          # if not, it is a internal conversion, so an electron
                   particle_vec.append("electron")               # !!!!!!!!! it is OK for our model? Does the electron leave with the kinetic enegy of the transition 
                   energy_vec.append(e_trans[index_rad][iDaughter][i_level][index_t])    # Update the energy vector
                   if transitionType[index_rad][iDaughter][i_level][index_t] == "EK":
                      particle_vec.append(["Atom_K",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the K shell of the dughter nucleus
                      energy_vec.append(0)
                   if transitionType[index_rad][iDaughter][i_level][index_t] == "EL1":
                      particle_vec.append(["Atom_L1",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the L1 shell of the dughter nucleus
                      energy_vec.append(0)
                   if transitionType[index_rad][iDaughter][i_level][index_t] == "EL2":
                      particle_vec.append(["Atom_L2",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the L2 shell of the dughter nucleus
                      energy_vec.append(0)
                   if transitionType[index_rad][iDaughter][i_level][index_t] == "EL3":
                      particle_vec.append(["Atom_L3",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the L3 shell of the dughter nucleus
                      energy_vec.append(0)
                   if transitionType[index_rad][iDaughter][i_level][index_t] == "EM":
                      particle_vec.append(["Atom_M",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the M shell of the dughter nucleus
                      energy_vec.append(0)
                   if transitionType[index_rad][iDaughter][i_level][index_t] == "EN":
                      particle_vec.append(["Atom_N",DaughterVec[index_rad][iDaughter]]) # record that an electron is missing on the N shell of the dughter nucleus
                      energy_vec.append(0)
                 e_sum += e_trans[index_rad][iDaughter][i_level][index_t]              # Energy summary
        
                 levelOftheDaughter = next_level[index_rad][iDaughter][i_level][index_t]-1   # set the next level

                
             
        
           # Finish with the daughter Nucleus
           if Display: print("\t Summary of the nuclear decay")
           if Display: print("\t\t particles : ", particle_vec)
           if Display: print("\t\t energy : ", energy_vec, "keV")
           # if Display: print("\t\t remaing energy : ", round(Q_value[index_rad][iDaughter]-e_sum,3), " keV")
    
    
           """
           II. LA RELAXATION ATOMIQUE
           """
           
           if Display: print("\t Summary of the atomic relaxation")
           
           lenElement = [] # pour détecter la présence de lacunes atomiques
           for element in particle_vec:
               lenElement.append(type(element))
           
           while list in lenElement:  # tant qu'il y a une lacune atomique
               lenElement = []         # pour détecter la présence de lacunes atomiques
               for element in particle_vec:
                   lenElement.append(type(element))
                   
               for i_part, part in enumerate(particle_vec):
                   if type(part) == list: # indice de la lacune dans le vecteur particle
                       tf,ef = tl.relaxation_atom(part[1],Rad[index_rad],part[0])   # tirage de la transition atomique
                       if type(tf) != int:
                           if tf[0] == 'X':                               # cas des rayons XK
                               if tf == 'XKA':                                          # cas des rayons XK_alpha
                                   particle_vec.append(["Atom_L", part[1]])             # ajout d'un lacune dans la couche L
                                   energy_vec.append(0)                                 # initialisation du vecteur energie
                                   particle_vec[i_part]='xKA'                           # mise à jour du vecteur particle avec le rayon x
                               elif tf == 'XKB':                                        # cas des rayons XK_beta
                                   particle_vec.append(["Atom_M", part[1]])             # ajout d'un lacune dans la couche M
                                   energy_vec.append(0)                                 # initialisation du vecteur energie
                                   particle_vec[i_part]='xKB'                           # mise à jour du vecteur particle avec le rayon x
                               elif tf == 'XL':
                                   particle_vec[i_part]='xL'                           # mise à jour du vecteur particle avec le rayon x
                                   # particle_vec.append(["Atom_M", part[1]])             # ajout d'un lacune dans la couche M
                                   # energy_vec.append(0)                                 # initialisation du vecteur energie
                               else:
                                   print("undetermined x rays type")
                               
                               energy_vec[i_part]=ef                                    # mise à jour du vecteur energie avec l'énergie du rayon x
                               e_sum += ef                                              # mise à jour du bilan energétique
                           if tf[0] == 'A':
                               if tf == 'Auger K':
                                   particle_vec.append(["Atom_L", part[1]])             # ajout de deux lacunes dans la couche L
                                   particle_vec.append(["Atom_L", part[1]])
                                   energy_vec.append(0)                                 # initialisation du vecteur energie
                                   energy_vec.append(0)
                                   particle_vec[i_part]='Auger K'                          # mise à jour du vecteur particle avec l'électron Auger'
                               elif tf == 'Auger L':
                                   particle_vec[i_part]='Auger L'                          # mise à jour du vecteur particle avec l'électron Auger'
                                   # particle_vec.append(["Atom_M", part[1]])             # ajout de deux lacunes dans la couche M
                                   # particle_vec.append(["Atom_M", part[1]])
                                   # energy_vec.append(0)                                 # initialisation du vecteur energie
                                   # energy_vec.append(0)         
                               else:
                                   print("undetermined Auger type")
                               energy_vec[i_part]=ef                                    # mise à jour du vecteur energie avec l'énergie de l'électron Auger
                               e_sum += ef                                              # mise à jour du bilan energétique
                       else:
                           particle_vec[i_part]="void"
    
           if Display: print("\t\t particles : ", particle_vec)            
           if Display: print("\t\t energy : ", energy_vec, "keV")
           # if Display: print("\t\t remaing energy : ", round(Q_value[index_rad][iDaughter]-e_sum,3), " keV")
               
                           
           # if me_M" in particle_vec): 
           #    print("OK")
#                for ip, p in enumerate(particle_vec):
#                  if ("Atom_K" in p) or ("Atom_L" in p) or ("Atom_M" in p):
#                     # appelle fonction() => Electron ou photon # energy
#                     # particle_vec[ip] = "electron"
#                     # energy_vec[ip] = Eout
#                     #
    
           """
           III. INTERACTION RAYONNEMENT/MATIERE + SPECTRES D'EMISSION
           """
           for i, p in enumerate(particle_vec):
             if p == "beta":
                 n_branch = len(e_branch[index_rad][iDaughter])
                 index_beta_energy = tl.sampling(p_beta[index_rad][iDaughter][-(1+index_branch)])
                 particle_vec[i] = "electron"
                 energy_vec[i] = e_beta[index_rad][iDaughter][-(1+index_branch)][index_beta_energy]
             
             if p == "beta+":
                 index_beta_energy = tl.sampling(p_beta[index_rad][iDaughter][-(1+index_branch)])
                 particle_vec[i] = "positron"
                 energy_vec[i] = e_beta[index_rad][iDaughter][-(1+index_branch)][index_beta_energy]
                 particle_vec.append("gamma")
                 particle_vec.append("gamma")
                 energy_vec.append(511)
                 energy_vec.append(511)
    
             if p == "gamma" or p == "xKA" or p == "xKB" or p == "xL":
                 energy_vec[i] = tl.energie_dep_gamma(energy_vec[i])
                 particle_vec[i] = "electron"
             
             if p == "Auger K" or p == "Auger L":
                 particle_vec[i] = "electron"
    
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
                if p == "electron" or p == "positron":
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
        
        print("\t TDCR calculation _ kB = ", kB_i, "cm/keV")
        if len(L) > 1: print("\t\t Progress = ", round(100*(L.tolist().index(L_i)+1)/len(L), 1), " %")
        print("radionuclide(s): ", Rad)
        # tl.toc()
        print("\t\t Free parameter = ", L_i, " keV-1")
        print("\t\t Efficiency of Triple coincident events = ", round(100*mean_efficiency_T[-1],3), "+/-", round(100*std_efficiency_T[-1],3), " %")
        print("\t\t Efficiency of Double coincident events = ", round(100*mean_efficiency_D[-1],3), "+/-", round(100*std_efficiency_D[-1],3), " %")
        print("\t\t Efficiency of Single events = ", round(100*mean_efficiency_S[-1],3), "+/-", round(100*std_efficiency_S[-1],3), " %")
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
    
    if len(mean_efficiency_S)>1:
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
    