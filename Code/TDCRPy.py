# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 16:49:09 2023

@author: romain.coulon
"""


## IMPORT PYTHON MODULES
import TDCR_model_lib as tl
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import scipy.optimize as opt

def TDCRPy(L,TD,Rad,pmf_1,N,kB,RHO,nE):
    
    """
    Read PenNuc File
    """
    out_PenNuc = []   
    particle = []          # Particle(s) from the Mother  --  indice 3
    p_branch = []          # Probablity of the different decay of branch -- indice 5
    e_branch = []          # Energy of the different decay of branch -- indice 4
    LevelDaughter = []     # Level of the Daughter nucleus just after the particle emission -- indice 6
    levelNumber = []       # The vector of level of the daughter to get information of all possible isomeric transitions -- indice 11
    prob_trans = []        # Probability for each transition -- indice 10
    prob_branch = []       # Probability for each branch -- indice 7
    levelEnergy = []       # Energy of each level -- indice 13
    transitionType = []    # type of each possible transitions (internal transitions or gamma emission) -- indice 8
    e_trans = []           # Energy of the transition -- indice 9
    next_level = []        # Next level on the daughter nucleus -- indice 12
    Q_value = []           # Energy of the reaction -- indice 2
    DaughterVec = []       # Daughters -- indice 0
    Pdaughter = []         # Probabiblity related to daughters -- indice 1
    Transition_prob_sum = []

    meta_vec = ["Am-242m","Pa-234m","Pm-148m","Pr-144m","Xe-133m","Te-127m","Ag-110m","Ag-108m","Tc-99m","Nb-95m","Y-90m","Mn-52m"]

    for rad_i in [Rad]:
        out_PenNuc = tl.readPenNuc2(rad_i)
        particle.append(out_PenNuc[3])       # Particle(s) from the Mother  --  indice 3
        p_branch.append(out_PenNuc[5])       # Probablity of the different decay of branch -- indice 5
        e_branch.append(out_PenNuc[4])       # Energy of the different decay of branch -- indice 4
        LevelDaughter.append(out_PenNuc[6])  # Level of the Daughter nucleus just after the particle emission -- indice 6
        levelNumber.append(out_PenNuc[11])   # The vector of level of the daughter to get information of all possible isomeric transitions -- indice 11
        prob_trans.append(out_PenNuc[10])    # Probability for each transition -- indice 10
        prob_branch.append(out_PenNuc[7])    # Probability for each branch -- indice 7
        levelEnergy.append(out_PenNuc[13])   # Energy of each level -- indice 13
        transitionType.append(out_PenNuc[8]) # type of each possible transitions (internal transitions or gamma emission) -- indice 8
        e_trans.append(out_PenNuc[9])        # Energy of the transition -- indice 9
        next_level.append(out_PenNuc[12])    # Next level on the daughter nucleus -- indice 12
        Q_value.append(out_PenNuc[2])        # Energy of the reaction -- indice 2
        DaughterVec.append(out_PenNuc[0])    # Daughters -- indice 0
        Pdaughter.append(out_PenNuc[1])      # Probabiblity related to daughters -- indice 1
        Transition_prob_sum.append(out_PenNuc[14])
    

    ## RUN THE MC CALCULATION
    efficiency_S = []        # results calculated efficiency for single events
    efficiency_D = []        # results calculated efficiency for double coincidence
    efficiency_T = []        # results calculated efficiency for triple coincidence
    
    for i in range(N): # Main Loop - Monte Carlo trials
    
       #tl.tic()
        particle_vec=[]
        energy_vec=[]
    
       #==============================
        # Sampling of the radionuclide
       #==============================
        index_rad = tl.sampling(pmf_1)
        rad_i = [Rad][index_rad]
        
        '''
        ===========================
        I. DESINTEGRATION NUCLEAIRE
        ===========================
        '''
    
        #=========================
        # Sampling of the daughter
        #=========================
        iDaughter=tl.sampling(np.asarray(Pdaughter[index_rad])/sum(np.asarray(Pdaughter[index_rad])))
        Daughter = DaughterVec[index_rad][iDaughter]
     
        #=============================
        # Sampling of the decay branch
        #=============================
        branch_i = tl.normalise(prob_branch[index_rad][iDaughter])   # normalise la proba de branch
        i_branch=tl.sampling(branch_i)                               # indice de la branche globale
        #if Display: print("132 branch:",prob_branch[index_rad][iDaughter])
    
        if p_branch[index_rad][iDaughter][i_branch] != []:
            branch_proba = tl.normalise(p_branch[index_rad][iDaughter][i_branch])
            index_subBranch = tl.sampling(branch_proba)                                            # indice de la branch precise
            particle_branch = particle[index_rad][iDaughter][i_branch][index_subBranch]            # sampled particle emitted by the mother
            energy_branch =  e_branch[index_rad][iDaughter][i_branch][index_subBranch]             # energy of the particle emitted by the mother
            probability_branch = p_branch[index_rad][iDaughter][i_branch][index_subBranch]         # probability of the sampled branch
            levelOftheDaughter = LevelDaughter[index_rad][iDaughter][i_branch][index_subBranch]    # Level of the daughter just after the particle emission from the mother
            level_before_trans = LevelDaughter[index_rad][iDaughter][i_branch][index_subBranch]

            
            #========
            # Scoring
            #========
            e_sum = energy_branch                               # Update the Energy summary
            particle_vec.append(particle_branch)                # Update of the particle vector
            energy_vec.append(energy_branch)                    # Update of the energy of the particle
        else:
            transition_prob = tl.normalise(Transition_prob_sum[index_rad][iDaughter])
            index_transition_level = tl.sampling(transition_prob)
            levelOftheDaughter = levelNumber[index_rad][iDaughter][index_transition_level][0]
            e_sum = 0
    
        '''
        ==============
        I-1 Transition
        ==============
        '''  
    
        while levelOftheDaughter > 0:                                                # Go on the loop while the daughter nucleus is a its fundamental level (energy 0)
            i_level = levelNumber[index_rad][iDaughter].index([levelOftheDaughter])  # Find the position in the daughter level vector
            if transitionType[index_rad][iDaughter][i_level] != []:
                
                #====================================================================
                # Sampling of the transition in energy levels of the daughter nucleus
                #====================================================================

                probability_tran = tl.normalise(prob_trans[index_rad][iDaughter][i_level])   # normaliser la proba de transition 
                index_t = tl.sampling(probability_tran)                                      # indice de la transition

                #========
                # Scoring
                #========
    
                if transitionType[index_rad][iDaughter][i_level][index_t] == "GA":            # if it is a gamma that has been emitted
                    particle_vec.append("gamma")                                              # Update of the particle vector
                    energy_vec.append(e_trans[index_rad][iDaughter][i_level][index_t])        # Update the energy vector
                else:                                                                         # if not, it is a internal conversion, so an electron
                    particle_vec.append("electron")                                           # !!!!!!!!! it is OK for our model? Does the electron leave with the kinetic enegy of the transition 
                    energy_vec.append(e_trans[index_rad][iDaughter][i_level][index_t])        # Update the energy vector
                    if transitionType[index_rad][iDaughter][i_level][index_t] == "EK":        # record that an electron is missing on the K shell of the dughter nucleus
                        particle_vec.append("Atom_K")
                        energy_vec.append(0)
    
                    if transitionType[index_rad][iDaughter][i_level][index_t] == "EL":       # record that an electron is missing on the L1 shell of the dughter nucleus
                        particle_vec.append("Atom_L")
                        energy_vec.append(0)
    
                    if transitionType[index_rad][iDaughter][i_level][index_t] == "EL1":       # record that an electron is missing on the L1 shell of the dughter nucleus
                        particle_vec.append("Atom_L1")
                        energy_vec.append(0)
    
                    if transitionType[index_rad][iDaughter][i_level][index_t] == "EL2":       # record that an electron is missing on the L2 shell of the dughter nucleus
                        particle_vec.append("Atom_L2")
                        energy_vec.append(0)
    
                    if transitionType[index_rad][iDaughter][i_level][index_t] == "EL3":       # record that an electron is missing on the L3 shell of the dughter nucleus
                        particle_vec.append("Atom_L3")
                        energy_vec.append(0)
    
                    if transitionType[index_rad][iDaughter][i_level][index_t] == "EM":        # record that an electron is missing on the M shell of the dughter nucleus
                        particle_vec.append("Atom_M")
                        energy_vec.append(0)
    
                    if transitionType[index_rad][iDaughter][i_level][index_t] == "EN":        # record that an electron is missing on the N shell of the dughter nucleus
                        particle_vec.append("Atom_N")
                        energy_vec.append(0)
    
                e_sum += e_trans[index_rad][iDaughter][i_level][index_t]                      # Energy summary
    
                levelOftheDaughter = next_level[index_rad][iDaughter][i_level][index_t]       # set the next level
             
            else:
                i_level = levelNumber[index_rad][iDaughter].index([levelOftheDaughter])
                print("warning:pas de données de transition:daughter,niveau,niveau d'énergie",DaughterVec[index_rad][iDaughter],levelOftheDaughter,levelEnergy[index_rad][iDaughter][i_level] )
                levelOftheDaughter = 0   # set the next level
    
        '''
        ==========================
        II. LA RELAXATION ATOMIQUE
        ==========================
        '''
      
        daughter_relax = DaughterVec[index_rad][iDaughter]
        for i_part in range(len(particle_vec)):
            relaxation = False
            if "Atom_K" in particle_vec[i_part] or "Atom_L" in particle_vec[i_part]:
                relaxation = True
            while relaxation:
                tf,ef = tl.relaxation_atom(daughter_relax,Rad[index_rad],particle_vec[i_part])
                if tf == "XKA":
                    particle_vec[i_part] = "Atom_L"
                    particle_vec.append(tf)
                    energy_vec.append(ef)
                    relaxation = True
                elif tf == "XKB":
                    particle_vec[i_part] = "Atom_M"
                    particle_vec.append(tf)
                    energy_vec.append(ef)
                    relaxation = False
                elif tf == "XL":
                    particle_vec[i_part] = "Atom_M"
                    particle_vec.append(tf)
                    energy_vec.append(ef)
                    relaxation = False
                elif tf == "Auger K":
                    particle_vec[i_part] = "Atom_L"
                    particle_vec.append(tf)
                    energy_vec.append(ef)
                    relaxation = True
                elif tf == "Auger L":
                    particle_vec[i_part] = "Atom_M"
                    particle_vec.append(tf)
                    energy_vec.append(ef)
                    relaxation = False
                else:
                    relaxation = False
                e_sum += ef
      
        '''
        ==========================================================
        III. INTERACTION RAYONNEMENT/MATIERE + SPECTRES D'EMISSION
        ==========================================================
        '''
        for i, p in enumerate(particle_vec):
            if p == "beta":
                e_b,p_b = tl.readBetaShape(rad_i,"beta-",level_before_trans)   # read the data of BetaShape
                index_beta_energy = tl.sampling(p_b)                           # sampling energy of beta
                particle_vec[i] = "electron"
                energy_vec[i] = e_b[index_beta_energy]
    
            if p == "beta+":
                e_b,p_b = tl.readBetaShape(rad_i,"beta+",level_before_trans)
                index_beta_energy = tl.sampling(p_b)
                particle_vec[i] = "positron"
                energy_vec[i] = e_b[index_beta_energy]
                particle_vec.append("gamma")
                particle_vec.append("gamma")
                energy_vec.append(511)
                energy_vec.append(511)
    
            if p == "gamma" or p == "XKA" or p == "XKB" or p == "XL":
                energy_vec[i] = tl.energie_dep_gamma(energy_vec[i])          # sampling energy free from photon
                particle_vec[i] = "photon"
            if p == "Auger K" or p == "Auger L":
                particle_vec[i] = "electron"
    
        '''
        ====================
        IV. LA SCINTILLATION
        ====================
        '''
   
       ## Calculation of the scintillation quenching with the Birks Model
        e_quenching=[]
        for i, p in enumerate(particle_vec):
            e_discrete = np.linspace(0,energy_vec[i],nE) # vector for the quenched  energy calculation keV
            delta_e = e_discrete[2]-e_discrete[1]  #keV
            if p == "alpha":
                energy_vec[i] = np.cumsum(delta_e/(1+kB*tl.stoppingpowerA(e_discrete)))[-1]
                e_quenching.append(energy_vec[i])
            elif p == "electron" or p == "positron":
                energy_vec[i] = tl.E_quench_e(energy_vec[i]*1e3,kB*1e3,nE)*1e-3
                e_quenching.append(energy_vec[i])
            else:
                e_quenching.append(0)   
    
    
        '''
        V. LE MESURE TDCR
        '''
        p_nosingle = np.exp(-L*np.sum(np.asarray(energy_vec))/3) # probability to have 0 electrons in a PMT
        p_single = 1-p_nosingle                                    # probability to have at least 1 electrons in a PMT
        efficiency_S.append(p_single)
        efficiency_T.append(p_single**3)
        efficiency_D.append(3*(p_single)**2-2*efficiency_T[-1])
       
    # We calculate the final estimator
    mean_efficiency_T=np.mean(efficiency_T) # average
    std_efficiency_T=np.std(efficiency_T)/np.sqrt(N)   # standard deviation
    mean_efficiency_D=np.mean(efficiency_D)
    std_efficiency_D=np.std(efficiency_D)/np.sqrt(N)
    mean_efficiency_S=np.mean(efficiency_S)
    std_efficiency_S=np.std(efficiency_S)/np.sqrt(N)
    TDCR_calcul=mean_efficiency_T/mean_efficiency_D
    
    RES=(TDCR_calcul-TD)**2
      
    return RES



## MESUREMENT DATA
TDCR_measure = 0.977784        # Measured TDCR value
u_TDCR_measure = 0.000711      # standard uncertainty

L=2.2
TD=TDCR_measure
Rad="H-3"
pmf_1=1
N=100
kB=1e-5
RHO=0.98
nE=1000


r=opt.minimize(TDCRPy, L, args=(TD,Rad,pmf_1,N,kB,RHO,nE), method='nelder-mead',options={'xatol': 1e-5, 'disp': True, 'maxiter':50})
L=r.x
print(L)

