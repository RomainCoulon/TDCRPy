# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:01:49 2023

A Monte-Carlo code to calculate detection efficiency in TDCR measurements

@author: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

## IMPORT PYTHON MODULES
import tdcrpy.TDCR_model_lib as tl
import importlib.resources
import configparser
import numpy as np
from tqdm import tqdm

def TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, mode, mode2, Display=False, barp=True):
    """
    This is a Monte-Carlo TDCR model

    Parameters
    ----------
    L : Float (if mode2="sym") or a tuple (if mode2="asym")
        Free parameter in keV-1.
    TD : float
        triple-to-double coincidence ratio. Not consider if mode2="asym". Not consider if mode2="asym".
    TAB : float
        triple-to-double coincidence ratio (coincidences between channel A and B). Not consider if mode2="sym".
    TBC : float
        triple-to-double coincidence ratio (coincidences between channel B and C). Not consider if mode2="sym".
    TAC : float
        triple-to-double coincidence ratio (coincidences between channel A and C). Not consider if mode2="sym".
    Rad : string
        List of radionuclides (eg. "H-3, Co-60").
    pmf_1 : string
        list of probability of each radionuclide (eg. "0.8, 0.2").
    N : integer
        Number of Monte-Carlo trials. recommanded N>10000. Not applied in the case of pure beta emitting radionuclides.
    kB : float
        Birks constant in cm/keV.
    mode : string
        "res" to return the residual, "eff" to return efficiencies.
    mode2 : string
        "sym" for symetrical model, "asym" for symetrical model.
    Display : Boolean, optional
        "True" to display details on the decay sampling. The default is False.
    barp : Boolean, optional
        "True" to display the calculation progress. The default is True.

    Returns
    -------
    Tuple
        if mode=="res", the residual (float).
        if mode=="eff", the efficiencies (list)

    """
    if barp: tl.display_header()
    config = configparser.ConfigParser()
    with importlib.resources.path('tdcrpy', 'config.toml') as data_path:
        file_conf = data_path       
    config.read(file_conf)
    Y=config["Inputs"].getboolean("Y")
    radListPureBeta=config["Inputs"].get("radListPureBeta")
    radListPureBeta=radListPureBeta.replace(" ","")
    radListPureBeta=radListPureBeta.split(',')
    X = Rad in radListPureBeta
    if X:
        nElist=config["Inputs"].get("nE")
        nElist=nElist.split(',')
        nElist = [int(i) for i in nElist]
    if X and Y:
        inE = radListPureBeta.index(Rad)
        nE = nElist[inE]
        out=tl.modelAnalytical(L,TD,TAB,TBC,TAC,Rad,kB,mode,mode2,nE)
        if mode == "res":
            return out
        if mode == "eff":
            return out[0], 0, out[1], 0, out[2], 0
    else:
        nE_electron = config["Inputs"].getint("nE_electron")
        nE_alpha = config["Inputs"].getint("nE_alpha")
        Rad=Rad.replace(" ","")
        Rad=Rad.split(",")
        pmf_1=pmf_1.split(",")
        pmf_1 = [float(x) for x in pmf_1]
        
        if np.size(pmf_1) > 1:
            if sum(pmf_1) !=1: print("warning p not equal to 1")
        elif pmf_1[0] != 1: print("warning")
    
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
    
        for rad_i in Rad:
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
    
        efficiency_S = []
        efficiency_D = []
        efficiency_T = []
        efficiency_AB = []
        efficiency_BC = []
        efficiency_AC = []
        
        if barp and not Display: NN = tqdm(range(N), desc="Processing", unit=" decays")
        else: NN = range(N)
        for i in NN: # Main Loop - Monte Carlo trials
            particle_vec=[]
            energy_vec=[]
            '''
            ===============================
            0. SAMPLING OF THE RADIONUCLIDE
            ===============================
            '''
            index_rad = tl.sampling(pmf_1)
            rad_i = Rad[index_rad]
            if Display: print("\n Sampled radionuclide = ", rad_i, "- L = ", L, ' keV-1 - kB = ', kB, ' cm/keV')
                  
            
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
            if Display: print("\t Sampled daughter:")
            if Display: print("\t\t Daughter = ", Daughter)           
            #=============================
            # Sampling of the decay branch
            #=============================
            branch_i = tl.normalise(prob_branch[index_rad][iDaughter])   # normalise la proba de branch
            i_branch=tl.sampling(branch_i)                               # indice de la branche globale
            if p_branch[index_rad][iDaughter][i_branch] != []:
                branch_proba = tl.normalise(p_branch[index_rad][iDaughter][i_branch])
                index_subBranch = tl.sampling(branch_proba)                                            # indice de la branch precise
                particle_branch = particle[index_rad][iDaughter][i_branch][index_subBranch]            # sampled particle emitted by the mother
                energy_branch =  e_branch[index_rad][iDaughter][i_branch][index_subBranch]             # energy of the particle emitted by the mother
                # probability_branch = p_branch[index_rad][iDaughter][i_branch][index_subBranch]         # probability of the sampled branch
                levelOftheDaughter = LevelDaughter[index_rad][iDaughter][i_branch][index_subBranch]    # Level of the daughter just after the particle emission from the mother
                level_before_trans = LevelDaughter[index_rad][iDaughter][i_branch][index_subBranch]
                if Display: print("\t Sampled decay branch:")
                if Display: print("\t\t Particle =               ", particle_branch)
                if Display: print("\t\t Energy of the particle = ", energy_branch, " keV")
                if Display: print("\t\t Level of the daughter nucleus = ", levelOftheDaughter)
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
                if Display: print("\t Sampled decay branch:")
                if Display: print("\t\t Particle = isomeric transition, no particle")
                if Display: print("\t\t Level of the nucleus : ",levelOftheDaughter)
                e_sum = 0
    
            '''
            ==============
            I-1 Transition
            ==============
            '''  
    
            if Display: print("\t Subsequent isomeric transition")                       # finish with the mother / now with the daughter
            while levelOftheDaughter > 0:                                                # Go on the loop while the daughter nucleus is a its fundamental level (energy 0)
                i_level = levelNumber[index_rad][iDaughter].index([levelOftheDaughter])  # Find the position in the daughter level vector
                if transitionType[index_rad][iDaughter][i_level] != []:
                    #====================================================================
                    # Sampling of the transition in energy levels of the daughter nucleus
                    #====================================================================
                    probability_tran = tl.normalise(prob_trans[index_rad][iDaughter][i_level])   # normaliser la proba de transition 
                    index_t = tl.sampling(probability_tran)                                      # indice de la transition
                    if Display: print("\t\t Energy of the level = ", levelEnergy[index_rad][iDaughter][i_level], " keV")
                    if Display: print("\t\t Transition type =          ", transitionType[index_rad][iDaughter][i_level][index_t])
                    if Display: print("\t\t Energy of the transition = ", e_trans[index_rad][iDaughter][i_level][index_t], "keV")
                    if Display: print("\t\t next level =               ", next_level[index_rad][iDaughter][i_level][index_t])
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
                    
            if Display: print("\t Summary of the nuclear decay")
            if Display: print("\t\t particles : ", particle_vec)
            if Display: print("\t\t energy    : ", energy_vec, "keV")
    
            '''
            ==========================
            II. LA RELAXATION ATOMIQUE
            ==========================
            '''
            daughter_relax = DaughterVec[index_rad][iDaughter]
            for i_part in range(len(particle_vec)):
                relaxation = False
                if "Atom_K" in particle_vec[i_part] or "Atom_L" in particle_vec[i_part] or "Atom_M" in particle_vec[i_part]:
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
                        if Display: print("untermined x or Auger")
                        relaxation = False
                    e_sum += ef
            if Display: print("\t Summary of the atomic relaxation")
            if Display: print("\t\t particles : ", particle_vec)            
            if Display: print("\t\t energy    : ", energy_vec, "keV")
    
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
                    energy_vec[i] = tl.energie_dep_beta(energy_vec[i])
    
                if p == "beta+":
                    e_b,p_b = tl.readBetaShape(rad_i,"beta+",level_before_trans)
                    index_beta_energy = tl.sampling(p_b)
                    particle_vec[i] = "positron"
                    energy_vec[i] = e_b[index_beta_energy]
                    energy_vec[i] = tl.energie_dep_beta(energy_vec[i]) # treated as beta-
                    particle_vec.append("gamma")
                    particle_vec.append("gamma")
                    energy_vec.append(511)
                    energy_vec.append(511)
    
                if p == "gamma" or p == "XKA" or p == "XKB" or p == "XL":
                    energy_vec[i] = tl.energie_dep_gamma(energy_vec[i],v=10)          # sampling energy free from photon
                    particle_vec[i] = "electron"
                if p == "Auger K" or p == "Auger L":
                    particle_vec[i] = "electron"
                    energy_vec[i] = tl.energie_dep_beta(energy_vec[i])
            if Display: print("\t Summary of the final charged particles")
            if Display: print("\t\t particles : ", particle_vec)
            if Display: print("\t\t energy after interaction : ", energy_vec, "keV")
            
            '''
            ====================
            IV. LA SCINTILLATION
            Calculation of the scintillation quenching with the Birks Model
            ====================
            '''
            if Display: print("\t Summary of the estimation of quenched energies")
            if Display: print("\t\t energy_vec      : ", energy_vec, "keV")
            e_quenching=[]
            for i, p in enumerate(particle_vec):
                if p == "alpha":
                    energy_vec[i] = tl.E_quench_a(energy_vec[i],kB,nE_alpha)
                    e_quenching.append(energy_vec[i])
                elif p == "electron" or p == "positron":
                    energy_vec[i] = tl.E_quench_e(energy_vec[i]*1e3,kB*1e3,nE_electron)*1e-3
                    e_quenching.append(energy_vec[i])
                else:
                    e_quenching.append(0)
            if Display: print("\t\t quenched energy : ", e_quenching, "keV")
    
            '''
            ====================
            V. LE MESURE TDCR
            ====================
            '''
            
            if mode2=="sym":
                p_nosingle = np.exp(-L*np.sum(np.asarray(e_quenching))/3) # probability to have 0 electrons in a PMT
                p_single = 1-p_nosingle                                    # probability to have at least 1 electrons in a PMT
                efficiency_S.append(p_single)
                efficiency_T.append(p_single**3)
                efficiency_D.append(3*(p_single)**2-2*efficiency_T[-1])
                if Display: print("\t Summary of TDCR measurement")
                if Display: print("\t\t Efficiency of single events: ", efficiency_S[-1])
                if Display: print("\t\t Efficiency of double events: ", efficiency_D[-1])
                if Display: print("\t\t Efficiency of triple events: ", efficiency_T[-1])
            elif mode2=="asym":
                pA_nosingle = np.exp(-L[0]*np.sum(np.asarray(e_quenching))/3) # probability to have 0 electrons in a PMT
                pA_single = 1-pA_nosingle                                    # probability to have at least 1 electrons in a PMT
                pB_nosingle = np.exp(-L[1]*np.sum(np.asarray(e_quenching))/3) # probability to have 0 electrons in a PMT
                pB_single = 1-pB_nosingle                                    # probability to have at least 1 electrons in a PMT
                pC_nosingle = np.exp(-L[2]*np.sum(np.asarray(e_quenching))/3) # probability to have 0 electrons in a PMT
                pC_single = 1-pC_nosingle                                    # probability to have at least 1 electrons in a PMT
                
                efficiency_AB.append(pA_single*pB_single)
                efficiency_BC.append(pB_single*pC_single)
                efficiency_AC.append(pA_single*pC_single)
                efficiency_T.append(pA_single*pB_single*pC_single)
                efficiency_D.append(efficiency_AB[-1]+efficiency_BC[-1]+efficiency_AC[-1]-2*efficiency_T[-1])
                efficiency_S.append(pA_single+pB_single+pC_single-efficiency_D[-1]-efficiency_T[-1])
                
                
    
        '''
        ====================
        VI. CALCULATION OF THE FINAL ESTIMATORS
        ====================
        '''
        mean_efficiency_T = np.mean(efficiency_T) # average
        std_efficiency_T = np.std(efficiency_T)/np.sqrt(N)   # standard deviation
        std_efficiency_T = np.sqrt(std_efficiency_T**2+1e-8) # combined with uncertainty due to quenching calculation
        mean_efficiency_D = np.mean(efficiency_D)
        std_efficiency_D = np.std(efficiency_D)/np.sqrt(N)
        std_efficiency_D = np.sqrt(std_efficiency_D**2+1e-8)
        mean_efficiency_S = np.mean(efficiency_S)
        std_efficiency_S = np.std(efficiency_S)/np.sqrt(N)
        std_efficiency_S = np.sqrt(std_efficiency_S**2+1e-8)
        if mode2=="sym":
            TDCR_calcul = mean_efficiency_T/mean_efficiency_D
        elif mode2=="asym":
            mean_efficiency_AB = np.mean(efficiency_AB)
            std_efficiency_AB = np.std(efficiency_AB)/np.sqrt(N)
            mean_efficiency_BC = np.mean(efficiency_BC)
            std_efficiency_BC = np.std(efficiency_BC)/np.sqrt(N)
            mean_efficiency_AC = np.mean(efficiency_AC)
            std_efficiency_AC = np.std(efficiency_AC)/np.sqrt(N)
            TDCR_calcul = mean_efficiency_T/mean_efficiency_D
            TABmodel = mean_efficiency_T/mean_efficiency_AB
            TBCmodel = mean_efficiency_T/mean_efficiency_BC
            TACmodel = mean_efficiency_T/mean_efficiency_AC
            
            
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
        
        if mode2=="sym":
            res=(TDCR_calcul-TD)**2
        elif mode2=="asym":
            res=(TAB-TABmodel)**2+(TBC-TBCmodel)**2+(TAC-TACmodel)**2
        
        if mode == "res":
            return res
        if mode == "eff":
            return mean_efficiency_S, std_efficiency_S, mean_efficiency_D, std_efficiency_D, mean_efficiency_T, std_efficiency_T