# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:01:49 2023

A Monte-Carlo code to calculate detection efficiency in TDCR measurements

@author: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

## IMPORT PYTHON MODULES
import tdcrpy.TDCR_model_lib as tl
# import TDCR_model_lib as tl
# import tdcrpy.TDCR_model_lib as tl
import importlib.resources
from importlib.resources import files
import configparser
import numpy as np
from tqdm import tqdm

def relaxAtom(daughter_relax,particle_vec,energy_vec,rad,Display=False,uncData=False):
         
    for i_part in range(len(particle_vec)):
        relaxation = False
        if "Atom_K" in particle_vec[i_part] or "Atom_L" in particle_vec[i_part] or "Atom_M" in particle_vec[i_part]:
            relaxation = True
        while relaxation:
            tf,ef = tl.relaxation_atom(daughter_relax,rad,particle_vec[i_part],uncData=uncData)
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
            elif tf == "Auger KLL":
                particle_vec[i_part] = "Atom_L"
                tf1,ef1 = tl.relaxation_atom(daughter_relax,rad,'Atom_L',uncData=uncData)
                particle_vec.append(tf)
                energy_vec.append(ef)
                particle_vec.append(tf1)
                energy_vec.append(ef1)
                if tf1 == 'Auger L':
                    particle_vec.append("Atom_M")
                    particle_vec.append("Atom_M")
                    energy_vec.append(0)
                    energy_vec.append(0)
                else:
                    particle_vec.append("Atom_M")
                    energy_vec.append(0)
                relaxation = True
            elif tf == "Auger KLX":
                particle_vec[i_part] = "Atom_L"
                particle_vec.append("Atom_M")
                particle_vec.append(tf)
                energy_vec.append(0)
                energy_vec.append(ef)
                relaxation = True    
            elif tf == "Auger KXY":
                particle_vec[i_part] = "Atom_M"
                particle_vec.append("Atom_M")
                particle_vec.append(tf)
                energy_vec.append(0)
                energy_vec.append(ef)
                relaxation = False    
            elif tf == "Auger L":
                particle_vec[i_part] = "Atom_M"
                particle_vec.append("Atom_M")
                particle_vec.append(tf)
                energy_vec.append(0)
                energy_vec.append(ef)
                relaxation = False
            else:
                if Display: print(f"\t\t untermined x or Auger = {tf}")
                relaxation = False
    return particle_vec, energy_vec

def TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V, mode, mode2, Display=False, barp=False,uncData=False):
    """
    This is the main function of the TDCRPy package running the Monte-Carlo Triple-to-Double Coincidence Ratio model.
    The computation is made for a given solution containing a radionuclide (or a mixture of radionuclides), a given volume of scintillator V and a given Birks constant kB. 
    
    It can operates in two modes:
       
       --> In mode="eff", it calculates the efficiency of the TDCR system as a function of a value (triplet) of free parameter(s) L, the measurement data is not used;
       
       --> In mode="res", it calculates the residual of the TDCR model parametrized by a value (or triplet) of free parameter(s) L and the measurement data TD, TAB, TBC, TAC.
    
    also, two configuration can be set:
       
        --> mode2="sym", where symmetry is considered between the 3 photomultiplier tubes - here L is a scalar and only the global TDCR value TD is used as measurement data.
        
        --> mode2="asym", where an asymmetry between the 3 photomultiplier tubes is possible - here L is a triplet and only the specific TDCR values TAB, TBC, TAC are used as measurement data.
    
    The parmeter N sets the number of Monte-Carlo trails used for the estimation. Each MC trial corresponds to a simulated radiactive decay.
    TDCRPY() used a set of fonctions from the tdcrpy.TDCR_model_lib module.
    
    Advanced settings can be configured in the config.toml file.
    
       --> By default Y = True so that the analytical model is applied for solution containing only pure beta emitting radionuclides. If you would like to apply the MC calculation also for these nuclides, set Y = False.
       
       --> If you would like to change the number of bins nE to discretize the linear energy space for quenching calculation, you can change nE_electron and nE_alpha parameters for respectively electrons and alpha particles.
       
       --> By default the calculation is set for Ultima-Gold cocktail mixed with a small amount of aqueous solution. You can adapt for a specific scintillator by changing the density, the mean charge number Z and the mean mass number A of the scintillator.
       
    
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
        Number of Monte-Carlo trials. recommanded N>10000 (see JCGM 101). Not applied in the case of pure beta emitting radionuclides.
    kB : float
        Birks constant in cm/keV.
    V : float
        volume of the scintillator in ml.
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
    res : float
        Residuals of the model compared the measurement data for (a) given free parmeters L. (only in mode="res")
    mean_efficiency_S : float
        Estimation of the efficiency of single counting events. (only in mode="eff")
    std_efficiency_S : float
        Standard uncertainty from calculation associated with the estimation of the efficiency of single counting events. (only in mode="eff")
    mean_efficiency_D : float
        Estimation of the efficiency of logic sum of double coincidences. (only in mode="eff")
    std_efficiency_D : float
        Standard uncertainty from calculation associated with the estimation of the efficiency of logic sum of double coincidences. (only in mode="eff")
    mean_efficiency_T : float
        Estimation of the efficiency of triple coincidences. (only in mode="eff")
    std_efficiency_T : float
        Standard uncertainty from calculation associated with the estimation of the efficiency of triple coincidences. (only in mode="eff")    
    """
    if barp: tl.display_header()
    config = configparser.ConfigParser()
    with importlib.resources.as_file(files('tdcrpy').joinpath('config.toml')) as data_path:
        file_conf = data_path       
    config.read(file_conf)
    tau=config["Inputs"].getfloat("tau")
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
        print(f"Analytical model used for {Rad}")
        out=tl.modelAnalytical(L,TD,TAB,TBC,TAC,Rad,kB,V,mode,mode2,nE)
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
        u_prob_trans = []        # uncertainty of the probabilities for each transition -- indice 10
        trans_halfLife = []      # half life of the transition
    
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
            u_prob_trans.append(out_PenNuc[16])
            trans_halfLife.append(out_PenNuc[15])
        # print("\n",trans_halfLife)
        
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
            if Display: print("\n\n Trial ",str(i+1),"- Sampled radionuclide: ", rad_i)
                  
            
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
            if Display:
                print("\t Sampled daughter:")
                print("\t\t Daughter = ", Daughter)           
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
                if Display:
                    print("\t Sampled decay branch:")
                    if particle_branch[:4]=="Atom":
                        if particle_branch=="Atom_K": print("\t\t Electron capture on K shell")
                        if particle_branch=="Atom_L": print("\t\t Electron capture on L shell")
                        if particle_branch=="Atom_L1": print("\t\t Electron capture on L1 shell")
                        if particle_branch=="Atom_L2": print("\t\t Electron capture on L2 shell")
                        if particle_branch=="Atom_L3": print("\t\t Electron capture on L3 shell")
                        if particle_branch=="Atom_M": print("\t\t Electron capture on M shell")
                        if particle_branch=="Atom_N": print("\t\t Electron capture on N shell")
                        if particle_branch=="Atom_O": print("\t\t Electron capture on O shell")
                    else:
                        print("\t\t Particle: ", particle_branch)
                        print("\t\t Energy of the branch transition = ", energy_branch, " keV")
                    print("\t\t Level of the daughter nucleus: ", levelOftheDaughter)
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
                if Display:
                    print("\t Sampled decay branch:")
                    print("\t\t Particle = isomeric transition, no particle")
                    print("\t\t Level of the nucleus : ",levelOftheDaughter)
                e_sum = 0
    
            '''
            ==============
            I-1 Transition
            ==============
            '''  
            if Display: print("\t Subsequent isomeric transition(s)")                       # finish with the mother / now with the daughter
            evenement = 1
            e_sum2 = 0
            particle_vec2 = []
            energy_vec2 = []
            while levelOftheDaughter > 0:                                                # Go on the loop while the daughter nucleus is a its fundamental level (energy 0)
                i_level = levelNumber[index_rad][iDaughter].index([levelOftheDaughter])  # Find the position in the daughter level vector
                
                t1 = np.random.exponential(trans_halfLife[index_rad][iDaughter][i_level][0], size=1)[0]

                # test whether the decay occurs within the coincidence resolving time or not
                if t1 > tau*1e-9: 
                    evenement = evenement + 1
                    if Display: print(f"\t\t Transition time from decay {round(t1*1e9,2)} ns > {tau} ns \n\t\t (half-life = {round(trans_halfLife[index_rad][iDaughter][i_level][0]*1e9,2)} ns)")
                                
                if transitionType[index_rad][iDaughter][i_level] != []:
                    #====================================================================
                    # Sampling of the transition in energy levels of the daughter nucleus
                    #====================================================================
                    
                    if uncData:                                                           # uncertainty
                        prob_trans_s=[]
                        for ipt, xpt in enumerate(prob_trans[index_rad][iDaughter][i_level]):
                            prob_trans_s.append(np.random.normal(xpt, u_prob_trans[index_rad][iDaughter][i_level][ipt], 1)[0])
                            
                        probability_tran = tl.normalise(prob_trans_s)   # normaliser la proba de transition
                    else:
                        probability_tran = tl.normalise(prob_trans[index_rad][iDaughter][i_level])   # normaliser la proba de transition 
                
                    index_t = tl.sampling(probability_tran)                                          # indice de la transition
                    if Display:
                        print("\t\t Energy of the level = ", levelEnergy[index_rad][iDaughter][i_level][0], " keV")
                        trans = transitionType[index_rad][iDaughter][i_level][index_t]
                        if trans == "GA":
                            print("\t\t Energy of the gamma ray = ", e_trans[index_rad][iDaughter][i_level][index_t], "keV")
                        elif trans == "EK":
                            print("\t\t Energy of the conversion electron from K shell = ", e_trans[index_rad][iDaughter][i_level][index_t], "keV")
                        elif trans == "EL":
                            print("\t\t Energy of the conversion electron from L shell = ", e_trans[index_rad][iDaughter][i_level][index_t], "keV")
                        elif trans == "EL1":
                            print("\t\t Energy of the conversion electron from L1 shell = ", e_trans[index_rad][iDaughter][i_level][index_t], "keV") 
                        elif trans == "EL2":
                            print("\t\t Energy of the conversion electron from L2 shell = ", e_trans[index_rad][iDaughter][i_level][index_t], "keV")
                        elif trans == "EL3":
                            print("\t\t Energy of the conversion electron from L3 shell = ", e_trans[index_rad][iDaughter][i_level][index_t], "keV")    
                        elif trans == "EM":
                            print("\t\t Energy of the conversion electron from M shell = ", e_trans[index_rad][iDaughter][i_level][index_t], "keV")
                        elif trans == "EN":
                            print("\t\t Energy of the conversion electron from N shell = ", e_trans[index_rad][iDaughter][i_level][index_t], "keV")
                        elif trans == "EO":
                            print("\t\t Energy of the conversion electron from O shell = ", e_trans[index_rad][iDaughter][i_level][index_t], "keV")     
                        print("\t\t next level = ", next_level[index_rad][iDaughter][i_level][index_t])
                   
                    #========
                    # Scoring
                    #========
                    
                    ## evenement retardé
                    
                    if evenement != 1:
                        if transitionType[index_rad][iDaughter][i_level][index_t] == "GA":             # if it is a gamma that has been emitted
                            particle_vec2.append("gamma")                                              # Update of the particle vector
                            energy_vec2.append(e_trans[index_rad][iDaughter][i_level][index_t])        # Update the energy vector
                        else:                                                                          # if not, it is a internal conversion, so an electron
                            particle_vec2.append("electron")                                           # !!!!!!!!! it is OK for our model? Does the electron leave with the kinetic enegy of the transition 
                            energy_vec2.append(e_trans[index_rad][iDaughter][i_level][index_t])        # Update the energy vector
                            if transitionType[index_rad][iDaughter][i_level][index_t] == "EK":         # record that an electron is missing on the K shell of the dughter nucleus
                                particle_vec2.append("Atom_K")
                                energy_vec2.append(0)
        
                            if transitionType[index_rad][iDaughter][i_level][index_t] == "EL":       # record that an electron is missing on the L1 shell of the dughter nucleus
                                particle_vec2.append("Atom_L")
                                energy_vec2.append(0)
        
                            if transitionType[index_rad][iDaughter][i_level][index_t] == "EL1":       # record that an electron is missing on the L1 shell of the dughter nucleus
                                particle_vec2.append("Atom_L1")
                                energy_vec2.append(0)
        
                            if transitionType[index_rad][iDaughter][i_level][index_t] == "EL2":       # record that an electron is missing on the L2 shell of the dughter nucleus
                                particle_vec2.append("Atom_L2")
                                energy_vec2.append(0)
        
                            if transitionType[index_rad][iDaughter][i_level][index_t] == "EL3":       # record that an electron is missing on the L3 shell of the dughter nucleus
                                particle_vec2.append("Atom_L3")
                                energy_vec2.append(0)
        
                            if transitionType[index_rad][iDaughter][i_level][index_t] == "EM":        # record that an electron is missing on the M shell of the dughter nucleus
                                particle_vec2.append("Atom_M")
                                energy_vec2.append(0)
        
                            if transitionType[index_rad][iDaughter][i_level][index_t] == "EN":        # record that an electron is missing on the N shell of the dughter nucleus
                                particle_vec2.append("Atom_N")
                                energy_vec2.append(0)
                        e_sum2 += e_trans[index_rad][iDaughter][i_level][index_t]                      # Energy summary  
                    
                    
                    ## evenement normal    
                    else:
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
                            
                    
                    
                    levelOftheDaughter = next_level[index_rad][iDaughter][i_level][index_t]           # set the next level
                
                else:
                    i_level = levelNumber[index_rad][iDaughter].index([levelOftheDaughter])
                    print("warning:pas de données de transition:daughter,niveau,niveau d'énergie",DaughterVec[index_rad][iDaughter],levelOftheDaughter,levelEnergy[index_rad][iDaughter][i_level] )
                    levelOftheDaughter = 0   # set the next level
            
            if Display:
                print("\n\t NUCLEAR DECAY--Prompt")
                if "Atom" in particle_vec[0]:
                    print(f'\t\t capture of an electron from the {particle_vec[0][5:]} shell')
                else:
                    print(f'\t\t {particle_vec[0]} transition of energy = {energy_vec[0]}, keV')
                
                for i in range(1,len(particle_vec)):
                    if "Atom" in particle_vec[i]:
                        print(f'\t\t an electron of intern conversion from the {particle_vec[i][5:]} shell')
                    else:
                        print(f'\t\t emitted {particle_vec[i]} of energy = {energy_vec[i]}, keV')
                        
                if evenement != 1:
                    print("\n\t NUCLEAR DECAY--Delayed")
                    for i, p in enumerate(particle_vec2):
                        if p[:4] != "Atom":
                            print(f'\t\t\t {p} transition of energy = {energy_vec2[i]}, keV')
                        else:
                            print(f'\t\t\t an electron of intern conversion from the {p[5:]} shell')
    
            '''
            ==========================
            II. LA RELAXATION ATOMIQUE
            ==========================
            '''

            daughter_relax = DaughterVec[index_rad][iDaughter]
            particle_vec, energy_vec = relaxAtom(daughter_relax,particle_vec,energy_vec,Rad[index_rad],Display=Display,uncData=uncData)


            ## evenement normal

            if Display:
                print("\n\t ATOMIC RECOMBINATION--Prompt")
                for i, p in enumerate(particle_vec):
                    if p[:4] != "Atom":
                        if p=="beta" or p=="beta+":
                            print(f'\t\t {p} transition of energy = {energy_vec[i]}, keV')
                        else:
                            print(f"\t\t emitted {p} of energy = {round(energy_vec[i],3)} keV")
                    else:
                        print(f'\t\t an electron left the {p[5:]} shell')  
            
            ## evenement retardee             
            if evenement != 1:
                
                if Display:print("\n\t ATOMIC RECOMBINATION--Delay\n\t Summary of the atomic relaxation")
                
                particle_vec2, energy_vec2 = relaxAtom(daughter_relax,particle_vec2,energy_vec2,Rad[index_rad],Display=Display,uncData=uncData)
                

                if Display:
                    print("\n\t ATOMIC RECOMBINATION--Delay")
                    for i, p in enumerate(particle_vec2):
                        if p[:4] != "Atom":
                            if p=="beta" or p=="beta+":
                                print(f'\t\t {p} transition of energy = {energy_vec2[i]}, keV')
                            else:
                                print(f"\t\t emitted {p} of energy = {round(energy_vec2[i],3)} keV")
                        else:
                            print(f'\t\t an electron left the {p[5:]} shell')
    
                '''
                ==========================================================
                III.a SPECTRES D'EMISSION
                ==========================================================
                '''
            if ("beta" in particle_vec) or ("beta+" in particle_vec):
                if Display: print("\n\t EMISSION OF BETA PARTICLES")   
            for i, p in enumerate(particle_vec):
                if p == "beta":
                    e_b,p_b = tl.readBetaShape(rad_i,"beta-",level_before_trans)   # read the data of BetaShape
                    index_beta_energy = tl.sampling(p_b)                           # sampling energy of beta
                    particle_vec[i] = "electron"
                    energy_vec[i] = e_b[index_beta_energy]
                    if Display: print(f"\t\t emitted {p} of energy = {round(energy_vec[i],3)} keV")
                        
                if p == "beta+":
                    e_b,p_b = tl.readBetaShape(rad_i,"beta+",level_before_trans)
                    index_beta_energy = tl.sampling(p_b)
                    particle_vec[i] = "positron"
                    energy_vec[i] = e_b[index_beta_energy]
                    particle_vec.append("gamma")
                    particle_vec.append("gamma")
                    energy_vec.append(511)
                    energy_vec.append(511)
                    if Display: print(f"\t\t emitted {p} of energy = {round(energy_vec[i],3)} keV")
            #energy_vec_initial = energy_vec    
            energy_vec_initial = energy_vec.copy()
            '''
            ==========================================================
                III.b INTERACTION RAYONNEMENT/MATIERE
            ==========================================================
            '''

            for i, p in enumerate(particle_vec):
                if p == "electron":
                    energy_vec[i] = tl.energie_dep_beta2(energy_vec[i],v=V)
        
                if p == "beta+":
                    energy_vec[i] = tl.energie_dep_beta2(energy_vec[i],v=V)
        
                if p == "gamma" or p == "XKA" or p == "XKB" or p == "XL":
                    p0 = particle_vec[i]
                    Ei = energy_vec[i]
                    Ed = tl.energie_dep_gamma2(Ei,v=V)          # sampling energy free from photon
                    if Ei == Ed: # effet photoelectrique
                        energie_ele_emis,lacune_ph,element_ph = tl.interaction_scintillation(Ed)
                        particule_emise_ph,energie_par_emise_ph,posi_lacune_ph,par_emise_ph = tl.relaxation_atom_ph(lacune_ph,element_ph,v=V)
                        particle_vec = particle_vec + par_emise_ph
                        energy_vec_initial[i]=energie_ele_emis
                        energy_vec[i]=energie_ele_emis # energie du photoélectron primaire
                        energy_vec_initial = energy_vec_initial + energie_par_emise_ph
                        energy_vec = energy_vec + energie_par_emise_ph
                    else: # diffusion Compton
                        energy_vec[i]=Ed
                    particle_vec[i] = "electron"
                        
                if "Auger" in p:
                    particle_vec[i] = "electron"
                    energy_vec[i] = tl.energie_dep_beta2(energy_vec[i],v=V)
            
            if Display:
                print("\n\t INTERACTION--Prompt \n\t Summary of the energy deposited by charged particles")
                for i, p  in enumerate(particle_vec):
                    if p[:4] != "Atom" and energy_vec[i]!=0:
                        if p == "gamma" or p == "XKA" or p == "XKB" or p == "XL": (f"\t\t the {p0} gives {energy_vec[i]} keV to a recoil electron")
                        else: print(f"\t\t {p} of energy = {round(energy_vec[i],3)} keV")
                        
            
            if evenement!=1:
                energy_vec_initial2 = energy_vec2.copy()
                for i, p in enumerate(particle_vec2):
                    if p == "electron":
                        energy_vec2[i] = tl.energie_dep_beta2(energy_vec2[i],v=V)
        
                    if p == "gamma" or p == "XKA" or p == "XKB" or p == "XL":
                        p0 = particle_vec2[i]
                        Ei_2 = energy_vec2[i]
                        Ed_2 = tl.energie_dep_gamma2(Ei_2,v=V)          # sampling energy free from photon
                        if Ei_2 == Ed_2:
                            energie_ele_emis2,lacune_ph2,element_ph2 = tl.interaction_scintillation(Ed_2)
                            particule_emise_ph2,energie_par_emise_ph2,posi_lacune_ph2,par_emise_ph2 = tl.relaxation_atom_ph(lacune_ph2,element_ph2,v=V)
                            energy_vec2[i]=energie_ele_emis2
                            energy_vec_initial2[i]=energie_ele_emis2
                            energy_vec2 = energy_vec2 + energie_par_emise_ph2
                            energy_vec_initial2 = energy_vec_initial2 + energie_par_emise_ph2
                            particle_vec2 = particle_vec2 + par_emise_ph2
                        else: # diffusion Compton
                            energy_vec2[i]=Ed_2    
                        particle_vec2[i] = "electron"
                        if Display:
                            print(f"\t\t {p0} give energy {energy_vec2[i]} keV to electron")
                            
                        
                    if "Auger" in p:
                        particle_vec2[i] = "electron"
                        energy_vec2[i] = tl.energie_dep_beta2(energy_vec2[i],v=V)
                
                if Display:
                    print("\n\t INTERACTION--Delay \n\t Summary of the energy deposited by charged particles")
                    for i, p  in enumerate(particle_vec2):
                        if p[:4] != "Atom" and energy_vec2[i]!=0:
                            if p == "gamma" or p == "XKA" or p == "XKB" or p == "XL": (f"\t\t the {p0} gives {energy_vec2[i]} keV to a recoil electron")
                            else: print(f"\t\t {p} of energy = {round(energy_vec2[i],3)} keV")
                
                '''
                ====================
                IV. LA SCINTILLATION
                Calculation of the scintillation quenching with the Birks Model
                ====================
                '''
            if Display: print(f"\n\t SCINTILLATION--Prompt \n\t\t Birks constant = {kB} cm/keV\n\t Summary of the estimation of quenched energies")
            e_quenching=[]
            for i, p in enumerate(particle_vec):
                if p == "alpha":
                    energy_vec[i] = tl.Em_a(energy_vec[i],kB,nE_alpha)
                    e_quenching.append(energy_vec[i])
                elif p == "electron" or p == "positron":
                    energy_vec[i] = tl.Em_e(energy_vec_initial[i]*1e3,energy_vec[i]*1e3,kB*1e3,nE_electron)*1e-3
                    e_quenching.append(energy_vec[i])
                else:
                    e_quenching.append(0)
            if Display: print("\t\t Birks constant = ", kB, ' cm/keV')
            if Display:
                for i, p in enumerate(particle_vec):
                    #print(e_quenching[i])
                    if p[:4] != "Atom": print(f"\t\t quenched energy of {p} = ", np.round(e_quenching[i],3), "keV")
            
            if evenement!=1:
                if Display: print(f"\n\t SCINTILLATION--Delayed \n\t\t Birks constant = {kB} cm/keV\n\t Summary of the estimation of quenched energies")
                e_quenching2=[]
                for i, p in enumerate(particle_vec2):
                    if p == "alpha":
                        energy_vec2[i] = tl.Em_a(energy_vec2[i],kB,nE_alpha)
                        e_quenching2.append(energy_vec2[i])
                    elif p == "electron" or p == "positron":
                        energy_vec2[i] = tl.Em_e(energy_vec_initial2[i]*1e3,energy_vec2[i]*1e3,kB*1e3,nE_electron)*1e-3
                        e_quenching2.append(energy_vec2[i])
                    else:
                        e_quenching2.append(0) 
                if Display:
                    for i, p in enumerate(particle_vec2):
                        if p[:4] != "Atom": print(f"\t\t quenched energy of {p} = ", round(e_quenching2[i],3), "keV")       
                
            '''
            ====================
            V. LE MESURE TDCR
            ====================
            '''            
            if mode2=="sym":
                if evenement !=1:
                    p_nosingle = np.exp(-L*np.sum(np.asarray(e_quenching))/3) # probability to have 0 electrons in a PMT
                    p_single = 1-p_nosingle                                    # probability to have at least 1 electrons in a PMT
                    p_nosingle2 = np.exp(-L*np.sum(np.asarray(e_quenching2))/3) # probability to have 0 electrons in a PMT
                    p_single2 = 1-p_nosingle2 
                    efficiency_S.append(p_single+p_single2)
                    efficiency_T.append(p_single**3+p_single2**3)
                    efficiency_D.append(3*(p_single)**2-2*p_single**3+(3*(p_single2)**2-2*p_single2**3))
                    if Display: print(f"\n\t COUNTING--Sym \n\t\t Free parameter = {L} keV-1 \n\t Summary of TDCR measurement (prompt)")
                    if Display: print("\t\t Free parameter = ", L, "keV-1")
                    if Display: print("\t\t Efficiency of single events = ", round(p_single,5))
                    if Display: print("\t\t Efficiency of double events = ", round(3*(p_single)**2-2*p_single**3,5))
                    if Display: print("\t\t Efficiency of triple events = ", round(p_single**3,5))
                    if Display: print("\t Summary of TDCR measurement (delayed)")
                    if Display: print("\t\t Efficiency of single events = ", round(p_single2,5))
                    if Display: print("\t\t Efficiency of double events = ", round(3*(p_single2)**2-2*p_single2**3,5))
                    if Display: print("\t\t Efficiency of triple events = ", round(p_single2**3,5))
                    if Display: print("\t Summary of TDCR measurement (prompt + delayed)")
                    if Display: print("\t\t Efficiency of single events = ", round(p_single+p_single2,5))
                    if Display: print("\t\t Efficiency of double events = ", round(3*(p_single)**2-2*p_single**3+(3*(p_single2)**2-2*p_single2**3),5))
                    if Display: print("\t\t Efficiency of triple events = ", round(p_single**3+p_single2**3,5))
                else:
                    p_nosingle = np.exp(-L*np.sum(np.asarray(e_quenching))/3) # probability to have 0 electrons in a PMT
                    p_single = 1-p_nosingle                                    # probability to have at least 1 electrons in a PMT
                    efficiency_S.append(p_single)
                    efficiency_T.append(p_single**3)
                    efficiency_D.append(3*(p_single)**2-2*efficiency_T[-1])
                    if Display: print(f"\n\t COUNTING--Sym \n\t\t Free parameter = {L} keV-1 \n\t Summary of TDCR measurement (prompt)")
                    if Display: print("\t\t Efficiency of single events = ", round(efficiency_S[-1],5))
                    if Display: print("\t\t Efficiency of double events = ", round(efficiency_D[-1],5))
                    if Display: print("\t\t Efficiency of triple events = ", round(efficiency_T[-1],5))                    
                                    
            elif mode2=="asym":
                if evenement !=1:
                    pA_nosingle = np.exp(-L[0]*np.sum(np.asarray(e_quenching))/3) # probability to have 0 electrons in a PMT
                    pA_single = 1-pA_nosingle                                    # probability to have at least 1 electrons in a PMT
                    pB_nosingle = np.exp(-L[1]*np.sum(np.asarray(e_quenching))/3) # probability to have 0 electrons in a PMT
                    pB_single = 1-pB_nosingle                                    # probability to have at least 1 electrons in a PMT
                    pC_nosingle = np.exp(-L[2]*np.sum(np.asarray(e_quenching))/3) # probability to have 0 electrons in a PMT
                    pC_single = 1-pC_nosingle                                    # probability to have at least 1 electrons in a PMT
                    
                    pA_nosingle2 = np.exp(-L[0]*np.sum(np.asarray(e_quenching2))/3) # probability to have 0 electrons in a PMT
                    pA_single2 = 1-pA_nosingle2                                    # probability to have at least 1 electrons in a PMT
                    pB_nosingle2 = np.exp(-L[1]*np.sum(np.asarray(e_quenching2))/3) # probability to have 0 electrons in a PMT
                    pB_single2 = 1-pB_nosingle2                                    # probability to have at least 1 electrons in a PMT
                    pC_nosingle2 = np.exp(-L[2]*np.sum(np.asarray(e_quenching2))/3) # probability to have 0 electrons in a PMT
                    pC_single2 = 1-pC_nosingle2                                    # probability to have at least 1 electrons in a PMT
                    
                    efficiency_AB.append(pA_single*pB_single+pA_single2*pB_single2)
                    efficiency_BC.append(pB_single*pC_single+pB_single2*pC_single2)
                    efficiency_AC.append(pA_single*pC_single+pA_single2*pC_single2)
                    efficiency_T.append(pA_single*pB_single*pC_single+pA_single2*pB_single2*pC_single2)
                    efficiency_D.append(pA_single*pB_single+pB_single*pC_single+pA_single*pC_single-2*pA_single*pB_single*pC_single+(pA_single2*pB_single2+pB_single2*pC_single2+pA_single2*pC_single2-2*pA_single2*pB_single2*pC_single2))
                    efficiency_S.append(pA_single+pB_single+pC_single-pA_single*pB_single+pB_single*pC_single+pA_single*pC_single-2*pA_single*pB_single*pC_single-pA_single*pB_single*pC_single+(pA_single2+pB_single2+pC_single2-pA_single2*pB_single2+pB_single2*pC_single2+pA_single2*pC_single2-2*pA_single2*pB_single2*pC_single2-pA_single2*pB_single2*pC_single2))
                    
                    if Display: print(f"\n\t COUNTING--Asym \n\t\t Free parameters (A,B,C) = {L[0]},{L[1]},{L[2]} keV-1 \n\t Summary of TDCR measurement (prompt)")
                    #if Display: print("\t Summary of TDCR measurement (prompt)")
                    if Display: print("\t\t Efficiency of single events: ", round(pA_single+pB_single+pC_single-pA_single*pB_single+pB_single*pC_single+pA_single*pC_single-2*pA_single*pB_single*pC_single-pA_single*pB_single*pC_single,5))
                    if Display: print("\t\t Efficiency of double events: ", round(pA_single*pB_single+pB_single*pC_single+pA_single*pC_single-2*pA_single*pB_single*pC_single,5))
                    if Display: print("\t\t Efficiency of triple events: ", round(pA_single*pB_single*pC_single,5))
                    if Display: print("\t Summary of TDCR measurement (delayed)")
                    if Display: print("\t\t Efficiency of single events: ", round(pA_single2+pB_single2+pC_single2-pA_single2*pB_single2+pB_single2*pC_single2+pA_single2*pC_single2-2*pA_single2*pB_single2*pC_single2-pA_single2*pB_single2*pC_single2,5))
                    if Display: print("\t\t Efficiency of double events: ", round(pA_single2*pB_single2+pB_single2*pC_single2+pA_single2*pC_single2-2*pA_single2*pB_single2*pC_single2,5))
                    if Display: print("\t\t Efficiency of triple events: ", round(pA_single2*pB_single2*pC_single2,5))
                else:
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
                    if Display: print(f"\n\t COUNTING--Asym \n\t\t Free parameters (A,B,C) = {L[0]},{L[1]},{L[2]} keV-1 \n\t Summary of TDCR measurement (prompt)")
                    if Display: print("\t\t Free parameter PMT A: ", L[0], "keV-1")
                    if Display: print("\t\t Free parameter PMT B: ", L[1], "keV-1")
                    if Display: print("\t\t Free parameter PMT C: ", L[2], "keV-1")
                    if Display: print("\t\t Efficiency of single events: ", round(efficiency_S[-1],5))
                    if Display: print("\t\t Efficiency of double events: ", round(efficiency_D[-1],5))
                    if Display: print("\t\t Efficiency of triple events: ", round(efficiency_T[-1],5))                    


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
        
        if mode2=="sym":
            res=(TDCR_calcul-TD)**2
        elif mode2=="asym":
            res=(TAB-TABmodel)**2+(TBC-TBCmodel)**2+(TAC-TACmodel)**2
        
        if mode == "res":
            return res
        if mode == "eff":
            if N<200:
                print("Warning. too low number of MC trials - inaccurate estimation")
                return mean_efficiency_S, 1, mean_efficiency_D, 1, mean_efficiency_T, 1
            else:
                return mean_efficiency_S, std_efficiency_S, mean_efficiency_D, std_efficiency_D, mean_efficiency_T, std_efficiency_T
        if mode =="dis":
            return efficiency_S, efficiency_D, efficiency_T    

# L = 1
# TD = 0.977667386529166
# TAB = 0.992232838598821
# TBC = 0.992343419459002
# TAC = 0.99275350064608
# Rad="Fe-55"
# pmf_1="1"
# N = 1000
# kB =1.0e-5
# V = 10
# mode = "eff"
# mode2 = "sym"

# out = TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V, mode, mode2, Display=True, barp=False,uncData=False)
# print("TDCR", out[4]/out[2])
# print("Eff D", out[2])