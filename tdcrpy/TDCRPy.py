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
from importlib.resources import files
import configparser
import numpy as np
from tqdm import tqdm
import tempfile
import os
import scipy.optimize as opt

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

def TDCRPy(L, Rad, pmf_1, N, kB, V, mode="eff", Display=False, barp=False, Smodel=True, record = False, readRecHist = False, uncData = False, fullMC = False):
    """
    This is the main function of the TDCRPy package.
    The computation is made for a given solution containing a radionuclide (or a mixture of radionuclides), a given volume of scintillator V and a given Birks constant kB. 
       
    The parmeter N sets the number of Monte-Carlo trails used for the estimation. Each MC trial corresponds to a simulated radiactive decay.
    TDCRPY() used a set of fonctions from the tdcrpy.TDCR_model_lib module.
    
    Advanced settings can be configured in the config.toml file.
    
    Parameters
    ----------
    L : Float or tuple
        If L is float, then L is the global free parameter. If L is tuple, then L is a triplet of free parameters. unit keV-1
    Rad : string
        List of radionuclides (eg. "H-3, Co-60").
    pmf_1 : string
        list of probability of each radionuclide (eg. "0.8, 0.2").
    N : integer
        Number of Monte-Carlo trials. recommanded N>10000 (see JCGM 101). Not applied in the case of the analytical model.
    kB : float
        Birks constant in cm/keV.
    V : float
        volume of the scintillator in ml.
    mode : string
        "eff" to return efficiencies, "dis" to return list of decay events.
    Display : Boolean, optional
        "True" to display details on the decay sampling. The default is False.
    barp : Boolean, optional
        "True" to display the calculation progress. The default is True.
    record : Boolean, optional
        "True" to record details of decay events. The default is False.
    Smodel : Boolean, optional
        "True" to run the stochastic TDCR model. False to run the analytical caclulation (available only for pure beta emitters).
    
    Returns
    -------
    mean_efficiency_S : float
        Estimation of the efficiency of single counting events.
    std_efficiency_S : float
        Standard uncertainty from calculation associated with the estimation of the efficiency of single counting events.
    mean_efficiency_D : float
        Estimation of the efficiency of logic sum of double coincidences.
    std_efficiency_D : float
        Standard uncertainty from calculation associated with the estimation of the efficiency of logic sum of double coincidences.
    mean_efficiency_T : float
        Estimation of the efficiency of triple coincidences.
    std_efficiency_T : float
        Standard uncertainty from calculation associated with the estimation of the efficiency of triple coincidences.
    mean_efficiency_AB : float
        detection efficiency of coincidences between channels A and B.
    std_efficiency_AB : float
        standard uncertainty of detection efficiency of coincidences between channels A and B.
    mean_efficiency_BC : float
        detection efficiency of coincidences between channels B and C.
    std_efficiency_BC : float
        standard uncertainty of Ddetection efficiency of coincidences between channels B and C.
    mean_efficiency_AC : float
        detection efficiency of coincidences between channels A and C.
    std_efficiency_AC : float
        standard uncertainty of detection efficiency of coincidences between channels A and C.
    mean_efficiency_D2 : float
        detection efficiency of coincidences in a C/N system.
    std_efficiency_D2 : float
        standard uncertainty of detection efficiency of coincidences in a C/N system.
    """
    
    if isinstance(L, (tuple, list)):
        symm = False
    else:
        symm = True
    
    
    if record:
        temp_dir = tempfile.gettempdir()
        recfile1 = os.path.join(temp_dir, "Temp_E0.txt")
        header_content1 = """# TDCRPy output: inital energies from nuclear decays
# Column 1: KPAR (1=electron, 2=photon, 3=positron, 4=alpha)
# Column 2: Energy in eV  (> ECNUC)
# Column 3: Decay number mod 100, two digits
# Column 4: Cascade number mod 10, one digit
# Column 5: Particle age in seconds (since decay started)
# 2            3  4 5
"""
        with open(recfile1, 'w') as file: file.write(header_content1)

        recfile2 = os.path.join(temp_dir, "Temp_E1.txt")
        header_content2 = """# TDCRPy output: deposited energies from nuclear decays
# Column 1: KPAR (1=electron, 2=photon, 3=positron, 4=alpha)
# Column 2: Energy in eV  (> ECNUC)
# Column 3: Decay number mod 100, two digits
# Column 4: Cascade number mod 10, one digit
# Column 5: Particle age in seconds (since decay started)
# 2            3  4 5
"""
        with open(recfile2, "w") as file: file.write(header_content2)

        recfile3 = os.path.join(temp_dir, "Temp_E2.txt")
        header_content3 = """# TDCRPy output: deposited quenched energies from nuclear decays
# Column 1: KPAR (1=electron, 2=photon, 3=positron, 4=alpha)
# Column 2: Energy in eV  (> ECNUC)
# Column 3: Decay number mod 100, two digits
# Column 4: Cascade number mod 10, one digit
# Column 5: Particle age in seconds (since decay started)
# 2            3  4 5
"""
        with open(recfile3, "w") as file: file.write(header_content3)
        
        header_content4 = """# TDCRPy output: detection probabilities
# Column 1: Decay number mod 100, two digits
# Column 2: detection probability of single events
# Column 3: detection probability if double coincidences
# Column 4: detection probability if triple coincidences
# 2 3 4
"""
        recfile4 = os.path.join(temp_dir, "Temp_E3.txt")
        with open(recfile4, "w") as file: file.write(header_content4)
            
    if barp: tl.display_header()
    config = configparser.ConfigParser()
    with importlib.resources.as_file(files('tdcrpy').joinpath('config.toml')) as data_path:
        file_conf = data_path       
    config.read(file_conf)
    tau=config["Inputs"].getfloat("tau")
    extDT=config["Inputs"].getfloat("extDT")
    measTime=config["Inputs"].getfloat("measTime")
    micCorr=config["Inputs"].getboolean("micCorr")
    radListPureBeta = ["H-3", "C-14", "S-35", "Ca-45", "Ni-63", "Sr-89", "Sr-90", "Tc-99", "Pm-147", "Pu-241"]
    if (Rad in radListPureBeta) and (not Smodel):
        nElist = [7000, 1000, 1000, 500, 2000, 500, 200, 500, 1000, 7000] # discretization
        inE = radListPureBeta.index(Rad)
        nE = nElist[inE]
        # print(f"Analytical model used for {Rad}")
        out=tl.modelAnalytical(L,1,1,1,1,Rad,kB,V,mode,symm,nE)
        if mode == "eff":
            return out[0], 0, out[1], 0, out[2], 0
    elif (not Smodel) and (not Rad in radListPureBeta):
        # print("cannot be processed by the analytical model.")
        # print(f"Analytical model used for {Rad}")
        out=tl.modelAnalytical(L,1,1,1,1,Rad,kB,V,mode,symm,1000)
        if mode == "eff":
            return out[0], 0, out[1], 0, out[2], 0
    
    elif readRecHist:
        efficiency_A2 = []
        efficiency_B2 = []
        efficiency_S = []
        efficiency_D = []
        efficiency_D2 = []
        efficiency_T = []
        efficiency_AB = []
        efficiency_BC = []
        efficiency_AC = []
        
        temp_dir = tempfile.gettempdir()
        recfile3 = os.path.join(temp_dir, "Temp_E2.txt")

        with open(recfile3, "r") as file:
   
            decaym = -1
            e_quenching = []; e_quenching2 = []; evenement=1; t1=0
            for line in file:
                if line[0] != "#":
                    line = line.split(' ')
                    line = [element for element in line if element != ""]
                    decay = int(line[2])
                    
                    if decay != decaym:
                        if decay>0:
                            # print(decay-1,e_quenching,e_quenching2, evenement)
                            if fullMC:
                                efficiency0_S, efficiency0_D, efficiency0_T, efficiency0_AB, efficiency0_BC, efficiency0_AC, efficiency0_D2 = tl.detectProbabilitiesMC(L, e_quenching, e_quenching2, t1, evenement, extDT, measTime)
                            else:
                                efficiency0_S, efficiency0_D, efficiency0_T, efficiency0_AB, efficiency0_BC, efficiency0_AC, efficiency0_D2 = tl.detectProbabilities(L, e_quenching, e_quenching2, t1, evenement, extDT, measTime)
                            efficiency_S.append(efficiency0_S)
                            efficiency_T.append(efficiency0_T)
                            efficiency_D.append(efficiency0_D)
                            efficiency_AB.append(efficiency0_AB)
                            efficiency_BC.append(efficiency0_BC)
                            efficiency_AC.append(efficiency0_AC)
                            efficiency_D2.append(efficiency0_D2)
                            
                            # print(efficiency0_D, "\n")
                        
                        energy = float(line[1])*1e-3
                        t1 = float(line[4])
                        decaym = decay
                        # print(decay, energy, t1, extDT)
                        e_quenching = []; e_quenching2 = []
                        evenement=1
                        e_quenching.append(energy)
                    else:
                        energy = float(line[1])*1e-3
                        t1 = float(line[4])
                        # print(decay, energy, t1, extDT)
                        if t1 > tau*1e-9:
                            evenement = evenement + 1
                            e_quenching2.append(energy)
                        else:
                            e_quenching.append(energy)
            
            if fullMC:
                efficiency0_S, efficiency0_D, efficiency0_T, efficiency0_AB, efficiency0_BC, efficiency0_AC, efficiency0_D2 = tl.detectProbabilitiesMC(L, e_quenching, e_quenching2, t1, evenement, extDT, measTime)
            else:
                efficiency0_S, efficiency0_D, efficiency0_T, efficiency0_AB, efficiency0_BC, efficiency0_AC, efficiency0_D2 = tl.detectProbabilities(L, e_quenching, e_quenching2, t1, evenement, extDT, measTime)
            
            efficiency_S.append(efficiency0_S)
            efficiency_T.append(efficiency0_T)
            efficiency_D.append(efficiency0_D)
            efficiency_AB.append(efficiency0_AB)
            efficiency_BC.append(efficiency0_BC)
            efficiency_AC.append(efficiency0_AC)
            efficiency_D2.append(efficiency0_D2)
                            
            # print(efficiency_D)
            outEff = tl.efficienciesEstimates(efficiency_S, efficiency_D, efficiency_T, efficiency_AB, efficiency_BC, efficiency_AC, efficiency_D2, N)
            if mode == "eff": return outEff
            if mode == "dis": return efficiency_S, efficiency_D, efficiency_T, efficiency_D2
                

            # temp1, temp2, temp3, temp4 = tl.read_temp_files()
            # ee_vec = tl.energyVectors3(temp3)
            # # to be continued
    
    else:
        efficiency_A2 = []
        efficiency_B2 = []
        efficiency_S = []
        efficiency_D = []
        efficiency_D2 = []
        efficiency_T = []
        efficiency_AB = []
        efficiency_BC = []
        efficiency_AC = []
        
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
                
        if barp and not Display: NN = tqdm(range(N), desc="Processing", unit=" decays")
        else: NN = range(N)
        for idec in NN: # Main Loop - Monte Carlo trials
            particle_vec=[]
            energy_vec=[]
            '''
            ===============================
            0. SAMPLING OF THE RADIONUCLIDE
            ===============================
            '''
            index_rad = tl.sampling(pmf_1)
            rad_i = Rad[index_rad]
            if Display: print("\n\n Trial ",str(idec+1),"- Sampled radionuclide: ", rad_i)
                  
            
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
                
                t1 = np.random.exponential(trans_halfLife[index_rad][iDaughter][i_level][0]/np.log(2), size=1)[0]

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
                for idisp, p in enumerate(particle_vec):
                    if p[:4] != "Atom":
                        if p=="beta" or p=="beta+":
                            print(f'\t\t {p} transition of energy = {energy_vec[idisp]}, keV')
                        else:
                            print(f"\t\t emitted {p} of energy = {round(energy_vec[idisp],3)} keV")
                    else:
                        print(f'\t\t an electron left the {p[5:]} shell')  
            
            ## evenement retardee             
            if evenement != 1:
                
                if Display:print("\n\t ATOMIC RECOMBINATION--Delay\n\t Summary of the atomic relaxation")
                
                particle_vec2, energy_vec2 = relaxAtom(daughter_relax,particle_vec2,energy_vec2,Rad[index_rad],Display=Display,uncData=uncData)
                

                if Display:
                    print("\n\t ATOMIC RECOMBINATION--Delay")
                    for idisp, p in enumerate(particle_vec2):
                        if p[:4] != "Atom":
                            if p=="beta" or p=="beta+":
                                print(f'\t\t {p} transition of energy = {energy_vec2[idisp]}, keV')
                            else:
                                print(f"\t\t emitted {p} of energy = {round(energy_vec2[idisp],3)} keV")
                        else:
                            print(f'\t\t an electron left the {p[5:]} shell')
    
                '''
                ==========================================================
                III.a SPECTRES D'EMISSION
                ==========================================================
                '''
            if ("beta" in particle_vec) or ("beta+" in particle_vec):
                if Display: print("\n\t EMISSION OF BETA PARTICLES")   
            for ipart, p in enumerate(particle_vec):
                if p == "beta":
                    e_b,p_b = tl.readBetaShape(rad_i,"beta-",level_before_trans)   # read the data of BetaShape
                    index_beta_energy = tl.sampling(p_b)                           # sampling in PDF
                    particle_vec[ipart] = "electron"
                    energy_vec[ipart] = e_b[index_beta_energy]
                    if Display: print(f"\t\t emitted {p} of energy = {round(energy_vec[ipart],3)} keV")
                        
                if p == "beta+":
                    e_b,p_b = tl.readBetaShape(rad_i,"beta+",level_before_trans)
                    index_beta_energy = tl.sampling(p_b)                          # sampling in PDF
                    particle_vec[ipart] = "positron"
                    energy_vec[ipart] = e_b[index_beta_energy]
                    particle_vec.append("gamma")
                    particle_vec.append("gamma")
                    energy_vec.append(511)
                    energy_vec.append(511)
                    if Display: print(f"\t\t emitted {p} of energy = {round(energy_vec[ipart],3)} keV")
            #energy_vec_initial = energy_vec    
            energy_vec_initial = energy_vec.copy()
            
            
            
            """
            ==========================================================
            III.a.1 record initial energies in file
            ==========================================================
            """
            if record:
                with open(recfile1, "a") as file:
                    # t1w=0
                    for irec, p in enumerate(particle_vec):
                        writeOn = True
                        
                        if p == "electron" or ("Auger" in p): col1 = "1"
                        elif p == "gamma" or ("X" in p): col1 = "2"
                        elif p == "positron": col1 = "3"
                        elif p == "alpha": col1 = "4"
                        else: writeOn = False

                        # if ('t1' in globals()) or ('t1' in locals()): t1w+=t1
                        
                        if writeOn:
                            if idec<10:
                                file.write(f"{col1} {energy_vec[irec]*1e3:.6E}  {idec} 1 0\n")
                            else:
                                file.write(f"{col1} {energy_vec[irec]*1e3:.6E} {idec} 1 0\n")
                    
                    if evenement != 1:
                        for irec, p in enumerate(particle_vec2):
                            writeOn = True
                            
                            if p == "electron" or ("Auger" in p): col1 = "1"
                            elif p == "gamma" or ("X" in p): col1 = "2"
                            elif p == "positron": col1 = "3"
                            elif p == "alpha": col1 = "4"
                            else: writeOn = False
    
                            # if ('t1' in globals()) or ('t1' in locals()): t1w+=t1
                            
                            if writeOn:
                                if idec<10:
                                    file.write(f"{col1} {energy_vec2[irec]*1e3:.6E}  {idec} 1 {t1:.6E}\n")
                                else:
                                    file.write(f"{col1} {energy_vec2[irec]*1e3:.6E} {idec} 1 {t1:.6E}\n")
                           
            
            
            '''
            ==========================================================
                III.b INTERACTION RAYONNEMENT/MATIERE
            ==========================================================
            '''

            for ipart, p in enumerate(particle_vec):
                if p == "electron":
                    energy_vec[ipart] = tl.energie_dep_beta2(energy_vec[ipart],v=V)
        
                if p == "beta+":
                    energy_vec[ipart] = tl.energie_dep_beta2(energy_vec[ipart],v=V)
        
                if p == "gamma" or p == "XKA" or p == "XKB" or p == "XL":
                    p0 = particle_vec[ipart]
                    Ei = energy_vec[ipart]
                    Ed = tl.energie_dep_gamma2(Ei,v=V)          # sampling energy free from photon
                    if Ei == Ed: # effet photoelectrique
                        energie_ele_emis,lacune_ph,element_ph = tl.interaction_scintillation(Ed)
                        particule_emise_ph,energie_par_emise_ph,posi_lacune_ph,par_emise_ph = tl.relaxation_atom_ph(lacune_ph,element_ph,v=V)
                        particle_vec = particle_vec + par_emise_ph
                        energy_vec_initial[ipart]=energie_ele_emis
                        energy_vec[ipart]=energie_ele_emis # energie du photoélectron primaire
                        energy_vec_initial = energy_vec_initial + energie_par_emise_ph
                        energy_vec = energy_vec + energie_par_emise_ph
                    elif Ed == Ei - 1022:           # creation de paire
                        #particle_vec[i] = 'electron'
                        E_e = (Ei-1022)/2
                        energy_vec[ipart] =  E_e#tl.energie_dep_beta2(E_e,v=V)
                        particle_vec.append("positron")
                        energy_vec.append(E_e)#tl.energie_dep_beta2(E_e,v=V))
                    else: # diffusion Compton
                        energy_vec[ipart]=Ed
                    particle_vec[ipart] = "electron"
                        
                if "Auger" in p:
                    particle_vec[ipart] = "electron"
                    energy_vec[ipart] = tl.energie_dep_beta2(energy_vec[ipart],v=V)
            
            if Display:
                print("\n\t INTERACTION--Prompt \n\t Summary of the energy deposited by charged particles")
                for ipart, p  in enumerate(particle_vec):
                    if p[:4] != "Atom" and energy_vec[ipart]!=0:
                        if p == "gamma" or p == "XKA" or p == "XKB" or p == "XL": (f"\t\t the {p0} gives {energy_vec[ipart]} keV to a recoil electron")
                        else: print(f"\t\t {p} of energy = {round(energy_vec[ipart],3)} keV")
                        
            
            if evenement!=1:
                energy_vec_initial2 = energy_vec2.copy()
                for ipart, p in enumerate(particle_vec2):
                    if p == "electron":
                        energy_vec2[ipart] = tl.energie_dep_beta2(energy_vec2[ipart],v=V)
        
                    if p == "gamma" or p == "XKA" or p == "XKB" or p == "XL":
                        p0 = particle_vec2[ipart]
                        Ei_2 = energy_vec2[ipart]
                        Ed_2 = tl.energie_dep_gamma2(Ei_2,v=V)          # sampling energy free from photon
                        if Ei_2 == Ed_2:        # effet photon-electrique
                            energie_ele_emis2,lacune_ph2,element_ph2 = tl.interaction_scintillation(Ed_2)
                            particule_emise_ph2,energie_par_emise_ph2,posi_lacune_ph2,par_emise_ph2 = tl.relaxation_atom_ph(lacune_ph2,element_ph2,v=V)
                            energy_vec2[ipart]=energie_ele_emis2
                            energy_vec_initial2[ipart]=energie_ele_emis2
                            energy_vec2 = energy_vec2 + energie_par_emise_ph2
                            energy_vec_initial2 = energy_vec_initial2 + energie_par_emise_ph2
                            particle_vec2 = particle_vec2 + par_emise_ph2
                        elif Ed_2 == Ei_2 - 1022:   # creation de paire
                            #particle_vec2[ipart] = 'electron'
                            E_e = (Ei_2-1022)/2
                            energy_vec2[ipart] =  E_e#tl.energie_dep_beta2(E_e,v=V)
                            particle_vec2.append("positron")
                            energy_vec2.append(E_e)#tl.energie_dep_beta2(E_e,v=V))
                        else: # diffusion Compton
                            energy_vec2[ipart]=Ed_2    
                        particle_vec2[ipart] = "electron"
                        if Display:
                            print(f"\t\t {p0} give energy {energy_vec2[ipart]} keV to electron")
                            
                        
                    if "Auger" in p:
                        particle_vec2[ipart] = "electron"
                        energy_vec2[ipart] = tl.energie_dep_beta2(energy_vec2[ipart],v=V)
                
                if Display:
                    print("\n\t INTERACTION--Delay \n\t Summary of the energy deposited by charged particles")
                    for ipart, p  in enumerate(particle_vec2):
                        if p[:4] != "Atom" and energy_vec2[ipart]!=0:
                            if p == "gamma" or p == "XKA" or p == "XKB" or p == "XL": (f"\t\t the {p0} gives {energy_vec2[ipart]} keV to a recoil electron")
                            else: print(f"\t\t {p} of energy = {round(energy_vec2[ipart],3)} keV")
                
            """
            ==========================================================
            III.b.1 record deposited energies in file
            ==========================================================
            """
            if record:
                with open(recfile2, "a") as file:
                    # t1w=0
                    for irec, p in enumerate(particle_vec):
                        writeOn = True
                        
                        if p == "electron" or ("Auger" in p): col1 = "1"
                        elif p == "gamma" or ("X" in p): col1 = "2"
                        elif p == "positron": col1 = "3"
                        elif p == "alpha": col1 = "4"
                        else: writeOn = False

                        # if ('t1' in globals()) or ('t1' in locals()): t1w+=t1
                        
                        if writeOn:
                            if idec<10:
                                file.write(f"{col1} {energy_vec[irec]*1e3:.6E}  {idec} 1 0\n")
                            else:
                                file.write(f"{col1} {energy_vec[irec]*1e3:.6E} {idec} 1 0\n")
                    
                    if evenement != 1:
                        for irec, p in enumerate(particle_vec2):
                            writeOn = True
                            
                            if p == "electron" or ("Auger" in p): col1 = "1"
                            elif p == "gamma" or ("X" in p): col1 = "2"
                            elif p == "positron": col1 = "3"
                            elif p == "alpha": col1 = "4"
                            else: writeOn = False
    
                            # if ('t1' in globals()) or ('t1' in locals()): t1w+=t1
                            
                            if writeOn:
                                if idec<10:
                                    file.write(f"{col1} {energy_vec2[irec]*1e3:.6E}  {idec} 1 {t1:.6E}\n")
                                else:
                                    file.write(f"{col1} {energy_vec2[irec]*1e3:.6E} {idec} 1 {t1:.6E}\n")
                
                
            '''
            ====================
            IV. LA SCINTILLATION
            Calculation of the scintillation quenching with the Birks Model
            ====================
            '''
            if Display: print(f"\n\t SCINTILLATION--Prompt \n\t\t Birks constant = {kB} cm/keV\n\t Summary of the estimation of quenched energies")
            e_quenching=[]
            for ipart, p in enumerate(particle_vec):
                if p == "alpha":
                    energy_vec[ipart] = tl.Em_a(energy_vec[ipart],kB,nE_alpha)
                    e_quenching.append(energy_vec[ipart])
                elif p == "electron" or p == "positron":
                    energy_vec[ipart] = tl.Em_e(energy_vec_initial[ipart]*1e3,energy_vec[ipart]*1e3,kB*1e3,nE_electron)*1e-3
                    if micCorr: energy_vec[ipart] = energy_vec[ipart]*tl.micelleLoss(energy_vec_initial[ipart])
                    e_quenching.append(energy_vec[ipart])
                else:
                    e_quenching.append(0)
            if Display: print("\t\t Birks constant = ", kB, ' cm/keV')
            if Display:
                for ipart, p in enumerate(particle_vec):
                    #print(e_quenching[i])
                    if p[:4] != "Atom": print(f"\t\t quenched energy of {p} = ", np.round(e_quenching[ipart],3), "keV")
            
            if evenement!=1:
                if Display: print(f"\n\t SCINTILLATION--Delayed \n\t\t Birks constant = {kB} cm/keV\n\t Summary of the estimation of quenched energies")
                e_quenching2=[]
                for ipart, p in enumerate(particle_vec2):
                    if p == "alpha":
                        energy_vec2[ipart] = tl.Em_a(energy_vec2[ipart],kB,nE_alpha)
                        e_quenching2.append(energy_vec2[ipart])
                    elif p == "electron" or p == "positron":
                        energy_vec2[ipart] = tl.Em_e(energy_vec_initial2[ipart]*1e3,energy_vec2[ipart]*1e3,kB*1e3,nE_electron)*1e-3
                        if micCorr: energy_vec2[ipart] = energy_vec2[ipart]*tl.micelleLoss(energy_vec_initial2[ipart])
                        e_quenching2.append(energy_vec2[ipart])
                    else:
                        e_quenching2.append(0) 
                if Display:
                    for ipart, p in enumerate(particle_vec2):
                        if p[:4] != "Atom": print(f"\t\t quenched energy of {p} = ", round(e_quenching2[ipart],3), "keV")
            
            
            """
            ==========================================================
            IV.1 record deposited quenched energies in file
            ==========================================================
            """
            if record:
                with open(recfile3, "a") as file:
                    # t1w=0
                    for irec, p in enumerate(particle_vec):
                        writeOn = True
                        
                        if p == "electron" or ("Auger" in p): col1 = "1"
                        elif p == "gamma" or ("X" in p): col1 = "2"
                        elif p == "positron": col1 = "3"
                        elif p == "alpha": col1 = "4"
                        else: writeOn = False

                        # if ('t1' in globals()) or ('t1' in locals()): t1w+=t1
                        
                        if writeOn:
                            if idec<10:
                                file.write(f"{col1} {e_quenching[irec]*1e3:.6E}  {idec} 1 0\n")
                            else:
                                file.write(f"{col1} {e_quenching[irec]*1e3:.6E} {idec} 1 0\n")
                    
                    if evenement != 1:
                        for irec, p in enumerate(particle_vec2):
                            writeOn = True
                            
                            if p == "electron" or ("Auger" in p): col1 = "1"
                            elif p == "gamma" or ("X" in p): col1 = "2"
                            elif p == "positron": col1 = "3"
                            elif p == "alpha": col1 = "4"
                            else: writeOn = False
    
                            # if ('t1' in globals()) or ('t1' in locals()): t1w+=t1
                            
                            if writeOn:
                                if idec<10:
                                    file.write(f"{col1} {e_quenching2[irec]*1e3:.6E}  {idec} 1 {t1:.6E}\n")
                                else:
                                    file.write(f"{col1} {e_quenching2[irec]*1e3:.6E} {idec} 1 {t1:.6E}\n")
            
            
                
            '''
            ====================
            V. LE MESURE TDCR
            ====================
            '''
            if evenement == 1: e_quenching2 = 0; t1=0
            if fullMC:
                efficiency0_S, efficiency0_D, efficiency0_T, efficiency0_AB, efficiency0_BC, efficiency0_AC, efficiency0_D2 = tl.detectProbabilitiesMC(L, e_quenching, e_quenching2, t1, evenement, extDT, measTime)
            else:
                efficiency0_S, efficiency0_D, efficiency0_T, efficiency0_AB, efficiency0_BC, efficiency0_AC, efficiency0_D2 = tl.detectProbabilities(L, e_quenching, e_quenching2, t1, evenement, extDT, measTime)
            efficiency_S.append(efficiency0_S)
            efficiency_T.append(efficiency0_T)
            efficiency_D.append(efficiency0_D)
            efficiency_AB.append(efficiency0_AB)
            efficiency_BC.append(efficiency0_BC)
            efficiency_AC.append(efficiency0_AC)            
            efficiency_D2.append(efficiency0_D2)
         

            """
            ==========================================================
            V.1 record detection probabilities in file
            ==========================================================
            """
            if record:
                with open(recfile4, "a") as file:
                    if symm:
                        file.write(f"{idec} {efficiency_S[-1]} {efficiency_D[-1]} {efficiency_T[-1]}\n")
                    else:
                        file.write(f"{idec} {efficiency_S[-1]} {efficiency_D[-1]} {efficiency_T[-1]} {efficiency_AB[-1]} {efficiency_BC[-1]} {efficiency_AC[-1]}\n")
            
        '''
        ====================
        VI. CALCULATION OF THE FINAL ESTIMATORS
        ====================
        '''
        mean_efficiency_S, std_efficiency_S, mean_efficiency_D, std_efficiency_D, mean_efficiency_T, std_efficiency_T, mean_efficiency_AB, std_efficiency_AB, mean_efficiency_BC, std_efficiency_BC, mean_efficiency_AC, std_efficiency_AC, mean_efficiency_D2, std_efficiency_D2 = tl.efficienciesEstimates(efficiency_S, efficiency_D, efficiency_T, efficiency_AB, efficiency_BC, efficiency_AC, efficiency_D2, N)
                        
        if mode == "eff":
            return mean_efficiency_S, std_efficiency_S, mean_efficiency_D, std_efficiency_D, mean_efficiency_T, std_efficiency_T, mean_efficiency_AB, std_efficiency_AB, mean_efficiency_BC, std_efficiency_BC, mean_efficiency_AC, std_efficiency_AC, mean_efficiency_D2, std_efficiency_D2
        if mode =="dis":
            return efficiency_S, efficiency_D, efficiency_T, efficiency_D2

def objectFct(L, TD, Rad, pmf_1, N, kB, V):
    """
    Objective function to minimize in order the estimate the detection efficiencies based on the measurements.

    Parameters
    ----------
    L : float or tuple
        If L is float, then L is the global free parameter. If L is tuple, then L is a triplet of free parameters. unit keV-1
    TD : float or tuple
        measurements. If TD is float, then TD is the measured TDCR parameter. If TD is tuple, then TD must contain the global TDCR parameter followed by specific ones (T/B, T/AB, T/BC, T/AC)
    Rad : string
        List of radionuclides (eg. "H-3, Co-60").
    pmf_1 : string
        list of probability of each radionuclide (eg. "0.8, 0.2").
    N : integer
        Number of Monte-Carlo trials. recommanded N>10000 (see JCGM 101). Not applied in the case of the analytical model.
    kB : float
        Birks constant in cm/keV.
    V : float
        volume of the scintillator in ml.

    Returns
    -------
    res : float
        The residual value.

    """
    
    if isinstance(TD, (tuple, list)):
        symm = False
    else:
        symm = True
        
    if symm:
        eff_model = TDCRPy(L, Rad, pmf_1, N, kB, V, readRecHist = True)
        TDCR_calcul = eff_model[4]/eff_model[2]
        res=(TDCR_calcul-TD)**2
    else:
        eff_model = TDCRPy(L, Rad, pmf_1, N, kB, V, readRecHist = True)
        TAB_calcul = eff_model[4]/eff_model[6]
        TBC_calcul = eff_model[4]/eff_model[8]
        TAC_calcul = eff_model[4]/eff_model[10]
        res=(TD[1]-TAB_calcul)**2+(TD[2]-TBC_calcul)**2+(TD[3]-TAC_calcul)**2        
    
    return res

def eff(TD, Rad, pmf_1, kB, V, N=10000, L=1, maxiter=20, xatol=1e-7, disp=False, Lbounds=[0.1, 10]):
    """
    Caclulation of the efficiency of a TDCR system based on the model TDCRPy.
    This function includes optimization procedures from scipy.

    Parameters
    ----------
    TD : float or tuple
        measurements. If TD is float, then TD is the measured TDCR parameter. If TD is tuple, then TD must contain the global TDCR parameter followed by specific ones (T/B, T/AB, T/BC, T/AC)
    Rad : string
        List of radionuclides (eg. "H-3, Co-60").
    pmf_1 : string
        list of probability of each radionuclide (eg. "0.8, 0.2").
    kB : float
        Birks constant in cm/keV.
    V : float
        volume of the scintillator in ml.
    N : interger, optional
        number of Monte-Carlo trials. The default is 10000.
    maxiter : interger, optional
        maximum number of iterations of the optimization procedures
    xatol : float
        convergence parameter of the Nelder Mead optimisation
    disp : Boolean
        to display detailed results of the procedure. Default is False.

    Returns
    -------
    L0 : float
        global free parameter.
    L : tuple
        free parameters (relevant for the asymetric model).
    eff_S : float
        counting efficiency of single events.
    u_eff_S : float
        standard uncertainty of eff_S.
    eff_D : float
        counting efficiency of double coincidences.
    u_eff_D : float
        standard uncertainty of eff_D.
    eff_T : float
        counting efficiency of triple coincidences.
    u_eff_T : float
        standard uncertainty of eff_T.
    eff_AB : float
        counting efficiency of coincidences AB.
    u_eff_AB : float
        standard uncertainty of eff_AB.
    eff_BC : float
        counting efficiency of coincidences BC.
    u_eff_BC : float
        standard uncertainty of eff_BC.
    eff_AC : float
        counting efficiency of coincidences AC.
    u_eff_AC : float
        standard uncertainty of eff_AC.    
    eff_D : float
        counting efficiency of double coincidences in C/N configuation (not relevant).
    u_eff_D : float
        standard uncertainty of eff_D in C/N configuation (not relevant).

    """
    if isinstance(TD, (tuple, list)):
        symm = False
    else:
        symm = True
    
    TDCRPy(L, Rad, pmf_1, N, kB, V, record = True)
       
    if symm: r=opt.minimize_scalar(objectFct, args=(TD, Rad, pmf_1, N, kB, V), method='bounded', bounds = (Lbounds[0], Lbounds[1]), options={'disp': disp, 'maxiter':maxiter})
    else: r=opt.minimize_scalar(objectFct, args=(TD[0], Rad, pmf_1, N, kB, V), method='bounded', bounds = (Lbounds[0], Lbounds[1]), options={'disp': disp, 'maxiter':maxiter})
    L0=r.x
    L=(L0, L0, L0)
    print(f"global free parameter = {L0} keV-1")
    
    if not symm:
        r=opt.minimize(objectFct, L, args=(TD, Rad, pmf_1, N, kB, V), method='nelder-mead',options={'xatol': xatol, 'disp': disp, 'maxiter':maxiter})
        L=r.x
        print(f"free parameters = {L} keV-1")   

    if symm: out=TDCRPy(L0, Rad, pmf_1, N, kB, V, readRecHist = True)
    else: out=TDCRPy(L, Rad, pmf_1, N, kB, V, readRecHist = True)
    eff_S = out[0]
    u_eff_S = out[1]
    eff_D = out[2]
    u_eff_D = out[3]
    eff_T = out[4]
    u_eff_T = out[5]
    eff_AB = out[6]
    u_eff_AB = out[7]
    eff_BC = out[8]
    u_eff_BC = out[9]
    eff_AC = out[10]
    u_eff_AC = out[11]
    eff_D2 = out[12]
    u_eff_D2 = out[13]
    
    return L0, L, eff_S, u_eff_S, eff_D, u_eff_D, eff_T, u_eff_T, eff_AB, u_eff_AB, eff_BC, u_eff_BC, eff_AC, u_eff_AC, eff_D2, u_eff_D2


def effA(TD, Rad, pmf_1, kB, V, L=1, maxiter=20, xatol=1e-7, disp=False, Lbounds=[0.1, 10]):
    """
    Caclulation of the efficiency of a TDCR system based on the model TDCRPy (analytical model).
    This function includes optimization procedures from scipy.

    Parameters
    ----------
    TD : float or tuple
        measurements. If TD is float, then TD is the measured TDCR parameter. If TD is tuple, then TD must contain the global TDCR parameter followed by specific ones (T/D, T/AB, T/BC, T/AC)
    Rad : string
        List of radionuclides (eg. "H-3, Co-60").
    pmf_1 : string
        list of probability of each radionuclide (eg. "0.8, 0.2").
    kB : float
        Birks constant in cm/keV.
    V : float
        volume of the scintillator in ml.
    maxiter : interger, optional
        maximum number of iterations of the optimization procedures
    xatol : float
        convergence parameter of the Nelder Mead optimisation
    disp : Boolean
        to display detailed results of the procedure. Default is False.

    Returns
    -------
    L0 : float
        global free parameter.
    L : tuple
        free parameters (relevant for the asymetric model).
    eff_S : float
        counting efficiency of single events.
    u_eff_S : float
        standard uncertainty of eff_S.
    eff_D : float
        counting efficiency of double coincidences.
    u_eff_D : float
        standard uncertainty of eff_D.
    eff_T : float
        counting efficiency of triple coincidences.
    u_eff_T : float
        standard uncertainty of eff_T.
    eff_AB : float
        counting efficiency of coincidences AB.
    u_eff_AB : float
        standard uncertainty of eff_AB.
    eff_BC : float
        counting efficiency of coincidences BC.
    u_eff_BC : float
        standard uncertainty of eff_BC.
    eff_AC : float
        counting efficiency of coincidences AC.
    u_eff_AC : float
        standard uncertainty of eff_AC.    
    eff_D : float
        counting efficiency of double coincidences in C/N configuation (not relevant).
    u_eff_D : float
        standard uncertainty of eff_D in C/N configuation (not relevant).

    """
    if isinstance(TD, (tuple, list)):
        symm = False
    else:
        symm = True

    if symm: r=opt.minimize_scalar(tl.modelAnalytical, args=(TD, TD, TD, TD, Rad, kB, V, "res", 1e3), method='bounded', bounds = (Lbounds[0], Lbounds[1]), options={'disp': disp, 'maxiter':maxiter})
    else: r=opt.minimize_scalar(tl.modelAnalytical, args=(TD[0], TD[1], TD[2], TD[3], Rad, kB, V, "res", 1e3), method='bounded', bounds = (Lbounds[0], Lbounds[1]), options={'disp': disp, 'maxiter':maxiter})
    L0=r.x
    L=(L0, L0, L0)
    print(f"global free parameter = {L0} keV-1")
    
    if not symm:
        r=opt.minimize(tl.modelAnalytical, L, args=(TD[0], TD[1], TD[2], TD[3], Rad, kB, V, "res", 1e3), method='nelder-mead',options={'xatol': xatol, 'disp': disp, 'maxiter':maxiter})
        L=r.x
        print(f"free parameters = {L} keV-1")   

    if symm: out=tl.modelAnalytical(L, TD, TD, TD, TD, Rad, kB, V, "eff", 1e3)
    else: out=tl.modelAnalytical(L, TD[0], TD[1], TD[2], TD[3], Rad, kB, V, "eff", 1e3)
    eff_S = out[0]
    eff_D = out[1]
    eff_T = out[2]
    # u_eff_S = out[1]
    # eff_D = out[2]
    # u_eff_D = out[3]
    # eff_T = out[4]
    # u_eff_T = out[5]
    # eff_AB = out[6]
    # u_eff_AB = out[7]
    # eff_BC = out[8]
    # u_eff_BC = out[9]
    # eff_AC = out[10]
    # u_eff_AC = out[11]
    
    return L0, L, eff_S, eff_D, eff_T #, u_eff_T, eff_AB, u_eff_AB, eff_BC, u_eff_BC, eff_AC, u_eff_AC, u_eff_S, u_eff_D, 



# mode = "eff"                # ask for efficiency calculation
# Rad="Fe-55"                 # radionuclides
# pmf_1="1"                   # relatives fractions of the radionulides
# N = 1000                    # number of Monte Carlo trials
# kB =1.0e-5                  # Birks constant in cm keV-1
# V = 10                      # volume of scintillator in mL
# L=np.logspace(-1,2,num=100) # free parameter in keV-1

# # Record decay histories in temporary files
# TDCRPy(L[0], Rad, pmf_1, N, kB, V, mode, barp=False, record=True, fullMC=False)
# # TDCRPy(100, Rad, pmf_1, N, kB, V, mode, barp=False, record=True, fullMC=True)

# effS, u_effS, effD, u_effD, effT, u_effT, effD2, u_effD2 = [], [],[], [],[], [], [], []
# for l in tqdm(L, desc="free parameters ", unit=" iterations"):
#   out = TDCRPy(l, Rad, pmf_1, N, kB, V, mode, readRecHist=True, fullMC=False)
#   effS.append(out[2])
#   u_effS.append(out[3])
#   effD.append(out[2])
#   u_effD.append(out[3])
#   effT.append(out[4])
#   u_effT.append(out[5])
#   effD2.append(out[12])
#   u_effD2.append(out[13])

# effS=np.asarray(effS)
# effT=np.asarray(effT)
# effD=np.asarray(effD)
# effD2=np.asarray(effD2)
# u_effS=np.asarray(u_effS)
# u_effT=np.asarray(u_effT)
# u_effD=np.asarray(u_effD)
# u_effD2=np.asarray(u_effD2)

# tdcr=effT/effD
# u_tdcr=np.sqrt(u_effD**2*effT**2/effD**4+u_effT**2/effD**2)

# import matplotlib.pyplot as plt
# plt.figure("efficiency vs free parameter")
# plt.clf()
# plt.errorbar(L,effD,yerr=u_effD,fmt="-k",label="double coincidences")
# plt.errorbar(L,effT,yerr=u_effT,fmt="-r",label="triple coincidences")
# plt.errorbar(L,effD2,yerr=u_effD2,fmt="-g",label="double coincidences (CIEMAT/NIST)")
# plt.xscale('log')
# plt.xlabel(r'$L$ /keV$^{-1}$', fontsize=14)
# plt.ylabel(r'$\epsilon$', fontsize=14)
# plt.legend()

# plt.figure("efficiency vs TDCR")
# plt.clf()
# plt.errorbar(tdcr,effD,xerr=u_tdcr,yerr=u_effD,fmt="-k")
# #plt.xscale('log')
# plt.xlabel(r'$R_T/R_D$', fontsize=14)
# plt.ylabel(r'$\epsilon_{D}$', fontsize=14)
# plt.show()


# L = 1
# # L = (1.1, 1.05, 1.15)
# # TD = 0.977667386529166
# # TD = (0.9767359812638453, 0.9925429293804757, 0.991829757077315, 0.9919970813639295) # source 1
# # TD = (0.9768862920127371, 0.9928478299182348, 0.9912531441227223, 0.9924249578285456) # source 2
# # TD = (0.9769014488454436, 0.9918130431206161, 0.9920156754198314, 0.9927119011073454) # source 3
# TD = (0.9764032345164899, 0.9928417189012709, 0.9911455450383777, 0.9920402844839974) # source 4
# Rad="Tc-99"
# pmf_1="1"
# N = 10000
# kB =1.4e-5
# V = 16
# mode = "eff"


# # # out = TDCRPy(L, Rad, pmf_1, N, kB, V, Display = False, record = True, readRecHist = False)
# # # print("result", out)
# # # out = TDCRPy(L, Rad, pmf_1, N, kB, V, Display = False, record = False, readRecHist = True)
# # # print("result", out)

# outS = eff(TD, Rad, pmf_1, kB, V, N=10000, L=1, maxiter=20, xatol=1e-7)
# outA = effA(TD, Rad, pmf_1, kB, V, L=1, maxiter=20, xatol=1e-7)
# print(outS)
# print(outA)