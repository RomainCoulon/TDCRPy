# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:04:46 2023

Library of function of the TDCRpy code

@author: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

### IMPORT Python Module
import urllib.request as rq
import importlib.resources
import numpy as np
import zipfile as zf
import time
import re
absolutePath = False


def TicTocGenerator():
    """
    Generator that returns time differences
    """
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    """
    Prints the time difference yielded by generator instance TicToc
    """
    
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    """
    Records a time in TicToc, marks the beginning of a time interval
    """
    toc(False)

def normalise(p_x):
    p_array = np.array(p_x)
    if len(p_x)>1:
        p_somme = sum(p_array)
        if p_somme>0.0:
            p_array = p_array/p_somme
    else:
        p_somme = p_x[0]
        p_array = p_array/p_somme
    p = list(p_array)
    return p
#print(normalise([1, 0.9974641, 0]))

def sampling(p_x):
    """
    This function aims to sample in a pdf or a pmf

    Parameters
    ----------
    p_x : float vector
        Probability Density (or mass) Function (PDF or PMF) of the random variable x.


    Returns
    -------
    i : integer
        index in x pointing the sampled value of the random variable X.
    """

    cf = np.cumsum(p_x) # Cummulative Density (or mass) Function (CDF or CMF)
    trial = float(np.random.rand(1)) # trial ~ U(0,1)
    
    for i, p in enumerate(cf):
        if p> trial: break
    return i


def readPenNuc(rad):
    """
    This function reads the PenNuc files published in the DDEP web page
    http://www.lnhb.fr/donnees-nucleaires/donnees-nucleaires-tableau/
    
    Parameters
    ----------
    rad : string
        Indentifier of the radionuclide (e.g. Na-22).

    Returns
    -------
    out : list
        list of elements in list describing the decay of the radionuclide with,
        *particle* is the list of particle of the decay branches 
        *p_branch* is the list of probabilty of the decay branches
        *e_branch* is the list of energy of the decay branches
        *LevelDaughter* is the list of energy levels of the daughter nucleus following the decay branch
        *levelNumber* is the list of number of the energy level of the daughter nucleus following the decay branch
        *prob* is the list of probabilty of isomeric transitions
        *levelEnergy* is the energy levels of the isomeric transitions
        *transitionType* is the type of isomeric transitions
        *e_trans* is the energy of the isomeric transitions
        *next_level* is the possible next energy level following a given isomeric transtion 
        *Q_value_vec[dd]* is the Q-value of the decay when decaying to the daughter of index dd
        *Daughter_vec[dd]* is the vecteur of daugher nuclei
        *Pdaughter_vec[dd]* is the probability of daughter nuclei
    """
    
    url = "http://www.lnhb.fr/nuclides/"+rad+".PenNuc.txt"
    file = rq.urlopen(url)
    
    # Format the data 
    decayData = []
    for line in file:
      decayData.append(line.decode("utf-8"))
      if "NDA " in decayData[-1]:    decayData[-1] = decayData[-1].replace("NDA ", "NDA; ") # number of daughter
      if "DAU " in decayData[-1]:    decayData[-1] = decayData[-1].replace("DAU ", "DAU; ") # daughter
      if "DDE " in decayData[-1]:    decayData[-1] = decayData[-1].replace("DDE ", "DDE; ") # daughter description : probability of disintegration to NDA, uncertainty, number of excited levels of NDA; number of branche to NDA
      if "Q " in decayData[-1]:      decayData[-1] = decayData[-1].replace("Q ", "Q; ")     # total energy of the branch, uncertainty 
      if "ALP " in decayData[-1]:    decayData[-1] = decayData[-1].replace("ALP ", "ALP; ") # type of the disintegration = alpha
      if "CK " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CK ", "CK; ") # type of the disintegration = Electron Capture K shell
      if "CL " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CL ", "CL; ") # type of the disintegration = Electron Capture L shell 1
      if "CL1 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CL1 ", "CL1; ") # type of the disintegration = Electron Capture L shell 1
      if "CL2 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CL2 ", "CL2; ") # type of the disintegration = Electron Capture L shell 2
      if "CM " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CM ", "CM; ") # type of the disintegration = Electron Capture M shell
      if "CN " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CN ", "CN; ") # type of the disintegration = Electron Capture N shell
      if "BEM " in decayData[-1]:    decayData[-1] = decayData[-1].replace("BEM ", "BEM; ") # type disintegration-beta- :branching ratio, uncertainty, level fed in daughter, energy, uncertainty,prohibition factor for beta emission
      if "BEP " in decayData[-1]:    decayData[-1] = decayData[-1].replace("BEP ", "BEP; ") # type disintegration-beta+ :branching ratio, uncertainty, level fed in daughter, energy, uncertainty,prohibition factor for beta emission
      if "LED " in decayData[-1]:    decayData[-1] = decayData[-1].replace("LED ", "LED; ") # level description:energy of level,uncertainty, number of transitions that depopulate this level,level time,uncertainty,level number
      if "GA " in decayData[-1]:     decayData[-1] = decayData[-1].replace("GA ", "GA; ")   # type transition: gamma
      if "EK " in decayData[-1]:     decayData[-1] = decayData[-1].replace("EK ", "EK; ")   # type transition: internal conversion -- electron K
      if "EL1 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("EL1 ", "EL1; ") # type transition: internal conversion -- electron L1
      if "EL2 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("EL2 ", "EL2; ") # type transition: internal conversion -- electron L2
      if "EL3 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("EL3 ", "EL3; ") # type transition: internal conversion -- electron L3
      if "EM " in decayData[-1]:     decayData[-1] = decayData[-1].replace("EM ", "EM; ")   # type transition: internal conversion -- electron M
      if "EN " in decayData[-1]:     decayData[-1] = decayData[-1].replace("EN ", "EN; ")   # type transition: internal conversion -- electron N
      if "COM" in decayData[-1]:     decayData[-1] = decayData[-1].replace("COM ", "COM; ")   # commentaire
      decayData[-1]=decayData[-1].split(";")
      if "\r\n" in decayData[-1][-1]:  decayData[-1][-1] = decayData[-1][-1].replace("\r\n", " ")
    
   
      
    """
    LOOP IN NDAUGH (Daughters)
    """
    Daughter_vec = []; Pdaughter_vec = []; Q_value_vec = []
    for d in decayData:
      if d[0] == "NDA": Ndaughter = int(d[1])   # Read the number of daughter
      if d[0] == "DAU" :
          d[1] = d[1].replace(" ","")
          Daughter_vec.append(d[1])        # Read the daughter
      if d[0] == "DDE": Pdaughter_vec.append(float(d[1])) # Read the probability of disintegration to daughter
      if d[0] == "Q":   Q_value_vec.append(float(d[1]))   # Read the Q-value
    
    out = []
    for dd in range(Ndaughter):
        
        # index_Daughter = sampling(np.asarray(Pdaughter_vec)/sum(np.asarray(Pdaughter_vec)))
        # Daughter = Daughter_vec[index_Daughter]
        # Pdaughter =  Pdaughter_vec[index_Daughter]
        # Q_value = Q_value_vec[index_Daughter]    
         
        
        """
        LOOP IN NBRANCH (Branches of each daughter)
        """
        # Probablity vector of the decay branch
        particle=[]; LevelDaughter = []; p_branch = []; u_p_branch = []; e_branch = []; u_e_branch = []
        daughterFlag = True
        countBranch = 0
        for d in decayData:
            if d[0] == "DAU" and d[1] == Daughter_vec[dd]:
                daughterFlag = True
            if d[0] == "DAU" and d[1] != Daughter_vec[dd]:
                daughterFlag = False
            
            if daughterFlag:
                # particle_i = []; p_branch_i=[]; u_p_branch_i=[];  e_branch_i=[];   u_e_branch_i=[];  LevelDaughter_i=[]
                if d[0] =="COM":
                    if "Branch" in d[1] or "Level" in d[1]:
                        print(d[1])
                        if countBranch % 2 == 1:
                            p_branch.append(p_branch_i); u_p_branch.append(u_p_branch_i); # Branching ratio of the decay branch
                            e_branch.append(e_branch_i); u_e_branch.append(u_e_branch_i); # Energy of the decay branch (kinetic energy of the particle)
                            LevelDaughter.append(LevelDaughter_i); particle.append(particle_i)
                            
                            particle_i = []; p_branch_i=[]; u_p_branch_i=[];  e_branch_i=[];   u_e_branch_i=[];  LevelDaughter_i=[]
                        else:
                            particle_i = []; p_branch_i=[]; u_p_branch_i=[];  e_branch_i=[];   u_e_branch_i=[];  LevelDaughter_i=[]
                    
                    
                        countBranch += 1    
                    # p_branch=
                
                
                if d[0] == "ALP" or d[0] == "BEM" or d[0] == "BEP" or d[0] == "CK" or d[0] == "CL" or d[0] == "CL1" or d[0] == "CL2" or d[0] == "CM" or d[0] == "CN":                         # Read information on the decay branch
                    if d[0] == "ALP": particle_i.append("alpha")
                    if d[0] == "BEP": particle_i.append("beta+")
                    if d[0] == "BEM": particle_i.append("beta")
                    if d[0] == "CK": particle_i.append("Atom_K")
                    if d[0] == "CL": particle_i.append("Atom_L")
                    if d[0] == "CL1": particle_i.append("Atom_L1")
                    if d[0] == "CL2": particle_i.append("Atom_L2")
                    if d[0] == "CM": particle_i.append("Atom_M")
                    if d[0] == "CN": particle_i.append("Atom_N")
                    if d[2] == "  ": d[2]=0
            
                    p_branch_i.append(float(d[1])); u_p_branch_i.append(float(d[2])); # Branching ratio of the decay branch
                    e_branch_i.append(float(d[4])); u_e_branch_i.append(float(d[5])); # Energy of the decay branch (kinetic energy of the particle)
                    LevelDaughter_i.append(int(d[3]))                               # Level fed in daughter
                    
                
        
        """
        LOOP IN NLEVEL (Levels for each daughter, starting in NLEVEL, ending in 1)
        """
        levelEnergy = []; fromBranch = []; lineL=[]; listTran = []; levelNumber = []
        transitionType = []; prob = []; u_prob = []; e_trans = []; u_e_trans = []; next_level = []
        transitionType_i = []; prob_i = []; u_prob_i = []; e_trans_i = []; u_e_trans_i = []; next_level_i = []; levelEnergy_i = []
        for d in decayData:
          if d[0] == "DAU" and d[1] == Daughter_vec[dd]:
                daughterFlag = True
          if d[0] == "DAU" and d[1] != Daughter_vec[dd]:
                daughterFlag = False
          if daughterFlag:
              if d[0] == "LED":
                # record transition details of the previous level
                 transitionType.append(transitionType_i); prob.append(prob_i); u_prob.append(u_prob_i); e_trans.append(e_trans_i); u_e_trans.append(u_e_trans_i); next_level.append(next_level_i);levelEnergy.append(levelEnergy_i)
            #    if int(d[6]) in LevelDaughter :  # Read the information on the possible energy levels after the emission of the alpha particle 
                 levelEnergy_i.append(float(d[1]))   #levelEnergy.append(float(d[1]));     # Energie (rounded) of the level
                 listTran.append(int(d[3]))           # Number of transitions that depopulate this level
                 levelNumber.append(int(d[6]))        # Level number
                 transitionType_i = []; next_level_i = []; levelEnergy_i = [];prob_i = []; u_prob_i = []; e_trans_i = []; u_e_trans_i = [];
                  
              """
              LOOP IN NTRANS (Transitions depopulating this level)
              """
              if d[0] == "GA" or d[0] == "EK" or d[0] == "EL1" or d[0] == "EL2" or d[0] == "EL3" or d[0] == "EM" or d[0] == "EN":
                  if d[1] == '  ' or d[1] == '   ': d[1] = 0
                  if d[2] == '  ': d[2] = 0
                  if d[4] == '  ': d[4] = 0
                  transitionType_i.append(d[0])                                          # Read the type of transtion
                  prob_i.append(float(d[1])); u_prob_i.append(float(d[2]))                 # Read the emission probability of the transition
                  e_trans_i.append(float(d[3])); u_e_trans_i.append(float(d[4]))           # Read the energy of the transition
                  next_level_i.append(int(d[5]))                                         # Read the level fed by this transition
            
        transitionType.append(transitionType_i); prob.append(prob_i); u_prob.append(u_prob_i); e_trans.append(e_trans_i); u_e_trans.append(u_e_trans_i); next_level.append(next_level_i);levelEnergy.append(levelEnergy_i)
    
        out.append([particle, p_branch, e_branch, LevelDaughter, levelNumber, prob, levelEnergy, transitionType, e_trans, next_level, Q_value_vec[dd], Daughter_vec[dd], Pdaughter_vec[dd]])
    return out
#out = readPenNuc('Am-242')
#print(len(out),out)
# tic()
# readPenNuc("Co-60")
# toc() # 0.016 s

if absolutePath: file_pennuc = 'G:\Python_modules\Jialin\Code\decayData\\All-nuclides_PenNuc.zip'
else:
    # file_pennuc = "decayData//All-nuclides_PenNuc.zip"
    with importlib.resources.path('tdcrpy', 'decayData') as data_path:
        file_pennuc = data_path / "All-nuclides_PenNuc.zip"

z_PenNuc = zf.ZipFile(file_pennuc)

def readPenNuc2(rad,z1=z_PenNuc):
    '''
     =========
     PARAMETRE
     =========
     rad -- type: str (par exemple: "Am-241") -- radionucléide 

     ======
     RETURN
     ======
     daughter -- indice 0 -- des noyaux fils -- len = nb de noyaux fils
     prob_daug -- indice 1 -- des probabilités de noyaux fils -- len = nb de noyaux fils
     energy_Q -- indice 2 -- des énergies de désintégrations -- len = nb de noyaux fils

     desin_type_tot -- indice 3 -- des types de désintégrations/particules émis
         len = nb de noyaux fils 
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des désintégrations possibles de chaque branch -- len de sous-list de sous-list = nb de type de désintégrations de chaque branch

     desin_energy_tot -- indice 4 -- des énergies de désintégrations/énergies de patricules émis
         len = nb de noyaux fils 
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des énergies de désintégrations possibles de chaque branch -- len de sous-list de sous-list = nb de type de désintégrations de chaque branch

     desin_prob_tot -- indice 5 -- des probabilités de désintégrations
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des probabilités de désintégrations possibles de chaque branch -- len de sous-list de sous-list = nb de type de désintégrations de chaque branch
     
     desin_level_tot -- indice 6 -- des niveaux atteints après des désintégrations
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des niveaux après des désintégrations de chaque branch -- len de sous-list de sous-list = nb de type de désintégrations de chaque branch
     
     prob_branch_tot -- indice 7 -- probabilités de chaque branch
         len = nb de noyaux fils
         sous-list -- des probabilités de branchs de noyau fil -- len de sous-list = nb de branch de chaque fil
         
     tran_type_tot -- indice 8 -- transitions possibles 
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des transitions possibles de chaque branch -- len de sous-list de sous-list = nb de type de transitions de chaque branch
     
     tran_energy_tot -- indice 9 -- énergies de transitions
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des énergies de transitions possibles de chaque branch -- len de sous-list de sous-list = nb de type de transitions de chaque branch
     
     tran_prob_tot -- indice 10 -- probabilités de transitions
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des probabilités de transitions possibles de chaque branch -- len de sous-list de sous-list = nb de type de transitions de chaque branch
     
     tran_level_tot -- indice 11 -- niveaux de branch correspondants
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des niveaux de chaque branch avant des transitions -- len de sous-list de sous-list = 1
     
     tran_level_end_tot -- indice 12 -- niveaux après des transitions
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des niveaux après des transitions de chaque branch -- len de sous-list de sous-list = nb de type de transitions de chaque branch
     
     level_energy_tot -- indice 13 -- énergies de niveaux
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des énergies de niveaux de chaque branch -- len de sous-list de sous-list = 1

     prob_tran_tot -- indice 14 -- la somme de transition de chaque branch
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des énergies de niveaux de chaque branch -- len de sous-list de sous-list = 1
     '''
    doc = rad + ".PenNuc.txt"
    with z1.open(doc) as file_P:
        decayData = file_P.readlines()

    for i in range(np.size(decayData)):
        decayData[i] = str(decayData[i])
        decayData[i] = decayData[i].replace("b'","")
        decayData[i] = decayData[i].replace("\\r\\n","")
        decayData[i] = decayData[i].replace("'","") 

    for il in range(len(decayData)):
        if "NDA " in decayData[il]: decayData[il] = decayData[il].replace("NDA ","NDA; ") 
        if "DAU " in decayData[il]: decayData[il] = decayData[il].replace("DAU ","DAU; ") 
        if "DDE " in decayData[il]: decayData[il] = decayData[il].replace("DDE ","DDE; ")
        if "Q " in decayData[il]: decayData[il] = decayData[il].replace("Q ","Q; ")
        if "ALP " in decayData[il]: decayData[il] = decayData[il].replace("ALP ","ALP; ")
        if "CK " in decayData[il]: decayData[il] = decayData[il].replace("CK ","CK; ")
        if "CL " in decayData[il]: decayData[il] = decayData[il].replace("CL ","CL; ")
        if "CO " in decayData[il]: decayData[il] = decayData[il].replace("CO ","CO; ")
        if "CL1 " in decayData[il]: decayData[il] = decayData[il].replace("CL1 ","CL1; ")
        if "CL2 " in decayData[il]: decayData[il] = decayData[il].replace("CL2 ","CL2; ")
        if "CL3 " in decayData[il]: decayData[il] = decayData[il].replace("CL3 ","CL3; ")
        if "CN " in decayData[il]: decayData[il] = decayData[il].replace("CN ","CN; ")
        if "CM " in decayData[il]: decayData[il] = decayData[il].replace("CM ","CM; ")
        if "BEM " in decayData[il]: decayData[il] = decayData[il].replace("BEM ","BEM; ")
        if "BEP " in decayData[il]: decayData[il] = decayData[il].replace("BEP ","BEP; ")
        if "LED " in decayData[il]: decayData[il] = decayData[il].replace("LED ","LED; ")
        if "GA " in decayData[il]: decayData[il] = decayData[il].replace("GA ","GA; ")
        if "EK " in decayData[il]: decayData[il] = decayData[il].replace("EK ","EK; ")
        if "EL " in decayData[il]: decayData[il] = decayData[il].replace("EL ","EL; ")
        if "EL1 " in decayData[il]: decayData[il] = decayData[il].replace("EL1 ","EL1; ")
        if "EL2 " in decayData[il]: decayData[il] = decayData[il].replace("EL2 ","EL2; ")
        if "EL3 " in decayData[il]: decayData[il] = decayData[il].replace("EL3 ","EL3; ")
        if "EM " in decayData[il]: decayData[il] = decayData[il].replace("EM ","EM; ")
        if "EN " in decayData[il]: decayData[il] = decayData[il].replace("EN ","EN; ")
        if "COM " in decayData[il]: decayData[il] = decayData[il].replace("COM ","COM; ")
        decayData[il] = decayData[il].split(';')

    for a1 in decayData:
        for a2 in range(len(a1)):
            a1[a2] = a1[a2].strip()

    '''
     ========================
     Repérer chaque noyau fil
     ========================

     daughter -- noyau(x) fil(s)
     posi_daug -- l'indice de démarcation de noyau fil
     posi_branch -- l'indice de démarcation de chaque branch
     posi_tran -- l'indice de démarcation de transition
     prob_daug -- probabilité de produire des noyaux fils
     nb_branch -- nombre de branch possible au dessus de l'état fonda (n>0)
     energy_Q -- l'énergie de désintégration

    '''    
    daughter = [];posi_daug = [];prob_daug=[];nb_branch=[];energy_Q=[];
    end = len(decayData)
    for indice,line in enumerate(decayData):
        if "NDA" == line[0]:
            nb_daug = int(line[1])
        if "DAU" == line[0]:
            daughter.append(line[1].replace(" ",""))
        if "COM" == line[0] and "Daughter" in line[1]:
            posi_daug.append(indice)
        if "Q" == line[0]:
            energy_Q.append(float(line[1]))
        if "DDE" == line[0]:
            prob_daug.append(float(line[1]))
            nb_branch.append(int(line[-2]))
    '''
     ==========
     LOOP START
     ==========

    ''' 
    posi_end=[]
    desin_type_tot=[];desin_energy_tot=[];desin_prob_tot=[];desin_level_tot=[]
    tran_type_tot=[];tran_energy_tot=[];tran_prob_tot=[];tran_level_end_tot=[]; 
    tran_level_tot=[];level_energy_tot=[]
    prob_branch_tot=[];prob_tran_tot=[] 

    '''
     =============
     LOOP DAUGHTER 
     =============

    '''
    for i1 in range(nb_daug):
        start_p = posi_daug[i1]
        if i1+1 == nb_daug:
            end_p = end
        else:
            end_p = posi_daug[i1+1]

        posi_end_i = []
        for i2 in range(start_p,end_p):
            if "COM" == decayData[i2][0] and "Branch" in decayData[i2][1]:
                posi_end_i.append(i2)
            if "COM" == decayData[i2][0] and "Level" in decayData[i2][1]:
                posi_end_i.append(i2)
        if end_p == end:
            posi_end_i.append(end)
        else:
            posi_end_i.append(posi_daug[i1+1])
        posi_end.append(posi_end_i)

        '''
         ==========================
         LOOP Branch and Transition
         ==========================
        '''
        desin_type_daug=[];desin_energy_daug=[];desin_prob_daug=[];desin_level_daug=[];
        tran_type_daug=[];tran_energy_daug=[];tran_prob_daug=[];tran_level_end_daug=[];
        tran_level_daug=[];level_energy_daug=[]
        prob_branch_daug=[];prob_tran_daug=[]

        for i3 in range(len(posi_end_i)-1):
            start_p1 = posi_end_i[i3]
            end_p1 = posi_end_i[i3+1]
            branch = False
            transition = False
            if "COM" == decayData[start_p1][0] and "Branch" in decayData[start_p1][1]:
                branch=True
            if "COM" == decayData[start_p1][0] and "Level" in decayData[start_p1][1]:
                transition = True

            '''
             ====================================
             LOOP EACH BLOCK OF BRANCH/TRANSITION
             ====================================
            ''' 
            tran_type_b=[];tran_prob_b=[];tran_energy_b=[]; tran_level_end_b=[];
            tran_level_b=[];level_energy_b=[]
            desin_type_b=[];desin_energy_b=[];desin_prob_b=[];desin_level_b=[];
             
            for i4 in decayData[start_p1+1:end_p1]:
                if start_p1+1 == end_p1:
                    break
                if branch:
                    if "ALP" == i4[0]:
                        desin_type_b.append("alpha")
                    if "BEP" == i4[0]:
                        desin_type_b.append("beta+")
                    if "BEM" == i4[0]:
                        desin_type_b.append("beta")
                    if "CK" == i4[0]:
                        desin_type_b.append("Atom_K")
                    if "CL" == i4[0]:
                        desin_type_b.append("Atom_L")
                    if "CL1" == i4[0]:
                        desin_type_b.append("Atom_L1")
                    if "CL2" == i4[0]:
                        desin_type_b.append("Atom_L2")
                    if "CL3" == i4[0]:
                        desin_type_b.append("Atom_L3")
                    if "CM" == i4[0]:
                        desin_type_b.append("Atom_M")
                    if "CN" == i4[0]:
                        desin_type_b.append("Atom_N")
                    if "CO" == i4[0]:
                        desin_type_b.append("Atom_O")
                    desin_prob_b.append(float(i4[1]))
                    desin_level_b.append(int(i4[3]))
                    desin_energy_b.append(float(i4[4]))
                if transition:
                    if i4[1] == '' or i4[1] == ' ': i4[1] = 0
                    if len(i4)>2 and i4[2] == '': i4[2] = 0
                    if len(i4)>4 and i4[4] == '': i4[4] = 0
                    if len(i4)>5 and i4[5] == '': i4[5] = 0
                    if "LED" == i4[0]:
                        tran_level_b.append(int(i4[-1]))
                        level_energy_b.append(float(i4[1]))
                    if i4[0] == "GA" or i4[0] == "EK" or i4[0] == "EL" or i4[0] == "EL1" or i4[0] == "EL2" or i4[0] == "EL3" or i4[0] == "EM" or i4[0] == "EN":
                        #print(i4)
                        tran_type_b.append(i4[0])
                        tran_prob_b.append(float(i4[1]))
                        tran_energy_b.append(float(i4[3]))
                        tran_level_end_b.append(int(i4[5]))
            if branch:
                desin_type_daug.append(desin_type_b)
                desin_energy_daug.append(desin_energy_b)
                desin_prob_daug.append(desin_prob_b)
                desin_level_daug.append(desin_level_b)
             
            if transition:
                tran_type_daug.append(tran_type_b)
                tran_energy_daug.append(tran_energy_b)
                tran_prob_daug.append(tran_prob_b)
                tran_level_end_daug.append(tran_level_end_b)
                tran_level_daug.append(tran_level_b)
                level_energy_daug.append(level_energy_b)

            if len(desin_prob_b)>0:
                desin_prob_array = np.array(desin_prob_b)
                prob_branch_i = np.sum(desin_prob_array)
                if prob_branch_i >= 1:
                    prob_branch_i = 1
                prob_branch_daug.append(prob_branch_i)
            elif branch and len(desin_prob_b)==0:
                prob_branch_daug.append(0)

            if len(tran_prob_b)>0:
                tran_prob_array = np.array(tran_prob_b)
                prob_tran_i = np.sum(tran_prob_array)
                if prob_tran_i >= 1:
                    prob_tran_i = 1
                prob_tran_daug.append(prob_tran_i)
            elif transition and len(tran_prob_b)==0:
                prob_tran_daug.append(0)

        tran_type_daug.append([])
        tran_prob_daug.append([])
        tran_energy_daug.append([])
        tran_level_end_daug.append([])
        tran_level_daug.append([])
        level_energy_daug.append([])
        prob_tran_daug.append(0)

        desin_type_tot.append(desin_type_daug)
        desin_energy_tot.append(desin_energy_daug)
        desin_prob_tot.append(desin_prob_daug)
        desin_level_tot.append(desin_level_daug)

        tran_type_tot.append(tran_type_daug)
        tran_energy_tot.append(tran_energy_daug)
        tran_prob_tot.append(tran_prob_daug)
        tran_level_end_tot.append(tran_level_end_daug)
        tran_level_tot.append(tran_level_daug)
        level_energy_tot.append(level_energy_daug)
        prob_branch_tot.append(prob_branch_daug)
        prob_tran_tot.append(prob_tran_daug)
    out = [daughter,prob_daug,energy_Q,desin_type_tot,desin_energy_tot,desin_prob_tot,desin_level_tot,prob_branch_tot,tran_type_tot,tran_energy_tot,tran_prob_tot,tran_level_tot,tran_level_end_tot,level_energy_tot,prob_tran_tot]
    return out
#tic()
#o = readPenNuc2("Sb-127")
#toc()
#print(o[-4])


def readPenNuc1(rad):
     '''
     =========
     PARAMETRE
     =========
     rad -- type: str (par exemple: "Am-241") -- radionucléide 

     ======
     RETURN
     ======
     daughter -- indice 0 -- des noyaux fils -- len = nb de noyaux fils
     prob_daug -- indice 1 -- des probabilités de noyaux fils -- len = nb de noyaux fils
     energy_Q -- indice 2 -- des énergies de désintégrations -- len = nb de noyaux fils

     desin_type_tot -- indice 3 -- des types de désintégrations/particules émis
         len = nb de noyaux fils 
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des désintégrations possibles de chaque branch -- len de sous-list de sous-list = nb de type de désintégrations de chaque branch

     desin_energy_tot -- indice 4 -- des énergies de désintégrations/énergies de patricules émis
         len = nb de noyaux fils 
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des énergies de désintégrations possibles de chaque branch -- len de sous-list de sous-list = nb de type de désintégrations de chaque branch

     desin_prob_tot -- indice 5 -- des probabilités de désintégrations
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des probabilités de désintégrations possibles de chaque branch -- len de sous-list de sous-list = nb de type de désintégrations de chaque branch
     
     desin_level_tot -- indice 6 -- des niveaux atteints après des désintégrations
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des niveaux après des désintégrations de chaque branch -- len de sous-list de sous-list = nb de type de désintégrations de chaque branch
     
     prob_branch_tot -- indice 7 -- probabilités de chaque branch
         len = nb de noyaux fils
         sous-list -- des probabilités de branchs de noyau fil -- len de sous-list = nb de branch de chaque fil
         
     tran_type_tot -- indice 8 -- transitions possibles 
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des transitions possibles de chaque branch -- len de sous-list de sous-list = nb de type de transitions de chaque branch
     
     tran_energy_tot -- indice 9 -- énergies de transitions
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des énergies de transitions possibles de chaque branch -- len de sous-list de sous-list = nb de type de transitions de chaque branch
     
     tran_prob_tot -- indice 10 -- probabilités de transitions
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des probabilités de transitions possibles de chaque branch -- len de sous-list de sous-list = nb de type de transitions de chaque branch
     
     tran_level_tot -- indice 11 -- niveaux de branch correspondants
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des niveaux de chaque branch avant des transitions -- len de sous-list de sous-list = 1
     
     tran_level_end_tot -- indice 12 -- niveaux après des transitions
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des niveaux après des transitions de chaque branch -- len de sous-list de sous-list = nb de type de transitions de chaque branch
     
     level_energy_tot -- indice 13 -- énergies de niveaux
         len = nb de noyaux fils
         sous-list -- des branchs possibles de noyau fil -- len de sous-list = nb de branch de chaque fil
         sous-list de sous-list -- des énergies de niveaux de chaque branch -- len de sous-list de sous-list = 1
     '''

     url = "http://www.lnhb.fr/nuclides/"+rad+".PenNuc.txt"
     file = rq.urlopen(url)
     
     # Format the data 
     decayData = []
     for line in file:
        decayData.append(line.decode("utf-8"))
        if "NDA " in decayData[-1]:    decayData[-1] = decayData[-1].replace("NDA ", "NDA; ") # number of daughter
        if "DAU " in decayData[-1]:    decayData[-1] = decayData[-1].replace("DAU ", "DAU; ") # daughter
        if "DDE " in decayData[-1]:    decayData[-1] = decayData[-1].replace("DDE ", "DDE; ") # daughter description : probability of disintegration to NDA, uncertainty, number of excited levels of NDA; number of branche to NDA
        if "Q " in decayData[-1]:      decayData[-1] = decayData[-1].replace("Q ", "Q; ")     # total energy of the branch, uncertainty 
        if "ALP " in decayData[-1]:    decayData[-1] = decayData[-1].replace("ALP ", "ALP; ") # type of the disintegration = alpha
        if "CK " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CK ", "CK; ") # type of the disintegration = Electron Capture K shell
        if "CL " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CL ", "CL; ") # type of the disintegration = Electron Capture L shell
        if "CL1 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CL1 ", "CL1; ") # type of the disintegration = Electron Capture L shell 1
        if "CL2 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CL2 ", "CL2; ") # type of the disintegration = Electron Capture L shell 2
        if "CL3 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CL3 ", "CL3; ") # type of the disintegration = Electron Capture L shell 3
        if "CM " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CM ", "CM; ") # type of the disintegration = Electron Capture M shell
        if "CN " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CN ", "CN; ") # type of the disintegration = Electron Capture N shell
        if "BEM " in decayData[-1]:    decayData[-1] = decayData[-1].replace("BEM ", "BEM; ") # type disintegration-beta- :branching ratio, uncertainty, level fed in daughter, energy, uncertainty,prohibition factor for beta emission
        if "BEP " in decayData[-1]:    decayData[-1] = decayData[-1].replace("BEP ", "BEP; ") # type disintegration-beta+ :branching ratio, uncertainty, level fed in daughter, energy, uncertainty,prohibition factor for beta emission
        if "LED " in decayData[-1]:    decayData[-1] = decayData[-1].replace("LED ", "LED; ") # level description:energy of level,uncertainty, number of transitions that depopulate this level,level time,uncertainty,level number
        if "GA " in decayData[-1]:     decayData[-1] = decayData[-1].replace("GA ", "GA; ")   # type transition: gamma
        if "EK " in decayData[-1]:     decayData[-1] = decayData[-1].replace("EK ", "EK; ")   # type transition: internal conversion -- electron K
        if "EL " in decayData[-1]:    decayData[-1] = decayData[-1].replace("EL ", "EL; ") # type transition: internal conversion -- electron L
        if "EL1 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("EL1 ", "EL1; ") # type transition: internal conversion -- electron L1
        if "EL2 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("EL2 ", "EL2; ") # type transition: internal conversion -- electron L2
        if "EL3 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("EL3 ", "EL3; ") # type transition: internal conversion -- electron L3
        if "EM " in decayData[-1]:     decayData[-1] = decayData[-1].replace("EM ", "EM; ")   # type transition: internal conversion -- electron M
        if "EN " in decayData[-1]:     decayData[-1] = decayData[-1].replace("EN ", "EN; ")   # type transition: internal conversion -- electron N
        if "COM " in decayData[-1]:    decayData[-1] = decayData[-1].replace("COM ", "COM; ") 
        decayData[-1]=decayData[-1].split(";")
        if "\r\n" in decayData[-1][-1]:  decayData[-1][-1] = decayData[-1][-1].replace("\r\n", " ")

     '''
     ========================
     Repérer chaque noyau fil
     ========================

     daughter -- noyau(x) fil(s)
     posi_daug -- l'indice de démarcation de noyau fil
     posi_branch -- l'indice de démarcation de chaque branch
     posi_tran -- l'indice de démarcation de transition
     prob_daug -- probabilité de produire des noyaux fils
     nb_branch -- nombre de branch possible au dessus de l'état fonda (n>0)
     energy_Q -- l'énergie de désintégration

     '''    
     daughter = [];posi_daug = [];prob_daug=[];nb_branch=[];energy_Q=[];
     end = len(decayData)
     for indice,line in enumerate(decayData):
         if "NDA" == line[0]:
             nb_daug = int(line[1])
         if "DAU" == line[0]:
             daughter.append(line[1].replace(" ",""))
         if "COM" == line[0] and "Daughter" in line[1]:
             posi_daug.append(indice)
         if "Q" == line[0]:
             energy_Q.append(float(line[1]))
         if "DDE" == line[0]:
             prob_daug.append(float(line[1]))
             nb_branch.append(int(line[-2]))

     '''
     ==========
     LOOP START
     ==========

     '''  
     posi_end=[]
     desin_type_tot=[];desin_energy_tot=[];desin_prob_tot=[];desin_level_tot=[]
     tran_type_tot=[];tran_energy_tot=[];tran_prob_tot=[];tran_level_end_tot=[]; 
     tran_level_tot=[];level_energy_tot=[]
     prob_branch_tot=[]

     '''
     =============
     LOOP DAUGHTER 
     =============
     '''
     for i1 in range(nb_daug):
        start_p = posi_daug[i1]
        if i1+1 == nb_daug:
            end_p = end
        else:
            end_p = posi_daug[i1+1]

        posi_end_i = []
        for i2 in range(start_p,end_p):
            if "COM" == decayData[i2][0] and "Branch" in decayData[i2][1]:
                posi_end_i.append(i2)
            if "COM" == decayData[i2][0] and "Level" in decayData[i2][1]:
                posi_end_i.append(i2)
        if end_p == end:
            posi_end_i.append(end)
        else:
            posi_end_i.append(posi_daug[i1+1])
        posi_end.append(posi_end_i)
     
        '''
        ==========================
        LOOP Branch and Transition
        ==========================
        '''

        desin_type_daug=[];desin_energy_daug=[];desin_prob_daug=[];desin_level_daug=[];
        tran_type_daug=[];tran_energy_daug=[];tran_prob_daug=[];tran_level_end_daug=[];
        tran_level_daug=[];level_energy_daug=[]
        prob_branch_daug=[]

        for i3 in range(len(posi_end_i)-1):
            start_p1 = posi_end_i[i3]
            end_p1 = posi_end_i[i3+1]
            branch = False
            transition = False
            if "COM" == decayData[start_p1][0] and "Branch" in decayData[start_p1][1]:
                branch=True
            if "COM" == decayData[start_p1][0] and "Level" in decayData[start_p1][1]:
                transition = True

            '''
            ====================================
            LOOP EACH BLOCK OF BRANCH/TRANSITION
            ====================================
            ''' 
            tran_type_b=[];tran_prob_b=[];tran_energy_b=[]; tran_level_end_b=[];
            tran_level_b=[];level_energy_b=[]
            desin_type_b=[];desin_energy_b=[];desin_prob_b=[];desin_level_b=[];
             
            for i4 in decayData[start_p1+1:end_p1]:
                
                if start_p1+1 == end_p1:
                    break
                if branch:
                    if "ALP" == i4[0]:
                        desin_type_b.append("alpha")
                    if "BEP" == i4[0]:
                        desin_type_b.append("beta+")
                    if "BEM" == i4[0]:
                        desin_type_b.append("beta")
                    if "CK" == i4[0]:
                        desin_type_b.append("Atom_K")
                    if "CL" == i4[0]:
                        desin_type_b.append("Atom_L")
                    if "CL1" == i4[0]:
                        desin_type_b.append("Atom_L1")
                    if "CL2" == i4[0]:
                        desin_type_b.append("Atom_L2")
                    if "CL3" == i4[0]:
                        desin_type_b.append("Atom_L3")
                    if "CM" == i4[0]:
                        desin_type_b.append("Atom_M")
                    if "CN" == i4[0]:
                        desin_type_b.append("Atom_N")
                    desin_prob_b.append(float(i4[1]))
                    desin_level_b.append(int(i4[3]))
                    desin_energy_b.append(float(i4[4]))

                if transition:
                    if i4[1] == '  ' or i4[1] == '   ': i4[1] = 0
                    if len(i4)>2 and i4[2] == '  ': i4[2] = 0
                    if len(i4)>4 and i4[4] == '  ': i4[4] = 0
                    if len(i4)>5 and i4[5] == '  ': i4[5] = 0
                    if "LED" == i4[0]:
                        tran_level_b.append(int(i4[-1]))
                        level_energy_b.append(float(i4[1]))
                    if i4[0] == "GA" or i4[0] == "EK" or i4[0] == "EL" or i4[0] == "EL1" or i4[0] == "EL2" or i4[0] == "EL3" or i4[0] == "EM" or i4[0] == "EN":
                        tran_type_b.append(i4[0])
                        tran_prob_b.append(float(i4[1]))
                        tran_energy_b.append(float(i4[3]))
                        tran_level_end_b.append(float(i4[5]))

            if branch:
                desin_type_daug.append(desin_type_b)
                desin_energy_daug.append(desin_energy_b)
                desin_prob_daug.append(desin_prob_b)
                desin_level_daug.append(desin_level_b)
             
            if transition:
                tran_type_daug.append(tran_type_b)
                tran_energy_daug.append(tran_energy_b)
                tran_prob_daug.append(tran_prob_b)
                tran_level_end_daug.append(tran_level_end_b)
                tran_level_daug.append(tran_level_b)
                level_energy_daug.append(level_energy_b)

            if len(desin_prob_b)>0:
                desin_prob_array = np.array(desin_prob_b)
                prob_branch_i = np.sum(desin_prob_array)
                if prob_branch_i >= 1:
                    prob_branch_i = 1
                prob_branch_daug.append(prob_branch_i)
            elif branch and len(desin_prob_b)==0:
                prob_branch_daug.append(0)

        tran_type_daug.append([])
        tran_prob_daug.append([])
        tran_energy_daug.append([])
        tran_level_end_daug.append([])
        tran_level_daug.append([])
        level_energy_daug.append([])

        desin_type_tot.append(desin_type_daug)
        desin_energy_tot.append(desin_energy_daug)
        desin_prob_tot.append(desin_prob_daug)
        desin_level_tot.append(desin_level_daug)

        tran_type_tot.append(tran_type_daug)
        tran_energy_tot.append(tran_energy_daug)
        tran_prob_tot.append(tran_prob_daug)
        tran_level_end_tot.append(tran_level_end_daug)
        tran_level_tot.append(tran_level_daug)
        level_energy_tot.append(level_energy_daug)
        prob_branch_tot.append(prob_branch_daug)

     out = [daughter,prob_daug,energy_Q,desin_type_tot,desin_energy_tot,desin_prob_tot,desin_level_tot,prob_branch_tot,tran_type_tot,tran_energy_tot,tran_prob_tot,tran_level_tot,tran_level_end_tot,level_energy_tot]
     return out

#tic()
#o = readPenNuc1("Ag-110m")
#print(o[-4])
#toc()
#===============================================================================================================
'''
rad = "Am-244m"
out1 = readPenNuc1(rad)
print(rad)
print("    ")
print("DAU",out1[0])
print("    ")
print("prob_daug",out1[1])
print("    ")
print("Q",out1[2])
print("    ")
print("desin_type",out1[3],len(out1[3]))
print("    ")
print("desin-prob",out1[5],len(out1[5]))
print("    ")
print("desin-level",out1[6],len(out1[6]))
print("    ")
print("prob-branch",out1[7])
print("    ")
print("len-tran-type",len(out1[8]))
for i in range(len(out1[8])):
    print("fil",i)
    for i1,p1 in enumerate(out1[8][i]):
        print(p1,len(p1))
    print("energy")
    for i2,p2 in enumerate(out1[9][i]):
        print(p2,len(p2))
    print("tran-prob")
    for i4,p4 in enumerate(out1[10][i]):
        print(p4)
    print("end_level")
    for i3,p3 in enumerate(out1[12][i]):
        print(p3,len(p3))
'''
#print(readPenNuc1('Tc-99m'))
#out = readPenNuc1('Tc-99m')
#print(sampling(out[1]))
#================================== StoppingPower for alpha particle ===========================================

if absolutePath: f_alpha = open('G:\Python_modules\Jialin\Code\Quenching\\alpha_toulene.txt')
else:
#	f_alpha = open('Quenching/alpha_toulene.txt')
    with importlib.resources.path('tdcrpy', 'Quenching') as data_path:
        f_alpha = open(data_path / "alpha_toulene.txt")

data_ASTAR = f_alpha.readlines()
f_alpha.close()
energy_alph = []
dEdx_alph = []
for i in range(np.size(data_ASTAR)):
    data_ASTAR[i] = data_ASTAR[i].split()
    for j in range(2):
        data_ASTAR[i][j] = float(data_ASTAR[i][j])*1e3  # dEdx from MeV.cm2/g to keV.cm2/g; energy from MeV to keV
    energy_alph.append(data_ASTAR[i][0])
    dEdx_alph.append(data_ASTAR[i][1])

def stoppingpowerA(e,rho=0.96,energy_alpha=energy_alph,dEdx_alpha=dEdx_alph):
    """
    Estimation of the stopping power of alpha particles using tabulated values form the ASTAR code
    ref: https://dx.doi.org/10.18434/T4NC7P
    
    Parameters
    ----------
    e : float
        energy of the alpha particle in keV.
    rho : float, optional
        density of the source in g.cm-3. The default is 0.96.
    energy_alpha : list, optional
        the list of energy (in keV) for which the stopping power was calculated with ASTAR. The default is energy_alph.
    dEdx_alpha : list, optional
        the list of stopping powers (in keV.cm2/g) associated with the energy vector. The default is dEdx_alph.

    Returns
    -------
    float
        Interpolated ASTAR estimation of the stopping power.

    """

    energy_alpha = np.array(energy_alpha)
    dEdx_alpha = np.array(dEdx_alpha)
    dEdx = np.interp(e,energy_alpha ,dEdx_alpha)   
    return dEdx*rho                        #unit keV.cm-1

# tic()
# stoppingpowerA(5000,rho=0.96,energy_alpha=energy_alph,dEdx_alpha=dEdx_alph)
# toc() # 0 s

#===============================================================================================

#========================   Nouveau modèle pour calculer le pouvoir d'arrête d'électron ========


if absolutePath: file_TanXia=open('G:\Python_modules\Jialin\Code\Quenching\\TandataUG.txt', "r")
else:
    #file_TanXia=open('Quenching/TandataUG.txt', "r")
    with importlib.resources.path('tdcrpy', 'Quenching') as data_path:
        file_TanXia = open(data_path / "TandataUG.txt")



data_TanXia=file_TanXia.read(); file_TanXia.close()
data_TanXia=data_TanXia.split("\n"); data_TanXia_f = np.empty(len(data_TanXia))
for i, x in enumerate(data_TanXia):
  if i<len(data_TanXia)-1: data_TanXia_f[i]=float(x)

def stoppingpower(e,rho=0.96,Z=5.2,A=11.04,emin=0,file=data_TanXia_f):
    """
    The stopping power of electrons between 20 keV and 1000 keV is a mixture of a radiative loss model [1], and a collision model [2] that has been validated agaisnt the NIST model ESTAR [3] recommanded by the ICRU Report 37 [4].
    At low energy - between 10 eV and 20 keV - the model from Tan and Xia [5] is implemented.
    Refs:
        [1] https://doi.org/10.1016/0020-708x(82)90244-7
        [2] https://www.ijstr.org/final-print/jan2017/Calculations-Of-Stopping-Power-And-Range-Of-Electrons-Interaction-With-Different-Material-And-Human-Body-Parts.pdf
        [3] https://dx.doi.org/10.18434/T4NC7P
        [4] ICRU Report 37, Stopping Powers for Electrons and Positrons
        [5] https://doi.org/10.1016/j.apradiso.2011.08.012
        
    Parameters
    ----------
    e : float
        Energy of the electron in eV.
    rho : float, optional
        density of the source in g.cm-3. The default is 0.96.
    Z : float, optional
        mean charge number of the source. The default is 5.2.
    A : float, optional
        mean mass number of the source. The default is 11.04.
    emin : float, optional
        the minimal energy to consider. The default is 0.
    file : list, optional
        tabulated data form the Tan and Xia model. The default is data_TanXia_f.

    Returns
    -------
    dEdx : float
        Calculated stopping power in MeV.cm-1.

    """
    # e:eV ;rho: g.cm-3
    mc_2 = 0.511 #MeV
    I = 65e-6 #MeV
    NA = 6.02e23
    ahc = 1.437e-13   #MeV.cm
    if e>=20000:
        e1 = e*1e-6 #MeV
        gamma = (e1+mc_2)/mc_2
        gamma_2 = gamma*gamma
        beta = np.sqrt(1-(1/gamma_2))
        beta_2 = beta**2
        tau = e1/mc_2
        terma = np.log(tau**2*(tau+2)/2)
        termb = 1+tau*tau/8-(2*tau+1)*np.log(2)
        termc = (tau+1)**2
        B0 = terma + termb/termc
        sc = 0.1535/beta_2*Z/A*(B0-2*np.log(I/mc_2))  
        term3 = NA*(Z**2)*rho*(e1+mc_2)/(137*(mc_2**2)*A)
        term4 = 4* np.log(2*gamma) -4/3
        sr = (ahc**2)*term3*term4
        #if T<1:sr=0
        dEdx = (sc + sr)*rho  #MeV.cm-1
    else:
            if e>emin:
                dEdx=float(file[int(e)]) #MeV.cm-1
            else:
                dEdx=0
    if dEdx<0:
        dEdx=0
    return dEdx    

# tic()
# stoppingpower(100,rho=0.96,Z=5.2,A=11.04,emin=0,file=data_TanXia_f)
# toc() # 0 s
#===================================================================================================

#===================================================================================================
#============================tracer stopping power of electron and alpha ===========================

'''
e1 = np.linspace(10,1.99e4,10000)
e2 = np.linspace(1.99e4,1e8,20000)
sp1 = []
for i in e1:
    sp1.append(stoppingpower(i))
sp2=[]
for i in e2:
    sp2.append(stoppingpower(i))

ea = np.linspace(1,8e3,4000)
spa = []
for i in ea:
    spa.append(stoppingpowerA(i)*1e-3)
x = ea*1e3
plt.plot(e1,sp1,color='blue',label="pouvoir d'arrête électron")
plt.plot(e2,sp2,color='blue')
plt.plot(x,spa,label="pouvoir d'arrête alpha")
plt.xscale('log')
plt.yscale('log')
plt.title("pouvoir d'arrête d'electron et alpha")
plt.ylabel("pouvoir d'arrête/MeV.cm-1")
plt.xlabel('énergie cinétique/eV')
plt.legend(fontsize=10,loc='best')
plt.savefig('Quenching/stoppingpowerE_A.png')

'''

#=============================================================================================

#====================  Fonction pour lire BetaShape   ========================================
# file_betashape = "decayData//All-nuclides_BetaShape.zip"
with importlib.resources.path('tdcrpy', 'decayData') as data_path:
    file_betashape = data_path / "All-nuclides_BetaShape.zip"

z_betashape = zf.ZipFile(file_betashape)

def readBetaShape(rad,mode,level,z=z_betashape):
    """
    This funcion reads the beta spectra calculated by the code BetaShape and published in the DDEP web page.
    refs:
        https://doi.org/10.1103/PhysRevC.92.059902
        http://www.lnhb.fr/ddep_wg/

    Parameters
    ----------
    rad : string
        identifier of the radionuclide. e.g. 'Na-22'
    mode : string
        identifier of the decay mode. 'beta-' or 'beta+'
    level : int
        level of the daughter after decay.  0,1,2,3 ....
    Returns
    -------
    e : list
        the energy vector in keV.
    dNdx : list
        the probability density in keV-1.

    """

    Rad = rad.replace('-','')
    name_doc = Rad+'/'+mode+'_'+Rad+'_'+ "trans" + str(level) +'.bs'
    with z.open(name_doc) as file_trans:
        data = file_trans.readlines()

    for i in range(np.size(data)):
        data[i] = str(data[i])
        data[i] = data[i].replace("b'",'')
        data[i] = data[i].replace("\\r\\n",'')
        data[i] = data[i].replace("'",'')
    for i in range(np.size(data)):
        data[i] = data[i].split()
    e = []
    dNdx = []
    
    while [] in data:
        data.remove([])
    
    for i in range(len(data)):
        ind = i
        if data[i][0] == 'E(keV)':break
    
    for j in range(i+1,len(data)):
        e.append(float(data[j][0])) # convert to float
        dNdx.append(float(data[j][1])) # convert to float
    dNdx /= sum(np.asarray(dNdx)) # normalization
    dNdx = list(dNdx)
    return e, dNdx

#tic()
#e,p = readBetaShape('C-11','beta+',0)
#toc()
#print(p,type(p))

# tic()
# e,p = readBetaShape('C-11','beta+',0)
# toc() # 0.016 s > 0s

def readTDCR17spectra(rad):
    file = open("decayData/spectra/spectrumTDCR17_"+rad+".dat")
    data = file.read()
    file.close()
    data = data.split("\n")
    e = []; r = []
    for i in data:
        if "keV" not in i:
            a = i.split("  ")
            if len(a)>1:
                e.append(float(a[1])*1e-3)
                r.append(float(a[2]))
    return e, r


# tic()
# e,r =  readBetaShape("Nb-95m", "beta-", 1)
# # print(e[-1])
# toc() # 0.016 s > 0 s

#=====================================================================================

## Display beta spectra
# rplot = "S-35"
# out = readBetaShape(rplot, "beta-", "tot"); print(sum(out[1]))
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

#=======================================================================================

#============================  Fonction quenching  =====================================

def E_quench_e(e,kB,nE):
    """
    This function calculate the quenched energy of electrons according to the Birks model of scintillation quenching

    Parameters
    ----------
    e : float
        energy of the electron in eV.
    kB : float
        Birks constant in cm/MeV.

    Returns
    -------
    float
        Quenched energy in eV.

    """
    
    e_dis = np.linspace(0,e,nE)
    delta = e_dis[2] - e_dis[1]
    q = 0
    for i in e_dis:
        q += delta/(1+kB*stoppingpower(i))
    return q

def E_quench_a(e,kB,nE): 
    """
    This function calculate the quenched energy alpha particles according to  the Birks model of scintillation quenching

    Parameters
    ----------
    e : float
        energy of the alpha particle in keV.
    kB : float
        Birks constant in cm/keV.

    Returns
    -------
    float
        Quenched energy in keV.

    """

    e_dis = np.linspace(1,e,nE)
    delta = e_dis[2] - e_dis[1]
    q = 0
    for i in e_dis:
        q += delta/(1+kB*stoppingpowerA(i))
    return q

#=========================================================================================



#print(E_quench_e(3.996e5,1e-2,1000)*1e-3)

#======================================================================================

#========================= Tracer les courbes avec kB différents ======================

'''

# s1 = []
# s2 = []
# s3 = []
# x = np.linspace(5,8000,4000) 

# for i in x:
#     s1.append(E_quench_a(i,kB=7e-6)/i)
#     s2.append(E_quench_a(i,kB=1e-5)/i)
#     s3.append(E_quench_a(i,kB=1.4e-5)/i)

# plt.plot(x,s1,label='kB=0.007cm/MeV',color='green',lw=2)
# plt.plot(x,s2,label='kB=0.01cm/MeV',ls=':',color='red',lw=3)
# plt.plot(x,s3,label='kB=0.014cm/MeV')
# plt.xscale('log')
# #plt.yscale('log')
# plt.legend(fontsize=12,loc='best')
# plt.xlabel('E de particule/keV')
# plt.ylabel("énergie d'extinction/E")
# plt.savefig("Quenching/beta 100-10K E_Q sur E.png")

#'''
#============================================================================================

#============================================================================================

#========================= énergie gamma ===================================================
#'''
if absolutePath: 
    fp1 = 'G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/matrice_10ml-photon_1_200k.txt'       #gamma-10ml-1-200keV-niveau 0
    fp2 = 'G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/matrice_10ml-photon_200_2000k.txt'    #gamma-10ml-200-2000keV-niveau 1
    fp3 = 'G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/matrice_10ml-photon_2000_10000k.txt'  #gamma-10ml-2000-10000keV-niveau 2
    fp4 = 'G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/matrice_16ml-photon_1_200k.txt'       #gamma-10ml-1-200keV-niveau 0
    fp5 = 'G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/matrice_16ml-photon_200_2000k.txt'       #gamma-10ml-1-200keV-niveau 1
    fp6 = 'G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/matrice_16ml-photon_2000_10000k.txt'       #gamma-10ml-1-200keV-niveau 2
    fe = "G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/E_depose.txt"  
else:
    # fp1 = 'MCNP-MATRIX/matrice/fichier/matrice_10ml-photon_1_200k.txt'      #gamma-10ml-1-200keV-niveau 0
    # fp2 = 'MCNP-MATRIX/matrice/fichier/matrice_10ml-photon_200_2000k.txt'   #gamma-10ml-200-2000keV-niveau 0
    # fp3 = 'MCNP-MATRIX/matrice/fichier/matrice_10ml-photon_2000_10000k.txt' #gamma-10ml-2000-10000keV-niveau 0
    # fp4 = 'MCNP-MATRIX/matrice/fichier/matrice_16ml-photon_1_200k.txt'      #gamma-10ml-1-200keV-niveau 0
    # fp5 = 'MCNP-MATRIX/matrice/fichier/matrice_16ml-photon_200_2000k.txt'      #gamma-10ml-1-200keV-niveau 0
    # fe = "MCNP-MATRIX/matrice/fichier/E_depose.txt"
    with importlib.resources.path('tdcrpy', 'MCNP-MATRIX') as data_path:
        fp1 = data_path / 'matrice/fichier/matrice_10ml-photon_1_200k.txt'      #gamma-10ml-1-200keV-niveau 0
        fp2 = data_path / 'matrice/fichier/matrice_10ml-photon_200_2000k.txt'   #gamma-10ml-200-2000keV-niveau 1
        fp3 = data_path / 'matrice/fichier/matrice_10ml-photon_2000_10000k.txt' #gamma-10ml-2000-10000keV-niveau 2
        fp4 = data_path / 'matrice/fichier/matrice_16ml-photon_1_200k.txt'      #gamma-10ml-1-200keV-niveau 0
        fp5 = data_path / 'matrice/fichier/matrice_16ml-photon_200_2000k.txt'      #gamma-10ml-1-200keV-niveau 1
        fp6 = data_path / 'matrice/fichier/matrice_16ml-photon_2000_10000k.txt'      #gamma-10ml-1-200keV-niveau 2
        fe = data_path / 'matrice/fichier/E_depose.txt'

'''
data1 = f1.readlines()
data2 = f2.readlines()
data3 = f3.readlines()
data_e = fe.readlines()

Matrice1 = np.zeros((1003,200))
Matrice2 = np.zeros((1003,901))
Matrice3 = np.zeros((1003,801))
Matrice_e = np.zeros((1002,3))

for i in range(1002):
    data_e[i] = data_e[i].split()
    for j in range(3):
        Matrice_e[i][j] = float(data_e[i][j])
for i in range(1003):
    data1[i] = data1[i].split()
    data2[i] = data2[i].split()
    data3[i] = data3[i].split()
    for j in range(200):
        Matrice1[i][j] = float(data1[i][j])
    for k in range(901):
        Matrice2[i][k] = float(data2[i][k])
    for l in range(801):
        Matrice3[i][l] = float(data3[i][l])
'''

def read_matrice(path,niveau):
    f = open(path)
    data = f.readlines()
    if niveau == 0:
        taille_x = 200
        taille_y = 1003
    elif niveau == 1:
        taille_x = 901
        taille_y = 1003
    elif niveau == 2:
        taille_x = 801
        taille_y = 1003
    elif niveau=='e':
        taille_x = 3
        taille_y = 1002

    matrice = np.zeros((taille_y,taille_x))
    for i in range(taille_y):
        data[i] = data[i].split()
        for j in range(taille_x):
            matrice[i][j] = float(data[i][j])
    return matrice


Matrice10_p_1 = read_matrice(fp1,0)
Matrice10_p_2 = read_matrice(fp2,1)
Matrice10_p_3 = read_matrice(fp3,2)
Matrice16_p_1 = read_matrice(fp4,0)
Matrice16_p_2 = read_matrice(fp5,1)
Matrice16_p_3 = read_matrice(fp6,2)
Matrice_e = read_matrice(fe,'e')

def energie_dep_gamma(e_inci,v,matrice10_1=Matrice10_p_1,matrice10_2=Matrice10_p_2,matrice10_3=Matrice10_p_3,matrice16_1=Matrice16_p_1,matrice16_2=Matrice16_p_2,matrice16_3=Matrice16_p_3,ed=Matrice_e):
    """ 
    ----------
    Parameters
    ----------
    e_inci : float
        l'énergie incidente de particule.
    * : TYPE
        DESCRIPTION.
    matrice10_1 : matrix
        matrice de photon de 1-200keV de solution 10ml.
    matrice2 : TYPE, optional
        matrice de photon de 200-2000keV de solution 10ml.
    matrice3 : TYPE, optional
        matrice de photon de 2000-10000keV de solution 10ml.
    ed : TYPE, optional
        matrice de bins d'énergie. colonne 0: 1-200keV; colonne 1: 200-2000keV

    Returns
    -------
    result : float
        l'énergie déposée.

    """
    ## sort keV / entrée : keV
    if e_inci <= 200:
        if e_inci < 1:
            index = 0            # index de colonne de la matrice de l'énergie incidente la plus proche 
        else:
            index = int(e_inci)-1
            
        if v == 10: 
            matrice = matrice10_1
        elif v == 16:
            matrice = matrice16_1
        e = ed[:,0]
    
    elif e_inci <= 2000:
        index = int((e_inci-200)/2)
        if v == 10: 
            matrice = matrice10_2
        elif v == 16:
            matrice = matrice16_2
        e = ed[:,1]

    else:
        index = (int(e_inci)-2000)//10
        if v == 10: 
            matrice = matrice10_3
        elif v == 16:
            matrice = matrice16_3
        e = ed[:,2]
    
    inde = sampling(matrice[1:,index])
    if inde == 1 : result = 0
        #elif e_inci<25: result = e[inde-1]*1e3*e_inci/matrice[0][index]
    else: result = e[inde]*1e3*e_inci/matrice[0][index]
    if result  > e_inci: result = e_inci
    return result

#for i in range(10):
 #   print(energie_dep_gamma(511,16))


if absolutePath: 
    fe1 = 'G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/matrice_10ml-beta-_1_200k.txt' # electron-10ml-1-200keV-niveau 0
    fe2 = 'G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/matrice_10ml-beta-_200_2000k.txt' # electron-10ml-200-2000keV-niveau 1
    fe3 = 'G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/matrice_10ml-beta-_2000_10000k.txt' # electron-10ml-2000-10000keV-niveau 2
    fe = "G:\Python_modules\Jialin\Code\\MCNP-MATRIX/matrice/fichier/E_depose.txt"   # electron-10ml-énergie-niveau 'e'
else:
    # fe1 = 'MCNP-MATRIX/matrice/fichier/matrice_beta-_1_200k.txt' # electron-10ml-1-200keV-niveau 0
    # fe2 = 'MCNP-MATRIX/matrice/fichier/matrice_beta-_200_2000k.txt' # electron-10ml-200-2000keV-niveau 1
    # fe3 = 'MCNP-MATRIX/matrice/fichier/matrice_beta-_2000_10000k.txt' # electron-10ml-2000-10000keV-niveau 2
    # fe = "MCNP-MATRIX/matrice/fichier/E_depose.txt" # electron-10ml-énergie-niveau 'e'
    with importlib.resources.path('tdcrpy', 'MCNP-MATRIX') as data_path:
        fe1 = data_path / 'matrice/fichier/matrice_10ml-beta-_1_200k.txt' # electron-10ml-1-200keV-niveau 0
        fe2 = data_path / 'matrice/fichier/matrice_10ml-beta-_200_2000k.txt' # electron-10ml-200-2000keV-niveau 1
        fe3 = data_path / 'matrice/fichier/matrice_10ml-beta-_2000_10000k.txt' # electron-10ml-2000-10000keV-niveau 2
        fe = data_path / 'matrice/fichier/E_depose.txt' # electron-10ml-énergie-niveau 'e'   



Matrice10_e_1 = read_matrice(fe1,0)
Matrice10_e_2 = read_matrice(fe2,1)
Matrice10_e_3 = read_matrice(fe3,2)
#Matrice_e = read_matrice(fe,'e')


def energie_dep_beta(e_inci,*,matrice10_1=Matrice10_e_1,matrice10_2=Matrice10_e_2,matrice10_3=Matrice10_e_3,ed=Matrice_e):
    """ 
    ----------
    Parameters
    ----------
    e_inci : float
        l'énergie incidente de particule.
    * : TYPE
        DESCRIPTION.
    matrice10_1 : matrix
        matrice d'électrons de 1-200keV de solution 10ml.
    matrice2 : TYPE, optional
        matrice d'électrons de 200-2000keV de solution 10ml.
    matrice3 : TYPE, optional
        matrice d'électrons de 2000-10000keV de solution 10ml.
    ed : TYPE, optional
        matrice de bins d'énergie. colonne 0: 1-200keV; colonne 1: 200-2000keV

    Returns
    -------
    result : float
        l'énergie déposée.

    """
    ## sort keV / entrée : keV
    if e_inci <= 200:
        if e_inci < 1:
            index = 0            # index de colonne de la matrice de l'énergie incidente la plus proche 
        else:
            index = int(e_inci)-1
        matrice = matrice10_1
        e = ed[:,0]
    
    elif e_inci <= 2000:
        index = int((e_inci-200)/2)
        #doc = 'MCNP-MATRIX/matrice/matrice_p_200_2000k.txt'
        matrice = matrice10_2
        #taille_x = 901
        e = ed[:,1]

    else:
        index = (int(e_inci)-2000)//10
        #doc = 'MCNP-MATRIX/matrice/matrice_p_2000_10000k.txt'
        matrice = matrice10_3
        #taille_x = 801
        e = ed[:,2]
    
    inde = sampling(matrice[1:,index])
    if inde == 1 : result = 0
        #elif e_inci<25: result = e[inde-1]*1e3*e_inci/matrice[0][index]
    else: result = e[inde]*1e3*e_inci/matrice[0][index]
    if result  > e_inci: result = e_inci
    return result

'''
r = []   
for i in range(10):
    r.append(energie_dep_beta(17.055))
print(r)
'''
#print(Matrice_e[0:5,:])

# tic()
# energie_dep_gamma(800)
# toc()   # 0 s



## Non parametric regression

# def regress(x,y):
#     z = sm.nonparametric.lowess(y, x, return_sorted = True)
#     return z

    

def writeEffcurves(x,y,uy,rad,p,kB,SDT):
    if SDT == "S":
        file = open("EfficiencyCurves/"+''.join(rad)+"/EffS_"+''.join(rad)+'_'+''.join(str(p))+'_'+str(kB)+".txt","w")
    elif SDT == "D":
        file = open("EfficiencyCurves/"+''.join(rad)+"/EffD_"+''.join(rad)+'_'+''.join(str(p))+'_'+str(kB)+".txt","w")
    elif SDT == "T":
        file = open("EfficiencyCurves/"+''.join(rad)+"/EffT_"+''.join(rad)+'_'+''.join(str(p))+'_'+str(kB)+".txt","w")
    else:
        print("Warning: unknown profil type")
    for i, xi in enumerate(x):
        file.write(str(xi)+" "+str(y[i])+" "+str(uy[i])+"\n")
    file.close()

#======================== read ENSDF ============================================
def transf_name(rad):     #  transformer le nom de rad par exemple '11C' à 'C11'
    """
    ---------
    PARAMETRE
    ---------
    rad -- type: str par exemple '108AG'
    
    ------
    RETURN
    ------
    RAD -- type: str par exemple 'AG108' qui correspond à la structure de fille de PenNuc

    """
    # rad -- type : str par exemple '108AG'

    name_lis = re.split('(\d+)',rad)
    RAD = name_lis[2]+name_lis[1]
    return RAD

#print(transf_name('108PD'))

# file_ensdf = 'decayData//All-nuclides_Ensdf.zip'
with importlib.resources.path('tdcrpy', 'decayData') as data_path:
    file_ensdf = data_path / 'All-nuclides_Ensdf.zip'

z_ensdf = zf.ZipFile(file_ensdf)
#print(z.namelist())
def readEShape(rad, *, z=z_ensdf):
    """
    --------------------------------------------------
    pour lire les fichiers dans All-nuclides_Ensdf.zip
    --------------------------------------------------
    ---------
    PARAMETRE
    ---------
    rad -- type: str par exemple 'Ag-108'
    z -- ENSDF files
    
    ------
    RETURN
    ------
    daug_name -- type: list -- les filles de désintégration
    Energy ----- type: list -- chaque élément comprend toutes les énergies de transition de la fille de même indice
    Prob ------- type: list -- chaque élément comprend toutes les proba de transition de la fille de même indice
    Type ------- type: list -- chaque élément comprend touts les types de transition de la fille de même indice

    """
       
    name = rad + '.txt'
    with z.open(name) as f:
        data = f.readlines()
        nl = np.size(data)
        #print(data)
        for i in range(nl):
            data[i] = str(data[i])
            data[i] = data[i].replace("b",'')
            data[i] = data[i].replace("\\r\\n",'')
            data[i] = data[i].replace("'",'')
            if "\\n" in data[i]:data[i] = data[i].replace("\\n","")     # pour traiter "\\n" dans Mn-52 et Mn-52m (cas particuliers)
        #print(data)
        for i in range(nl):
            data[i] = data[i].split()
        
        for i in range(nl):
            if i>0 and ('L' in data[i]) and ("AUGER" in data[i]) and ("|]" in data[i-1]):
                data.insert(i,[data[i][0],'T'])
    #print(data)
    index_auger = []
    index_end = []
    daug_name = []
    posi = []
    for i,p in enumerate(data):
        if 'DECAY' in p:
            daug_name.append(transf_name(p[0]))
        if 'Auger' in p:
            index_auger.append(i)
        if len(p)==2 and 'T' in p:
            posi.append(i)
        if 'P' in p:
            index_end.append(i)
            posi.append(i)

    Energy = []           # enregistrer les résultats (énergie) complètes
    energy = []           # enregistrer les résultats (énergie) d'une fille
    Type = []             # enregistrer les résultats (type de transition) complètes
    type_ = []            # enregistrer les résultats (type de transition) d'une fille
    Prob = []             # enregistrer les résultats (proba de transition) complètes
    prob = []             # enregistrer les résultats (proba de transition) d'une fille

    for i in range(len(posi)-1):
        start = posi[i]+1
        end = posi[i+1]
        d = data[start:end]   # bloc
        e = []                # enregistrer les résultats (énergie) d'un bloc
        prob_b = []           # enregistrer les résultats (proba) d'un bloc
        type_b = []           # enregistrer les résultats (type) d'un bloc
        if start==end:        # sauter les lignes blaches et continues
            continue
        if start-1 in index_end:   # sauter le bloc entre deux filles
            continue

        for n,p1 in enumerate(d):   
            if '-' in p1[2]:                # calculer et remplacer les intervalles
                x = p1[2].split('-') 
                p1[2] = round((float(x[0])+float(x[1]))/2,3)

            if '(total)' in d[0]:           # traiter le bloc qui comprend (total) dans la première ligne 
                if '(total)' in p1:
                    prob_b.append(float(p1[3]))
                    e.append(p1[2])
                    type_b.append(p1[-2])
                continue 
            elif '|]' in p1:                # traiter un bloc qui comprend |]
                if len(p1)>6:               # repérer la ligne qui comprend la proba
                    prob_b.append(float(p1[4]))
                e.append(float(p1[2]))      # enregistrer les valeurs d'énergie
                if 'AUGER' in p1:           # traiter le cas d'auger et |] 
                    if 'K' in p1[-2]:       # Auger K
                        type_b.append('Auger K')
                    else:
                        print('erreur')
                elif 'X' in p1[-1]:         # Rayon X
                    type_b.append(p1[-1][0:3])
            else:                           # traiter le cas sans |] ni (total)
                if len(p1)==4 and 'X' in p1[-1]:   # le cas de rayon X sans |] ni (total) ni proba
                    continue                       # sauter cette ligne
                elif len(p1)==5 and 'L' in p1:     # le cas de Auger L sans |] ni (total) ni proba
                    continue                       # sauter cette ligne
                else:                       # traiter le cas sans |] ni (total) mais complet
                    e.append(float(p1[2]))         # enregistrer énergie
                    prob_b.append(float(p1[3]))    # enregistrer proba
                    if 'L' in p1:
                        type_b.append('Auger L')   # enregistrer type Auger L
                    else:
                        type_b.append(p1[-1][0:3]) # enregistrer type Rayon X

        if len(prob_b)==1 and len(e)>1:            # calculer la valeur moyenne et l'enregistrer au cas où |] compris et valeurs complètes
            energy.append(np.mean(e))
            prob.append(prob_b[0])
            type_.append(type_b[0])
        elif len(e)==len(prob_b) and len(e)>=1:    # enregistrer les valeurs au cas où sans |] et valeurs complètes
            for i in range(len(e)):
                energy.append(e[i])
                prob.append(prob_b[i])
                type_.append(type_b[i])
        if end in index_end or end+1 in index_end: # enregistrer les résultats à la fin d'une fille
            Energy.append(energy)
            Prob.append(prob)
            Type.append(type_)
            energy = []
            prob = []
            type_ = []
    #Energy =  np.array(Energy)
    #Prob =  np.array(Prob)         
    return  daug_name,Energy,Prob,Type        

# tic()
# d,e,p,t = readEShape('Cf-252')
# toc() # 0 s
#========  tester readEShape ==============
# tic()
# d,e,p,t = readEShape('Cf-252')
#print(d,e,p,t)
#d1,e1,p1,t1 = readEShape('At-218')
#print(d1,e1,p1,t1)
# toc()
# print(d,e[0][1],p[1][2],t)
# print('  ')
# for i in range(len(d)):
#     print(d[i],e[i],p[i],t[i])
#     print(' ')
#     print(len(e[i]),len(p[i]),len(t[i]))

#============  traiter la relaxation ===============
def relaxation_atom(daugther,rad,lacune='defaut'):
    """
    ---------
    PARAMETRE
    ---------
    daugther -- type: str -- la fille tirée dans cette itération (par exemple NB95,PD110 etc.)
    rad ------- type: str -- le radionucléide étudié (par exemple Am-241, C-11 etc.) 
    lacuen ---- type: str -- la lacune atomique (par exemple 'Atom_K','Atom_L' etc.)

    ------
    RETURN
    ------
    Type ---- type de transition: Auger L ou K ou Rayon X
    Energy -- énergie correspondante

    """
    daug_name,Energy,Prob,Type = readEShape(rad)  # tirer les vecteurs de rad d'Ensdf 

    index_daug = daug_name.index(daugther)        # repérer l'indice de fille correspondante
    
    Energie = np.array(Energy[index_daug])                  # tirer le vecteur d'énergie
    probability = np.array(Prob[index_daug])                # tirer le vecteur de proba
    type_transi = Type[index_daug]                # tirer le vecteur de type


    if len(probability) > 0:                      # le cas où le vecteur de proba/energie/type n'est pas vide
        '''
        posi_L = []
        posi_K = []
        
        for i, p in enumerate(type_transi):       # repérer les indices de transition L ou K
            if 'L' in p:
                posi_L.append(i)                  # enregistrer la posotion de transition L

            elif 'K' in p:
                posi_K.append(i)                  # enregistrer la posotion de transition L
        '''
        if 'L' in lacune:                         # traiter le transition de couche L
            prob_2 = []
            energy_2 = []
            type_2 = []
            for il, pl in enumerate(type_transi):
                if 'L' in pl:
                    prob_2.append(probability[il])    # enregistrer les proba de transition L
                    energy_2.append(Energie[il])      # enregistrer les energies de transition L
                    type_2.append(type_transi[il])    # enregistrer les types de transition L

        elif 'K' in lacune:                           # traiter le transition de couche K
            prob_2 = []
            energy_2 = []
            type_2 = []
            for ik, pk in enumerate(type_transi):     
                if 'K' in pk:
                    prob_2.append(probability[ik])    # enregistrer les proba de transition K
                    energy_2.append(Energie[ik])      # enregistrer les energie de transition K
                    type_2.append(type_transi[ik])    # enregistrer les type de transition K

        #elif lacune=='defaut':                        # traiter le cas particulier qui ne précise pas la lacune
            #prob_2 = probability
            #energy_2 = Energie
            #type_2 = type_transi
            
        else: # try to debugg
            # print("issue: ", lacune)
            prob_2 = 0   #probability
            energy_2 = 0  # Energie
            #if "M" in lacune:
                #type_2 = "Atom_M"  #type_transi
            #if "N" in lacune:
                #type_2 = "Atom_N"            
        
     # sampling
        #if len(probability)>1:               # le cas où la taille du vecteur de proba supérieur à 1

            #prob_somme = np.sum(prob_2)      # calculer la somme de proba
            #prob_2 /= prob_somme             # normaliser la proba
        if lacune != "Atom_M" and lacune != "Atom_N":
            if len(prob_2) != 0:
                prob_2 = np.array(prob_2)           # convert to array
                if len(probability)>1:
                    prob_somme = np.sum(prob_2)      # calculer la somme de proba
                    prob_2 /= prob_somme 
                index_fin = sampling(prob_2)        # sample in probability of transition
                type_fin = type_2[index_fin]        # type of transition     
                energie_fin = energy_2[index_fin]   # energy of the transition
            else:
                # print("pas de transition de rayon X ni d'électron Auger pour cette lacune: ",lacune)
                type_fin = 'NON'
                energie_fin = 0            
        else:
            type_fin = 'NON'
            energie_fin = 0 
    else:                                            # le cas où le vecteur de proba est vide 
        #print("pas de transition de rayon X ni d'électron Auger")
        type_fin = 'NON'
        energie_fin = 0
    return type_fin,energie_fin
# tic()
# tf,ef = relaxation_atom('CR52', 'Mn-52', 'Atom_K')
# toc() # 0 s
#print(tf,ef)
