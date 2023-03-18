# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:04:46 2023

Library of function for the TDCR model

@author: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

### IMPORT Python Module
import urllib.request as rq
import numpy as np
from numpy.core.multiarray import where
import zipfile as zf

def sampling(p_x):
    """
    This function aims to sample in a pdf or a pmf

    Parameters
    ----------
    p_x : float vector
        Probability Density (or mass) Function (PDF or PMF) of the random variable x.
    x : vector of any kind (Optional)
        values of the random varible X.

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
    This function reads the PenNuc file 
    """
    
    url = "http://www.lnhb.fr/nuclides/"+rad+".PenNuc.txt"
    file = rq.urlopen(url)
    
    # Format the data 
    decayData = []
    for line in file:
      decayData.append(line.decode("utf-8"))
      if "NDA " in decayData[-1]:    decayData[-1] = decayData[-1].replace("NDA ", "NDA; ") # number of daughter
      if "DDE " in decayData[-1]:    decayData[-1] = decayData[-1].replace("DDE ", "DDE; ") # daughter description : probability of disintegration to NDA, uncertainty, number of excited levels of NDA; number of branche to NDA
      if "Q " in decayData[-1]:      decayData[-1] = decayData[-1].replace("Q ", "Q; ")     # total energy of the branch, uncertainty 
      if "ALP " in decayData[-1]:    decayData[-1] = decayData[-1].replace("ALP ", "ALP; ") # type of the disintegration
      if "BEM " in decayData[-1]:    decayData[-1] = decayData[-1].replace("BEM ", "BEM; ") # type disintegration-beta- :branching ratio, uncertainty, level fed in daughter, energy, uncertainty,prohibition factor for beta emission
      if "LED " in decayData[-1]:    decayData[-1] = decayData[-1].replace("LED ", "LED; ") # level description:energy of level,uncertainty, number of transitions that depopulate this level,level time,uncertainty,level number
      if "GA " in decayData[-1]:     decayData[-1] = decayData[-1].replace("GA ", "GA; ")   # type transition: gamma
      if "EK " in decayData[-1]:     decayData[-1] = decayData[-1].replace("EK ", "EK; ")   # type transition: internal conversion -- electron K
      if "EL1 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("EL1 ", "EL1; ") # type transition: internal conversion -- electron L1
      if "EL2 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("EL2 ", "EL2; ") # type transition: internal conversion -- electron L2
      if "EL3 " in decayData[-1]:    decayData[-1] = decayData[-1].replace("EL3 ", "EL3; ") # type transition: internal conversion -- electron L3
      if "EM " in decayData[-1]:     decayData[-1] = decayData[-1].replace("EM ", "EM; ")   # type transition: internal conversion -- electron M
      if "EN " in decayData[-1]:     decayData[-1] = decayData[-1].replace("EN ", "EN; ")   # type transition: internal conversion -- electron N
      decayData[-1]=decayData[-1].split(";")
      if "\r\n" in decayData[-1][-1]:  decayData[-1][-1] = decayData[-1][-1].replace("\r\n", " ")
    
    # Probablity vector of the decay branch
    particle=[]; LevelDaughter = []; p_branch = []; u_p_branch = []; e_branch = []; u_e_branch = []
    for d in decayData:
      
      """
      LOOP IN NDAUGH (Daughters)
      """
      
      if d[0] == "NDA": Ndaughter = int(d[1])   # Read the number of daughter
      if d[0] == "DDE": Pdaughter = float(d[1]) # Read the probability of disintegration to daughter
      if d[0] == "Q":   Q_value = float(d[1])   # Read the Q-value
    
      """
      LOOP IN NBRANCH (Branches of each daughter)
      """
    
      if d[0] == "ALP" or d[0] == "BEM":                         # Read information on the decay branch
        if d[0] == "ALP": particle.append("alpha")
        if d[0] == "BEM": particle.append("beta")
        if d[2] == "  ": d[2]=0
        p_branch.append(float(d[1])); u_p_branch.append(float(d[2])); # Branching ratio of the decay branch
        e_branch.append(float(d[4])); u_e_branch.append(float(d[5])); # Energy of the decay branch (kinetic energy of the particle)
        LevelDaughter.append(int(d[3]))                               # Level fed in daughter
    
    """
    LOOP IN NLEVEL (Levels for each daughter, starting in NLEVEL, ending in 1)
    """
    levelEnergy = []; fromBranch = []; lineL=[]; listTran = []; levelNumber = []
    transitionType = []; prob = []; u_prob = []; e_trans = []; u_e_trans = []; next_level = []
    transitionType_i = []; prob_i = []; u_prob_i = []; e_trans_i = []; u_e_trans_i = []; next_level_i = []
    for d in decayData:
      if d[0] == "LED":
        # record transition details of the previous level
        transitionType.append(transitionType_i); prob.append(prob_i); u_prob.append(u_prob_i); e_trans.append(e_trans_i); u_e_trans.append(u_e_trans_i); next_level.append(next_level_i)
    #    if int(d[6]) in LevelDaughter :  # Read the information on the possible energy levels after the emission of the alpha particle 
        levelEnergy.append(float(d[1]));     # Energie (rounded) of the level
        listTran.append(int(d[3]))           # Number of transitions that depopulate this level
        levelNumber.append(int(d[6]))        # Level number
        transitionType_i = []; prob_i = []; u_prob_i = []; e_trans_i = []; u_e_trans_i = []; next_level_i = []
          
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
    
    transitionType.append(transitionType_i); prob.append(prob_i); u_prob.append(u_prob_i); e_trans.append(e_trans_i); u_e_trans.append(u_e_trans_i); next_level.append(next_level_i)
    return particle, p_branch, e_branch, LevelDaughter, levelNumber, prob, levelEnergy, transitionType, e_trans, next_level, Q_value



def stoppingpowerA(e,doc,rho): 
    # doc-data of ASTAR(.txt)(unit:MeV)
    # rho: density of the absorber (g.cm-3)
    f = open(doc)
    data = f.readlines()
    energy = []
    dEdx_data = []
    for i in range(np.size(data)):
        data[i] = data[i].split()
        for j in range(2):
            data[i][j] = float(data[i][j])*1e3    #  unit from MeV to keV 
        energy.append(data[i][0])
        dEdx_data.append(data[i][1])
    energy = np.array(energy)  
    dEdx_data = np.array(dEdx_data)        # unit:keV.cm^2.g-1
    dEdx = np.interp(e,energy,dEdx_data)   
    return dEdx*rho                        #unit keV.cm-1


file_TanXia=open('TandataUG.txt', "r"); data_TanXia=file_TanXia.read(); file_TanXia.close()
data_TanXia=data_TanXia.split("\n"); data_TanXia_f = np.empty(len(data_TanXia))
for i, x in enumerate(data_TanXia):
  if i<len(data_TanXia)-1: data_TanXia_f[i]=float(x)


def stoppingpowerE(e,*,za=0.5459,rho,I=60,spmodel="Beth_Tan_Xia",emin=0,file=data_TanXia_f): # valid
    c1=0.57
    c2=1.16
    Cte=0.307075*c1           # Relativist constant /(eV.cm^2) 
    mec2=511000               # Mass energy of an electron /eV
    dEdx=0                    # initial value of the stopping power
   
    if spmodel=="Beth":
        beta=np.sqrt(c2*e/mec2)            # relative velocity
        if e>54:
            dEdx=Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)    
        else:
            dEdx=0
    
    if spmodel=="Beth_extr_1":
        if e>400:
            beta=np.sqrt(c2*e/mec2)            # relative velocity
            if beta>=1:
                beta
            dEdx=Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)    
        else:
            if e>emin:
                beta=np.sqrt(c2*400/mec2)       # relative velocity
                dEdx=20*Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)*e**-0.5
            else:
                dEdx=0

    if spmodel=="Beth_extr_2":
        if e>100:
            beta=np.sqrt(c2*e/mec2)            # relative velocity
            dEdx=Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)    
        else:
            if e>emin:
                beta=np.sqrt(c2*100/mec2)       # relative velocity
                dEdx=0.01*Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)*e
            else:
                dEdx=0
      
    if spmodel=="Beth_extr_3":
        if e>1000:
            beta=np.sqrt(c2*e/mec2)            # relative velocity
            dEdx=Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)
        else:
            if e>emin:
                beta=np.sqrt(c2*1000/mec2)     # relative velocity
                dEdx=1/(1000**-1.1)*Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)*e**-1.1
            else:
                dEdx=0
    
    if spmodel=="Beth_Tan_Xia":
        if e>=20000:
            beta=np.sqrt(c2*e/mec2)            # relative velocity
            if beta>=1:
                beta=0.97
            dEdx=Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)
        else:
            if e>emin:
                dEdx=float(file[int(e)])
            else:
                dEdx=0
    if dEdx<0:
        dEdx=0
    return dEdx

def readBetaSpectrum(rad):
    f=open("spectrum_"+rad+".txt", "r") # open the file 
    ne = 1000   # number of energy bins
    e=np.empty(ne) # declare an empty array to store the energy values
    p=np.empty(ne) # declare an empty array to probability values
    i=0
    for line in f: # scan the line of the file
        e[i]=float(line[:9]) # read the energy value and convert in float 
        p[i]=float(line[9:]) # read the probability value and convert in float
        i+=1
    sumP=sum(p) # sum of probabilities
    p=p/sumP # re-normilize to make the sum of probabilities equal to 1
    f.close() #close the file
    return e, p, ne

def readBetaShape(rad):
    url = "http://www.lnhb.fr/nuclides/"+rad+".BetaShape.zip"
    file = rq.urlopen(url)
    z = zf.ZipFile(file)
    z_content = z.namelist()
    for index,p in enumerate(z_content):
        if "trans0" in p: i_trans0 = index
        if "trans1" in p: i_trans1 = index
    with z.open(z_content[i_trans0]) as file_tran0:
        data = file_trans.readlines()
    for i in range(np.size(data)):
        data[i] = str(data[i])
        data[i] = data[i].split("b'")

    del_term = data[1]
    e = []
    dNdx = []

    for i in range(np.size(data)):
        data[i] = data[i].split(del_term)
    for ind, p in enumerate(data):
        if p == ['E(keV)', 'dN/dE', 'calc.', 'unc.'] : break
    for j in range(ind,np.size(data)):
        e.append(data[j][0])
        dNdx.append(data[j][1])
    return e,dNdx
    