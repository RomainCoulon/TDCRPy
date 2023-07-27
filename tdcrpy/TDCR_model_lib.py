# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 16:04:46 2023

Library of function of the TDCRpy code

@author: Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

"""
======= Import Python Module =======
"""

import importlib.resources
import pkg_resources
import configparser
import numpy as np
import zipfile as zf
import time
import re
import os

"""
======= Import ressource data =======
"""

# import advanced configuration data
config = configparser.ConfigParser()
with importlib.resources.path('tdcrpy', 'config.toml') as data_path:
    file_conf = data_path       
config.read(file_conf)
RHO=config["Inputs"].getfloat("density")
Z=config["Inputs"].getfloat("Z")
A=config["Inputs"].getfloat("A")

# import PenNuc data
with importlib.resources.path('tdcrpy', 'decayData') as data_path:
    file_pennuc = data_path / "All-nuclides_PenNuc.zip"
z_PenNuc = zf.ZipFile(file_pennuc)

# import BetaShape data
with importlib.resources.path('tdcrpy', 'decayData') as data_path:
    file_betashape = data_path / "All-nuclides_BetaShape.zip"
z_betashape = zf.ZipFile(file_betashape)

# import ENSDF data
with importlib.resources.path('tdcrpy', 'decayData') as data_path:
    file_ensdf = data_path / 'All-nuclides_Ensdf.zip'
z_ensdf = zf.ZipFile(file_ensdf)

# import photon interaction data (MCNP6 calculation) 
with importlib.resources.path('tdcrpy', 'MCNP-MATRIX') as data_path:
    fp1 = data_path / 'matrice/fichier/matrice_10ml-photon_1_200k.txt'          #gamma-10ml-1-200keV-niveau 0
    fp2 = data_path / 'matrice/fichier/matrice_10ml-photon_200_2000k.txt'       #gamma-10ml-200-2000keV-niveau 1
    fp3 = data_path / 'matrice/fichier/matrice_10ml-photon_2000_10000k.txt'     #gamma-10ml-2000-10000keV-niveau 2
    fp4 = data_path / 'matrice/fichier/matrice_16ml-photon_1_200k.txt'          #gamma-10ml-1-200keV-niveau 0
    fp5 = data_path / 'matrice/fichier/matrice_16ml-photon_200_2000k.txt'       #gamma-10ml-1-200keV-niveau 1
    fp6 = data_path / 'matrice/fichier/matrice_16ml-photon_2000_10000k.txt'     #gamma-10ml-1-200keV-niveau 2
    fe = data_path / 'matrice/fichier/E_depose.txt'

# import electron interaction data (MCNP6 calculation) 
with importlib.resources.path('tdcrpy', 'MCNP-MATRIX') as data_path:
    fe1 = data_path / 'matrice/fichier/matrice_10ml-beta-_1_200k.txt' # electron-10ml-1-200keV-niveau 0
    fe2 = data_path / 'matrice/fichier/matrice_10ml-beta-_200_2000k.txt' # electron-10ml-200-2000keV-niveau 1
    fe3 = data_path / 'matrice/fichier/matrice_10ml-beta-_2000_10000k.txt' # electron-10ml-2000-10000keV-niveau 2
    fe4 = data_path / 'matrice/fichier/matrice_16ml-beta-_1_200k.txt' # electron-16ml-1-200keV-niveau 0
    fe = data_path / 'matrice/fichier/E_depose.txt' # electron-10ml-énergie-niveau 'e'   

# import stopping power data for electron
with importlib.resources.path('tdcrpy', 'Quenching') as data_path:
    file_TanXia = open(data_path / "TandataUG.txt")

data_TanXia=file_TanXia.read(); file_TanXia.close()
data_TanXia=data_TanXia.split("\n"); data_TanXia_f = np.empty(len(data_TanXia))
for i, x in enumerate(data_TanXia):
  if i<len(data_TanXia)-1: data_TanXia_f[i]=float(x)

# import stopping power data for electron for alpha particle (ASTAR data)
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

"""
======= Library of functions =======
"""

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
    """
    This function is used to ensure that the sum of probability is equal to 1.

    Parameters
    ----------
    p_x : list
        vector of probabilities.

    Returns
    -------
    p : list
        normalized probability vector.

    """
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

def readPenNuc2(rad,z1=z_PenNuc):
    '''
    This function is used to read PenNuc files to format the decay data in lists readable by TDCRPy.
    
    
    Parameters
    ----------
    rad : string
        name of the radionculide (for example: "Am-241"). 

    Returns
    -------
    daughter : list
        list of the daughter nucleus -- indice 0.
    prob_daug : list
        list of probabilities to produce daugter nuclei -- indice 1.
    energy_Q : list
        list of Q value for each transition to a given daughter nucleus -- indice 2.
    desin_type_tot : list[list]
        list of type of decay branch / emitted particules -- indice 3. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to possible decay mode of each branch.
    desin_energy_tot : list[list]
        list of the energies of decay transition or the emitted particles -- indice 4. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to possible decay mode of each branch.
    desin_prob_tot : list[list]
        list of the prabability of decay transition or the emitted particles -- indice 5. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to possible decay mode of each branch.
    desin_level_tot : list[list]
        list of energy level that the daughter nucleus can have just after the decay of the mother nucleus -- indice 6. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to possible decay mode of each branch.
    prob_branch_tot : list
        list of branch probabilities -- indice 7. It contains a sub-list for all possible branches of a given daughter nucleus.
    tran_type_tot : list[list]
        list of all possible transitions -- indice 8. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to possible decay mode of each branch.
    tran_energy_tot : list[list]
        list of energy associated with transitions -- indice 9. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to possible decay mode of each branch.
    tran_prob_tot : list[list]
        list of probability associated with transitions -- indice 10. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to possible decay mode of each branch.
    tran_level_tot : list[list]
        list of corresponding branch levels -- indice 11. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to the level before the transition.
    tran_level_end_tot : list[list]
        list of level following given transitions -- indice 12. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to the level after the transition.
    level_energy_tot : list[list]
        list of energy levels -- indice 13. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to possible decay mode of each branch.
    prob_tran_tot : list[list]
        list of sum of transition of each branches -- indice 14. It contains a sub-list for all possible branches of a given daughter nucleus and a sub-sub list related to possible decay mode of each branch.

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

#================================== StoppingPower for alpha particle ===========================================

def stoppingpowerA(e,rho=RHO,energy_alpha=energy_alph,dEdx_alpha=dEdx_alph):
    """
    Estimation of the stopping power of alpha particles using tabulated values form the ASTAR code
    
    ref:
     
        https://dx.doi.org/10.18434/T4NC7P
    
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


#===============================================================================================

#========================   Nouveau modèle pour calculer le pouvoir d'arrête d'électron ========

def stoppingpower(e,rho=RHO,Z=Z,A=A,emin=0,file=data_TanXia_f):
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

#=============================================================================================

#====================  Fonction pour lire BetaShape   ========================================

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
    level : int or string
        level of the daughter after decay.  0,1,2,3 .... or 'tot' in case of pure beta emitting radionuclides
    Returns
    -------
    e : list
        the energy vector in keV.
    dNdx : list
        the probability density in keV-1.
    
    """

    Rad = rad.replace('-','')
    if level == 'tot':
        name_doc = Rad+'/'+mode+'_'+Rad+'_tot.bs'
    else:
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

#============================================================================================

#============================================================================================

#========================= énergie gamma ===================================================
#'''

def read_matrice(path,niveau):
    """
    This function read the response matrix calculated by MCNP6 simulation.

    Parameters
    ----------
    path : string
        path to the response matrix.
    niveau : integer or string
        energy range of the response matrix. 0: [1-200] keV; 1: [200-2000] keV; 2: [2000-10000] keV. "e" for the input energy matrix.

    Returns
    -------
    matrice : list[list]
        formatted response matrix.

    """
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
    """ This function samples the energy deposited by a x or gamma rays in the scintillator using response calculated by the Monte-Carlo code MCNP6. 
    
    Parameters
    ----------
    e_inci : float
        energy of the photon in keV.
    v : float
        volume of the scintillator in ml.
    matrice10_1 : list[list], optional
        response matrix for photons in the range [1-200] keV and for a scintillator volume of 10 ml.
    matrice10_2 : list[list], optional
        response matrix for photons in the range [200-2000] keV and for a scintillator volume of 10 ml.
    matrice10_3 : list[list], optional
        response matrix for photons in the range [2000-10000] keV and for a scintillator volume of 10 ml.
    ed : list[list], optional
        matrix of input energies. column 0: [1-200] keV; column 1: [200-2000] keV; column 2: [2000-10000] keV

    Returns
    -------
    result : float
        deposited energy in keV.

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

def energie_dep_gamma2(e_inci,v,matrice10_1=Matrice10_p_1,matrice10_2=Matrice10_p_2,matrice10_3=Matrice10_p_3,matrice16_1=Matrice16_p_1,matrice16_2=Matrice16_p_2,matrice16_3=Matrice16_p_3,ed=Matrice_e):
    """ This function samples the energy deposited by a x or gamma rays in the scintillator using response calculated by the Monte-Carlo code MCNP6. 
    
    Parameters
    ----------
    e_inci : float
        energy of the photon in keV.
    v : float
        volume of the scintillator in ml.
    matrice10_1 : list[list], optional
        response matrix for photons in the range [1-200] keV and for a scintillator volume of 10 ml.
    matrice10_2 : list[list], optional
        response matrix for photons in the range [200-2000] keV and for a scintillator volume of 10 ml.
    matrice10_3 : list[list], optional
        response matrix for photons in the range [2000-10000] keV and for a scintillator volume of 10 ml.
    matrice16_1 : list[list], optional
        response matrix for photons in the range [1-200] keV and for a scintillator volume of 16 ml.
    matrice16_2 : list[list], optional
        response matrix for photons in the range [200-2000] keV and for a scintillator volume of 16 ml.
    matrice16_3 : list[list], optional
        response matrix for photons in the range [2000-10000] keV and for a scintillator volume of 16 ml.
    ed : list[list], optional
        matrix of input energies. column 0: [1-200] keV; column 1: [200-2000] keV; column 2: [2000-10000] keV

    Returns
    -------
    result : float
        deposited energy in keV.

    """
    
    ## sort keV / entrée : keV
    if e_inci <= 200:
        if e_inci < 1:
            index = 0            # index de colonne de la matrice de l'énergie incidente la plus proche 
        else:
            index = int(e_inci)-1
            
        if v == 10: 
            matrice = matrice10_1[1:,index]
            matrice0 = matrice10_1[0,index]
        elif v == 16:
            matrice = matrice16_1[1:,index]
            matrice0 = matrice16_1[0,index]
        else:
            matrice = (matrice16_1[1:,index]-matrice10_1[1:,index])*v/6 + (matrice10_1[1:,index]-(matrice16_1[1:,index]-matrice10_1[1:,index])*10/6)
            matrice0 = (matrice16_1[0,index]-matrice10_1[0,index])*v/6 + (matrice10_1[0,index]-(matrice16_1[0,index]-matrice10_1[0,index])*10/6)
        e = ed[:,0]
    
    elif e_inci <= 2000:
        index = int((e_inci-200)/2)
        if v == 10: 
            matrice = matrice10_2[1:,index]
            matrice0 = matrice10_2[0,index]
        elif v == 16:
            matrice = matrice16_2[1:,index]
            matrice0 = matrice16_2[0,index]
        else:
            matrice = (matrice16_2[1:,index]-matrice10_2[1:,index])*v/6 + (matrice10_2[1:,index]-(matrice16_2[1:,index]-matrice10_2[1:,index])*10/6) 
            matrice0 = (matrice16_2[0,index]-matrice10_2[0,index])*v/6 + (matrice10_2[0,index]-(matrice16_2[0,index]-matrice10_2[0,index])*10/6) 
        e = ed[:,1]

    else:
        index = (int(e_inci)-2000)//10
        if v == 10: 
            matrice = matrice10_3[1:,index]
            matrice0 = matrice10_3[0,index]
        elif v == 16:
            matrice = matrice16_3[1:,index]
            matrice0 = matrice16_3[0,index]
        else:
            matrice = (matrice16_3[1:,index]-matrice10_3[1:,index])*v/6 + (matrice10_3[1:,index]-(matrice16_3[1:,index]-matrice10_3[1:,index])*10/6) 
            matrice0 = (matrice16_3[0,index]-matrice10_3[0,index])*v/6 + (matrice10_3[0,index]-(matrice16_3[0,index]-matrice10_3[0,index])*10/6)
        e = ed[:,2]
    
    inde = sampling(matrice)
    if inde == 1 : result = 0
        #elif e_inci<25: result = e[inde-1]*1e3*e_inci/matrice[0][index]
    else: result = e[inde]*1e3*e_inci/matrice0
    if result  > e_inci: result = e_inci
    return result

Matrice10_e_1 = read_matrice(fe1,0)
Matrice10_e_2 = read_matrice(fe2,1)
Matrice10_e_3 = read_matrice(fe3,2)
Matrice16_e_3 = read_matrice(fe4,0)
#Matrice_e = read_matrice(fe,'e')

def energie_dep_beta(e_inci,*,matrice10_1=Matrice10_e_1,matrice10_2=Matrice10_e_2,matrice10_3=Matrice10_e_3,ed=Matrice_e):
    """ This function samples the energy deposited by an electron in the scintillator using response calculated by the Monte-Carlo code MCNP6. 
    
    Parameters
    ----------
    e_inci : float
        energy of the electron in keV.
    matrice10_1 : list[list], optional
        response matrix for electrons in the range [1-200] keV and for a scintillator volume of 10 ml.
    matrice10_2 : list[list], optional
        response matrix for electrons in the range [200-2000] keV and for a scintillator volume of 10 ml.
    matrice10_3 : list[list], optional
        response matrix for electrons in the range [2000-10000] keV and for a scintillator volume of 10 ml.
    ed : list[list], optional
        matrix of input energies. column 0: [1-200] keV; column 1: [200-2000] keV; column 2: [2000-10000] keV

    Returns
    -------
    result : float
        deposited energy in keV.

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


def energie_dep_beta2(e_inci,*,matrice10_1=Matrice10_e_1,matrice10_2=Matrice10_e_2,matrice10_3=Matrice10_e_3,ed=Matrice_e):
    """ This function samples the energy deposited by an electron in the scintillator using response calculated by the Monte-Carlo code MCNP6. 
    
    Parameters
    ----------
    e_inci : float
        energy of the electron in keV.
    matrice10_1 : list[list], optional
        response matrix for electrons in the range [1-200] keV and for a scintillator volume of 10 ml.
    matrice10_2 : list[list], optional
        response matrix for electrons in the range [200-2000] keV and for a scintillator volume of 10 ml.
    matrice10_3 : list[list], optional
        response matrix for electrons in the range [2000-10000] keV and for a scintillator volume of 10 ml.
    ed : list[list], optional
        matrix of input energies. column 0: [1-200] keV; column 1: [200-2000] keV; column 2: [2000-10000] keV

    Returns
    -------
    result : float
        deposited energy in keV.

    """
    ## sort keV / entrée : keV
    if e_inci <= 200:
        if e_inci < 1:
            index = 0            # index de colonne de la matrice de l'énergie incidente la plus proche 
        else:
            index = int(e_inci)-1
        matrice = matrice10_1[1:,index]
        matrice0 = matrice10_1[0,index]
        e = ed[:,0]
    
    elif e_inci <= 2000:
        index = int((e_inci-200)/2)
        #doc = 'MCNP-MATRIX/matrice/matrice_p_200_2000k.txt'
        matrice = matrice10_2[1:,index]
        matrice0 = matrice10_2[0,index]
        #taille_x = 901
        e = ed[:,1]

    else:
        index = (int(e_inci)-2000)//10
        #doc = 'MCNP-MATRIX/matrice/matrice_p_2000_10000k.txt'
        matrice = matrice10_3[1:,index]
        matrice0 = matrice10_3[0,index]
        #taille_x = 801
        e = ed[:,2]
    
    inde = sampling(matrice)
    if inde == 1 : result = 0
        #elif e_inci<25: result = e[inde-1]*1e3*e_inci/matrice[0][index]
    else: result = e[inde]*1e3*e_inci/matrice0
    if result  > e_inci: result = e_inci
    return result


def writeEffcurves(x,y,uy,rad,p,kB,SDT):
    """
    This function writes efficiency curves

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    uy : TYPE
        DESCRIPTION.
    rad : TYPE
        DESCRIPTION.
    p : TYPE
        DESCRIPTION.
    kB : TYPE
        DESCRIPTION.
    SDT : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
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
def transf_name(rad):
    """ This function format the name of the nuclide to match with the PenNuc format.
    
    Parameters 
    ----------
    rad : string
        name of the radionculdie such as '108AG'.
    
    Returns 
    -------
    RAD : string
        name of the radionuclide such as 'AG108' that match with PenNuc format.

    """
    name_lis = re.split('(\d+)',rad)
    RAD = name_lis[2]+name_lis[1]
    return RAD

def readEShape(rad, *, z=z_ensdf):
    """ This function reads the ENSDF zip files and format the data to be processed by TDCRPy. 

    Parameters
    ----------
    rad : string
        name of the radionuclide such as 'Ag-108'.
    z : ZipFile object
        zip ENSDF file.
    
    Returns
    -------
    daug_name : list
        daughter nucleus of the decay
    Energy : list
        comprise all transition energies of the daughter nucleus.
    Prob : list
        comprise all transtion probabilities of the daughter nucleus.
    Type : list
        comprise all type of transition of the daughter nucleus.

    """
       
    name = rad + '.txt'
    with z.open(name) as f:
        data = f.readlines()
        nl = np.size(data)
        for i in range(nl):
            data[i] = str(data[i])
            data[i] = data[i].replace("b",'')
            data[i] = data[i].replace("\\r\\n",'')
            data[i] = data[i].replace("'",'')
            if "\\n" in data[i]:data[i] = data[i].replace("\\n","")     # pour traiter "\\n" dans Mn-52 et Mn-52m (cas particuliers)
        for i in range(nl):
            data[i] = data[i].split()
        for i in range(nl):
            if i>0 and ('L' in data[i]) and ("AUGER" in data[i]) and ("|]" in data[i-1]):
                data.insert(i,[data[i][0],'T'])
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
    return  daug_name,Energy,Prob,Type        

#============  traiter la relaxation ===============
def relaxation_atom(daugther,rad,lacune='defaut'):
    """ This function simulates the atomic rearangement following a missing electron an inner shell of the daughter atom.
    
    Parameters
    ----------
    daugther : string
        The daughter nucleus (for example NB95,PD110 etc.)
    rad : string
        The mother nucleus (for exemple Am-241, C-11 etc.) 
    lacune  : string
        The shell where the electron is missing (for example 'Atom_K','Atom_L' etc.)

    Returns 
    -------
    Type : type of transition Auger L or K, or X Ray.
    Energy : corresponding energy in keV.

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


def modelAnalytical(L,TD,TAB,TBC,TAC,rad,kB,V,mode,mode2,ne):
    """
    TDCR analytical model that is used for pure beta emitting radionuclides
    
    Parameters
    ----------
    L : float or tuple
        free parameter(s).
    TD : float
        triple-to-double coincidence ratio that was measured (logic sum).
    TAB : float
        triple-to-double coincidence ratio that was measured (channels A and B).
    TBC : flat
        triple-to-double coincidence ratio that was measured (channels B and C).
    TAC : float
        triple-to-double coincidence ratio that was measured (channels A and C).
    rad : string
        radionuclide (eg. "Na-22").
    kB : float
        Birks constant in cm/keV.
    V : float
        volume of the scintillator in ml. run only for 10 ml
    mode : string
        "res" to return the residual, "eff" to return efficiencies.
    mode2 : string
        "sym" for symetrical model, "asym" for symetrical model.
    nE : integer
         Number of bins for the quenching function.
    
    
    Returns
    -------
    res : float
        Residuals of the model compared the measurement data for (a) given free parmeters L. (only in mode="res")
    mean_efficiency_S : float
        Estimation of the efficiency of single counting events. (only in mode="eff")
    mean_efficiency_D : float
        Estimation of the efficiency of logic sum of double coincidences. (only in mode="eff")
    mean_efficiency_T : float
        Estimation of the efficiency of triple coincidences. (only in mode="eff")
    
    """
    
    e, p = readBetaShape(rad, 'beta-', 'tot')
    em=np.empty(len(e))
    for i, ei in enumerate(e):
        ed = energie_dep_beta(ei)
        em[i] = E_quench_e(ed*1e3,kB*1e3,ne)*1e-3
        
        
    if mode2=="sym":
        eff_S = sum(p*(1-np.exp(-L*em/3)))
        eff_T = sum(p*(1-np.exp(-L*em/3))**3)
        eff_D = sum(p*(3*(1-np.exp(-L*em/3))**2-2*(1-np.exp(-L*em/3))**3))
        TDCR_calcul=eff_T/eff_D
        res=(TDCR_calcul-TD)**2

    if mode2=="asym":
        # eff_A = sum(p*(1-np.exp(-L[0]*em/3)))
        # eff_B = sum(p*(1-np.exp(-L[1]*em/3)))
        # eff_C = sum(p*(1-np.exp(-L[2]*em/3)))
        eff_AB = sum(p*(1-np.exp(-L[0]*em/3))*(1-np.exp(-L[1]*em/3))) 
        eff_BC = sum(p*(1-np.exp(-L[1]*em/3))*(1-np.exp(-L[2]*em/3))) 
        eff_AC = sum(p*(1-np.exp(-L[0]*em/3))*(1-np.exp(-L[2]*em/3))) 
        eff_T = sum(p*(1-np.exp(-L[0]*em/3))*(1-np.exp(-L[1]*em/3))*(1-np.exp(-L[2]*em/3)))
        eff_D = sum(p*((1-np.exp(-L[0]*em/3))+(1-np.exp(-L[1]*em/3))+(1-np.exp(-L[2]*em/3))-2*(1-np.exp(-L[0]*em/3))*(1-np.exp(-L[1]*em/3))*(1-np.exp(-L[2]*em/3))))
        eff_S = sum(p*((1-np.exp(-L[0]*em/3))+(1-np.exp(-L[1]*em/3))+(1-np.exp(-L[2]*em/3))-((1-np.exp(-L[0]*em/3))+(1-np.exp(-L[1]*em/3))+(1-np.exp(-L[2]*em/3))-2*(1-np.exp(-L[0]*em/3))*(1-np.exp(-L[1]*em/3))*(1-np.exp(-L[2]*em/3)))-(1-np.exp(-L[0]*em/3))*(1-np.exp(-L[1]*em/3))*(1-np.exp(-L[2]*em/3))))
        TABmodel = eff_T/eff_AB
        TBCmodel = eff_T/eff_BC
        TACmodel = eff_T/eff_AC
        res=(TAB-TABmodel)**2+(TBC-TBCmodel)**2+(TAC-TACmodel)**2
    
    if mode == "res":
        return res
    if mode == "eff":
        return eff_S, eff_D, eff_T
    
def clear_terminal():
    """Function to clear the terminal screen
    """
    if os.name == "posix":
        os.system("clear")  # For UNIX/Linux/MacOS
    else:
        os.system("cls")    # For Windows

def display_header():
    """ Function to display the header.
    """
    clear_terminal()
    version = pkg_resources.get_distribution("tdcrpy").version
    header_text = r'''
 ______  ______  ______ _______  ________
|__  __||  ___ \|  ___||  ___ | |  ____ |
  | |   | |  | || |    | |  | | | |___| |___     ___
  | |   | |  | || |    | |__| | |  _____|\  \   |  |
  | |   | |__| || |____|  __  \ | |       \  \  |  |
  |_|   |_____/ |_____||_|  \__\|_|        \  \_|  |
  +++++++++++++++++++++++++++++++++++++++++/      /
  ________________________________________/      /
 |______________________________________________/     

'''
    header_text2 = "version "+version+"\n\
BIPM 2023 - license MIT \n\
distribution: https://pypi.org/project/TDCRPy \n\
developement: https://github.com/RomainCoulon/TDCRPy \n\n\
start calculation..."
 
 # Start Calculation
    print(header_text)
    print(header_text2)
