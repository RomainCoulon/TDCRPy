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
# import statsmodels.api as sm
import matplotlib.pyplot as plt
import time

absolutePath = False


def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)


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
      if "ALP " in decayData[-1]:    decayData[-1] = decayData[-1].replace("ALP ", "ALP; ") # type of the disintegration = alpha
      if "CK " in decayData[-1]:    decayData[-1] = decayData[-1].replace("CK ", "CK; ") # type of the disintegration = Electron Capture K shell
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
    
      if d[0] == "ALP" or d[0] == "BEM" or d[0] == "BEP" or d[0] == "CK" or d[0] == "CL1" or d[0] == "CL2" or d[0] == "CM" or d[0] == "CN":                         # Read information on the decay branch
        if d[0] == "ALP": particle.append("alpha")
        if d[0] == "BEP": particle.append("beta+")
        if d[0] == "BEM": particle.append("beta")
        if d[0] == "CK": particle.append("Atom_K")
        if d[0] == "CL1": particle.append("Atom_L1")
        if d[0] == "CL2": particle.append("Atom_L2")
        if d[0] == "CM": particle.append("Atom_M")
        if d[0] == "CN": particle.append("Atom_N")
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

# tic()
# readPenNuc("Co-60")
# toc() # 0.016 s

#===============================================================================================================



#================================== StoppingPower for alpha particle ===========================================

if absolutePath: f_alpha = open('G:\Python_modules\Jialin\MCTDCRcode\\alpha_toulene.txt')
else: f_alpha = open('alpha_toulene.txt')
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
    # rho: density of the absorber (g.cm-3)
    # e keV
    # energy keV
    # dEdx: keV.cm2/g
    energy_alpha = np.array(energy_alpha)
    dEdx_alpha = np.array(dEdx_alpha)
    dEdx = np.interp(e,energy_alpha ,dEdx_alpha)   
    return dEdx*rho                        #unit keV.cm-1

# tic()
# stoppingpowerA(5000,rho=0.96,energy_alpha=energy_alph,dEdx_alpha=dEdx_alph)
# toc() # 0 s
#===================================================================================================================


#===================================================================================================================


if absolutePath: file_TanXia=open('G:\Python_modules\Jialin\MCTDCRcode\\TandataUG.txt', "r")
else: file_TanXia=open('TandataUG.txt', "r")
data_TanXia=file_TanXia.read(); file_TanXia.close()
data_TanXia=data_TanXia.split("\n"); data_TanXia_f = np.empty(len(data_TanXia))
for i, x in enumerate(data_TanXia):
  if i<len(data_TanXia)-1: data_TanXia_f[i]=float(x)


# def stoppingpowerE(e,*,za=0.5459,rho=0.96,I=65,spmodel="Beth_Tan_Xia",emin=0,file=data_TanXia_f): # valid
#     c1=0.57
#     c2=1.16
#     Cte=0.307075*c1           # Relativist constant /(eV.cm^2) 
#     mec2=511000               # Mass energy of an electron /eV
#     dEdx=0                    # initial value of the stopping power
   
#     if spmodel=="Beth":
#         beta=np.sqrt(c2*e/mec2)            # relative velocity
#         if e>54:
#             dEdx=Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)    
#         else:
#             dEdx=0
    
#     if spmodel=="Beth_extr_1":
#         if e>400:
#             beta=np.sqrt(c2*e/mec2)            # relative velocity
#             if beta>=1:
#                 beta
#             dEdx=Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)    
#         else:
#             if e>emin:
#                 beta=np.sqrt(c2*400/mec2)       # relative velocity
#                 dEdx=20*Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)*e**-0.5
#             else:
#                 dEdx=0

#     if spmodel=="Beth_extr_2":
#         if e>100:
#             beta=np.sqrt(c2*e/mec2)            # relative velocity
#             dEdx=Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)    
#         else:
#             if e>emin:
#                 beta=np.sqrt(c2*100/mec2)       # relative velocity
#                 dEdx=0.01*Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)*e
#             else:
#                 dEdx=0
      
#     if spmodel=="Beth_extr_3":
#         if e>1000:
#             beta=np.sqrt(c2*e/mec2)            # relative velocity
#             dEdx=Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)
#         else:
#             if e>emin:
#                 beta=np.sqrt(c2*1000/mec2)     # relative velocity
#                 dEdx=1/(1000**-1.1)*Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)*e**-1.1
#             else:
#                 dEdx=0
    
#     if spmodel=="Beth_Tan_Xia":
#         if e>=20000:
#             beta=np.sqrt(c2*e/mec2)            # relative velocity
#             if beta>=1:
#                 beta=0.97
#             dEdx=Cte*rho*za*(beta**-2)*(np.log(mec2*beta**2/(I*(1-beta**2)))-beta**2)
#         else:
#             if e>emin:
#                 dEdx=(float(file[int(e)]))*1e6
#             else:
#                 dEdx=0
#     if dEdx<0:
#         dEdx=0
#     return dEdx




#===============================================================================================

#========================   Nouveau modèle pour calculer le pouvoir d'arrête d'électron ========

def stoppingpower(e,rho=0.96,Z=5.2,A=11.04,emin=0,file=data_TanXia_f):
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
plt.savefig('stoppingpowerE_A.png')

'''

#===================================================================================================


# def readBetaSpectrum(rad):
#     f=open("spectrum_"+rad+".txt", "r") # open the file 
#     ne = 1000   # number of energy bins
#     e=np.empty(ne) # declare an empty array to store the energy values
#     p=np.empty(ne) # declare an empty array to probability values
#     i=0
#     for line in f: # scan the line of the file
#         e[i]=float(line[:9]) # read the energy value and convert in float 
#         p[i]=float(line[9:]) # read the probability value and convert in float
#         i+=1
#     sumP=sum(p) # sum of probabilities
#     p=p/sumP # re-normilize to make the sum of probabilities equal to 1
#     f.close() #close the file
#     return e, p, ne


#=============================================================================================

#====================  Fonction pour lire BetaShape   ========================================

def readBetaShape(rad,mode,trans):
    # mode(str): 'beta-','beta+'
    # trans(str):'trans0','trans1' ....
    file = "All-nuclides_BetaShape.zip"
    z = zf.ZipFile(file)
    Rad = rad.replace('-','')
    name_doc = Rad+'/'+mode+'_'+Rad+'_'+trans+'.bs'
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
    return e, dNdx

# tic()
# readBetaShape("Co-60", "beta-", "tot")
# toc() # 0.016 s



#=====================================================================================

## Display beta spectra
# out = readBetaShape("H-3", "beta-", "tot")
# print(sum(out[1]))
# plt.figure("beta spectrum")
# plt.plot(out[0],out[1],label='Beta spectrum',ls='-',color='red',lw=3)
# plt.xscale('linear')
# plt.legend(fontsize=12,loc='best')
# plt.xlabel(r'$E$ /keV', fontsize=12)
# plt.ylabel(r'd$N$/d$E$ /keV$^{-1}$',fontsize=12)
# plt.savefig("BetaSpectrum.png")


#=======================================================================================

#============================  Fonction quenching  =====================================

def E_quench_e(e,kB): # e : eV  kB:cm/MeV
    e_dis = np.linspace(0,e,10000)
    delta = e_dis[2] - e_dis[1]
    q = 0
    for i in e_dis:
        q += delta/(1+kB*stoppingpower(i))
    return q #eV

def E_quench_a(e,kB): # e : keV   kB:cm/keV
    e_dis = np.linspace(1,e,10000)
    delta = e_dis[2] - e_dis[1]
    q = 0
    for i in e_dis:
        q += delta/(1+kB*stoppingpowerA(i))
    return q #keV

#=========================================================================================



#print(E_quench_a(5.5e3,1e-5))

#======================================================================================

#========================= Tracer les courbes avec kB différents ======================

'''
s1 = []
s2 = []
s3 = []
x = np.linspace(10,1e4,1000) 

for i in x:
    s1.append(E_quench_e(i,kB=7e-3)*1e-3)
    s2.append(E_quench_e(i,kB=1e-2)*1e-3)
    s3.append(E_quench_e(i,kB=1.4e-2)*1e-3)

plt.plot(x,s2,label='E_quenched/E_0.01',ls=':',color='red',lw=3)
plt.plot(x,s1,label='E_quenched/E_0.007',color='green',lw=2)
plt.plot(x,s3,label='E_quenched/E_0.014')
#plt.xscale('log')
#plt.yscale('log')
plt.legend(fontsize=12,loc='best')
plt.xlabel('E_emitted/eV')
plt.ylabel('quenching energy/E_emitted /keV')
plt.savefig("quenching E_test 10-10k.png")
'''
#============================================================================================

#============================================================================================

#========================= énergie gamma ===================================================
#'''
if absolutePath: 
    f1 = open('G:\Python_modules\Jialin\MCTDCRcode\\MCNP-MATRIX/matrice/matrice_p_1_200k.txt')
    f2 = open('G:\Python_modules\Jialin\MCTDCRcode\\MCNP-MATRIX/matrice/matrice_p_200_2000k.txt')
    f3 = open('G:\Python_modules\Jialin\MCTDCRcode\\MCNP-MATRIX/matrice/matrice_p_2000_10000k.txt')
    fe = open("G:\Python_modules\Jialin\MCTDCRcode\\MCNP-MATRIX/matrice/E_depose.txt")   
else:
    f1 = open('MCNP-MATRIX/matrice/matrice_p_1_200k.txt')
    f2 = open('MCNP-MATRIX/matrice/matrice_p_200_2000k.txt')
    f3 = open('MCNP-MATRIX/matrice/matrice_p_2000_10000k.txt')
    fe = open("MCNP-MATRIX/matrice/E_depose.txt")

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

#'''

def energie_dep_gamma(e_inci,*,matrice1=Matrice1,matrice2=Matrice2,matrice3=Matrice3,ed=Matrice_e):
    ## sort keV / entrée : MeV
    if e_inci <= 200:
        index = int(e_inci)            # index de colonne de la matrice de l'énergie incidente la plus proche 
        #doc = 'MCNP-MATRIX/matrice/matrice_p_1_200k.txt'
        matrice = matrice1
        #taille_x = 200
        e = ed[:,0]
    
    elif e_inci <= 2000:
        index = int((e_inci-200)/2)
        #doc = 'MCNP-MATRIX/matrice/matrice_p_200_2000k.txt'
        matrice = matrice2
        #taille_x = 901
        e = ed[:,1]

    else:
        index = (int(e_inci)-2000)//10
        #doc = 'MCNP-MATRIX/matrice/matrice_p_2000_10000k.txt'
        matrice = matrice3
        #taille_x = 801
        e = ed[:,2]
    
    '''
    with open(doc) as f:
        data = f.readlines()
    
        matrice = np.zeros((1002,taille_x))

        for i in range(1002):
            data[i] = data[i].split()
            for j in range(taille_x):
                matrice[i][j] = float(data[i][j])
    '''
    inde = sampling(matrice[1:,index])
    if inde == 1 : result = 0
    else: result = e[inde]*1e3*e_inci/matrice[0][index]
    return result

'''
r = []   
for i in range(100):
    r.append(energie_dep_gamma(2568))
# print(r)
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
def readEShape(rad):
    file = 'All-nuclides_Ensdf.zip'
    z = zf.ZipFile(file)
    name = rad + '.txt'
    with z.open(name) as f:
        data = f.readlines()

        for i in range(np.size(data)):
            data[i] = str(data[i])
            data[i] = data[i].replace("b",'')
            data[i] = data[i].replace("\\r\\n",'')
            data[i] = data[i].replace("'",'')
        
        for i in range(np.size(data)):
            data[i] = data[i].split()
        
        '''
        index_decay = []
        index_auger = []
        index_end = []

        for i,p in enumerate(data):
            if 'DECAY' in p:
                i_decay = i
                index_decay.append(i_decay)
        
            if 'Auger' in p:
                i_auger = i
                index_auger.append(i_auger)
                p = data.index('P',i_auger,i_decay)
                index_end.append(p)
            else:
                print('no auger electrons')
                break

        '''
    index_decay = []
    index_auger = []
    index_end = []
    for i,p in enumerate(data):
        if 'IT' in p:
            print('isomeric transition')
        if 'DECAY' in p:
            index_decay.append(i)
        if 'Auger' in p:
            #i_auger = i
            index_auger.append(i)
        if 'P' in p:
            i_p = i
            index_end.append(i_p)
    proba = []
    Type = []
    #print(index_decay,index_auger,index_end)
    #'''
    for i in range(len(index_auger)):
        start = index_auger[i]
        end = index_end[i]
        pb = []
        ty = []
        for j in range(start+3,end-1):
            if len(data[j]) <=2:continue
            if 'AUGER' in data[j]:
                if len(data[j][2])>6:
                    d = data[j][2].split('-')
                    data[j][2] = round((float(d[1])+float(d[0]))/2,3)
                pb.append(data[j][2])
                ty.append(data[j][-2])
            else:
                pb.append(data[j][2])
                ty.append(data[j][-1])
                #print(data[j][2],data[j][-1])
        proba.append(pb)
        Type.append(ty)    
    #'''
    
    return proba,Type

pr,typ = readEShape('Ag-108m')
print(pr,typ)