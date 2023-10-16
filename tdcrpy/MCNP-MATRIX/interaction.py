# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:16:55 2023

@author: romain.coulon
"""
import numpy as np
import matrice as m
import matplotlib.pyplot as plt



def readMCNP(energy,niveau,par,v,npas=1000,mode='N'):
    '''
    ******************************************
    pour lire les fichiers de résultat de MCNP
    ******************************************
    ---------
    PARAMETRE
    ---------
    energy ------- type: int ---- énergie incidente (keV) 

    niveau ------- type: int ---- l'ordre de matrice 
        0: 1ère matrice -- 1-200 keV   
        1: 2ème matrice -- 200-2000 keV
        2: 3ème matrice -- 2000-10000 keV

    par ---------- type: str ---- type de particule
        p: photon
        b: béta- (électron)
        bp: béta+ (positron)

    npas --------- type: int ---- nombre de bins

    mode --------- type: str ---- noraliser ou pas
        N: non
        Y: yes

    ------
    RETURN
    ------
    e --------- type: array ----- le vecteur d'énergie déposée
    p --------- type: array ----- le vecteur de proba d'énergie déposée
    
    '''

    output = 'output/output'+str(v)+'ml/'
    if par == 'p':
        name1 = 'photon/'
    elif par=='b':
        name1 = 'beta-/'
    else:
        name1='beta+/'

    if niveau == 0:
        name_doc = output + name1 + 'output_p_1-200/output_'+str(int(energy))+'keV.o'
    elif niveau == 1:
        name_doc = output + name1 + 'output_p_200-2000/output_'+str(int(energy))+'keV.o'
    else:
        name_doc = output + name1 + 'output_p_2000-10000/output_'+str(int(energy))+'keV.o'

    e = []
    p = []
    f = open(name_doc)
    data = f.readlines()
    f.close()

    for i in range(len(data)):
        m = i
        if data[i].find('cell  5') == 1: break


    end_point = m + 2 + npas + 2
    for j in range(m+2,end_point):
        data[j] = data[j].split()
        e.append(float(data[j][0]))
        p.append(float(data[j][1]))
    
    if mode=='Y':    
        p /= sum(np.asarray(p)) # normaliser p
    return e,p

E=np.arange(1,197,1)
niv = 0
par = "p"
v=10

p_photoelec=[]
for Es in E:
    out  = readMCNP(Es,niv,par,v,npas=1000,mode='N')
    
    ve = np.asarray(out[0])*1000
    vp = np.asarray(out[1])
    
    i=np.where(ve > Es)[0][0]
    p_photoelec.append(vp[i])
    print(vp[i])

"""UG"""
# fraction relative des atomes
pH = 0.578772
pC = 0.338741
pN = 0.000302
pO = 0.082022
pP = 0.000092
pCl = 0.000071
          
# energie des rayons x
xK_H = 0
xK_C = 0.277
xK_N = 0.3924
xK_O = 0.5249
xK_P = 2.0137
xK_Cl = 2.62239

plt.figure('Photoelectric prob')
plt.clf()
plt.plot(E,p_photoelec,'-p',label=r' ')
plt.xlabel(r'$E$ /keV', fontsize=16)
plt.ylabel(r'$p$ /keV-1', fontsize=16)
plt.xticks(fontsize=14)
plt.legend(fontsize=12)


plt.figure('Spectrum')
plt.clf()
plt.plot(ve,vp,'-p',label=r' ')
plt.plot(ve[i],vp[i],"+k")
plt.xlabel(r'$E$ /keV', fontsize=16)
plt.ylabel(r'$p$ /keV-1', fontsize=16)
plt.yscale("log")
plt.xticks(fontsize=14)
plt.legend(fontsize=12)






