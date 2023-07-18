# -*- coding: utf-8 -*-
"""
Created on Mon May 15 15:41:41 2023

@author: romain.coulon
"""



import zipfile
import numpy as np

file = 'decayData/All-nuclides_Ensdf.zip'
z = zipfile.ZipFile(file)

with z.open('Ag-108m.txt') as f:
    data = f.readlines()
    nL = np.size(data)
    
    ## fomatage en list de vecteur pour chaque line du fichier ENSDF
    for i in range(nL):              # boucle dans data
        data[i] = str(data[i])                  # conversion en string
        data[i] = data[i].replace("b'",'')       # on enlève les caratères de formatage
        data[i] = data[i].replace("\\r\\n'",'')  # on enlève les caratères de formatage
    for i in range(nL):
        data[i] = data[i].split()
    
    for i in range(nL):
        if i>0 and ("L" in data[i]) and ("AUGER" in data[i]) and ("|]" in data[i-1]):
            data.insert(i,[data[i][0], 'T'])

    ## Définir des repères
    index_auger = []          # indice de
    daugther = []             # liste des filles
    index_end = []            # indice dernier valeur pour une fille
    posi = []                 # indice des blocs transitions 
    for i, p in enumerate(data): # boucle dans le fichier
        if 'DECAY' in p:
            daugther.append(p[0])   # liste des filles
        if 'Auger' in p:
            index_auger.append(i)   # repère de ligne pour ensuite définir l'ensemble des blocs
        if len(p)==2:
            posi.append(i)          # repère de ligne pour ensuite définir les blocs (2)
        # if 'L' in p
        #     posi:
        if 'P' in p:
            index_end.append(i)     # repère de fin pour ensuite définir l'ensemble des blocs
            posi.append(i)
            
    ## Filtrage des données utiles et formatage
    energy = []
    prob = []
    Type = []
    Fille = []
    for i in range(len(posi)-1): # bloucle dans les blocs
        start = posi[i]+1        # indice debut du bloc
        end = posi[i+1]          # indice fin du bloc
        d = data[start:end]      # Le bloc
        e_b = []                   # Energie des transitions du bloc
        prob_b = []                # Proba des transitions du bloc
        Type_b = []                # Types des transitions du bloc
        Fille_b = []               # Fille des transitions du bloc
        #e2 = []
        if start == end:
            continue
        if start-1 in index_end:
            #e = []
            continue
        for n, p1 in enumerate(d):
            #e2 = []
            if '-' in p1[2]: # traitement des intervals de nombres réels
                x = p1[2].split('-')
                p1[2] = round((float(x[1])+float(x[0]))/2,3) # Moyenne
            if '|]' in p1:        # traitment des accolades 
                if len(p1)>6:     # repère de la ligne portant la proba pour le groupe 
                    prob_b.append(float(p1[4]))
                e_b.append(float(p1[2]))
            else:
                e_b.append(float(p1[2]))
                prob_b.append(float(p1[3]))
        
        if len(prob_b)==1 and len(e_b)>1:
            e_b = [np.mean(e_b)]
        
        print(e_b)
        print(prob_b)
        print(" ")
