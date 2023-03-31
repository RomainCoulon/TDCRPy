# code pour créer un fichier pour la matrice d'énergie incidente et d'énergie déposée et la matrice CDF

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors

start_energy = 0      #keV --- debut de enrgie incidente
end_energy = 20      #keV --- fin de energie incidente

if end_energy <= 200:
     delta_incident = 1    #keV --- delta de energie incidente
     #npas = 1000  # nombre de pas 
     #delta_E = 0.2         #keV --- delta_E de l'energie deposee
elif end_energy <=2000:
     delta_incident = 1    #keV --- delta de energie incidente
     #npas = 2000  # nombre de pas 
     #delta_E = 1         #keV --- delta_E de l'energie deposee
else:
     delta_incident = 10    #keV --- delta de energie incidente
     #npas = 1000  # nombre de pas 
     #delta_E = 10         #keV --- delta_E de l'energie deposee

npas = 1000  # nombre de pas 

#========================================================================================

#===================== fonction pour lire output du MCNP ================================

# attention:éviter 0keV pour la fonction readMCNP
def readMCNP(energy,npas,mode='N'):
    e = []
    p = []
    f = open('output/output_'+str(int(energy))+'keV.o')
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



taille_x = int((end_energy - start_energy + 1)/delta_incident)
taille_y = npas+1
energy_inci = np.linspace(start_energy,end_energy,taille_x)
matrice_p = np.zeros((taille_y,taille_x))
matrice_cdf = np.zeros((taille_y,taille_x))
#print(taille_x)

for i in range(taille_x):
    energy = energy_inci[i]
    if energy==0.:continue
    e,p = readMCNP(energy,npas)
    cdf = np.cumsum(p)
    matrice_p[0][i] = energy
    matrice_cdf[0][i] = energy
    for j in range(1,taille_y):
        matrice_p[j][i] = p[j-1]
        matrice_cdf[j][i] = cdf[j-1]

if energy_inci[0] == 0.:
    matrice_p = np.delete(matrice_p,0,axis=1)
    matrice_cdf = np.delete(matrice_cdf,0,axis=1)
    taille_x = taille_x - 1
    #print(energy)
#print(matrice_p[0][:])
x = matrice_p.shape[1]
y = matrice_p.shape[0]
#print(x,y)

for i in range(x):
    x1 = matrice_p[0][i]
    #a = energy_inci[i]
    #x = np.linspace(a,a,1002)
    col = tuple(np.random.rand(3))
    for j in range(1,y):
        y1 = matrice_p[j][i]
        plt.scatter(x1,y1,color=col)  
        #print(x1,y1)
plt.yscale('log')
plt.xticks(np.arange(1,21,1))  
plt.savefig('matrice/matrice.png')



'''
with open('matrice/matrice_0_20k.txt','w') as file:
    file.write('# matrice energy\n')
    for i in range(taille_y):
        for j in range(taille_x):
            file.write("%e"%matrice_p[i][j])
            file.write('         ')
        file.write('\n')
    file.write('\n')
    file.write('\n')
    for i in range(taille_y):
        for j in range(taille_x):
            file.write("%e"%matrice_cdf[i][j])
            file.write('         ')
        file.write('\n')
'''   