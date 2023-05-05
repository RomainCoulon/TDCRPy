# code pour créer un fichier pour la matrice d'énergie incidente et d'énergie déposée et la matrice CDF

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from matplotlib.pyplot import MultipleLocator


start_energy = 0         # keV --- debut de enrgie incidente
end_energy = 20          # keV --- fin de energie incidente

'''
if end_energy <= 200:       # 1k - 200k
     delta_incident = 1     # keV --- delta de energie incidente
     #n_ei = 200            # nombre de pas de l'énergie incidente
     #npas = 1000           # nombre de pas de l'énergie déposée
     #delta_E = 0.2         # keV --- delta_E de l'energie deposee

elif end_energy <=2000:     # 200k - 2M
     delta_incident = 2     # keV --- delta de energie incidente
     #n_ei = 1800           # nombre de pas de l'énergie incidente (200k - 2M)
     #npas = 1000           # nombre de pas de l'énergie déposée
     #delta_E = 2           # keV --- delta_E de l'energie deposee

else:                       # 2M - 10M 
     delta_incident = 10    # keV --- delta de energie incidente
     #n_ei = 800            # nombre de pas de l'énergie incidente (2M - 10M)
     #npas = 1000           # nombre de pas de l'énergie déposée 
     #delta_E = 10          # keV --- delta_E de l'energie deposee
'''
#npas = 1000                 # nombre de pas de l'énergie déposée

#========================================================================================

#===================== fonction pour lire output du MCNP ================================

# attention:éviter 0keV pour la fonction readMCNP
def readMCNP(energy,niveau,par,npas=1000,mode='N'):
    if par == 'p':
        name1 = 'p/'
    else:
        name1 = 'b/'

    if niveau == 0:
        output = 'output1'
    elif niveau == 1:
        output = 'output2'
    else:
        output = 'output3'

    name_doc = output + name1 + 'output_' + str(int(energy))+'keV.o'
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


# ---------------------- créer la matrice de pdf et cdf -----------------------

def creat_matrice(niveau,par,mode='pdf',npas=1000):
    # niveau : 0 -- énergie basse (0-200k); 1 -- énergie moyenne (200k-2M); 2 -- haute énergie (2M-10M)
    # mode : pdf -- probabilité ; cdf -- cumsum
    PAR = par
    NIVEAU = niveau
    if niveau == 0:
        end_energy = 200
        start_energy = 1
        taille_x = 200      # 1k-200k où delta=1k

    elif niveau == 1:
        end_energy = 2000
        start_energy = 200
        taille_x = 901     # 200k-2000k où delta=2k

    else:
        end_energy = 10000
        start_energy = 2000
        taille_x = 801      #2M-10M où delta = 0.1M

    taille_y = npas+3
    energy_inci = np.linspace(start_energy,end_energy,taille_x)
    matrice_p = np.zeros((taille_y,taille_x))
    #matrice_cdf = np.zeros((taille_y,taille_x))

    for i in range(taille_x):
        energy = energy_inci[i]
        e,p = readMCNP(energy,niveau=NIVEAU,par=PAR)
        #if mode == "cdf":
         #   p = np.cumsum(p)
        matrice_p[0][i] = energy       # entete de matrice
        #matrice_cdf[0][i] = energy
        for j in range(1,taille_y):
            matrice_p[j][i] = p[j-1]
            #matrice_cdf[j][i] = cdf[j-1]

    return e,matrice_p

'''
if energy_inci[0] == 0.:
    matrice_p = np.delete(matrice_p,0,axis=1)
    matrice_cdf = np.delete(matrice_cdf,0,axis=1)
    taille_x = taille_x - 1
'''

#print(energy)
#print(matrice_p[0][:])
#x = matrice_p.shape[1]
#y = matrice_p.shape[0]
#print(x,y)

#'''
def matrice_fig(matrice_p,start,end,e,par): 
    # matrice_p : la matrice complète de chaque gamme
    # vecteur de l'énergie déposée pour chaque gamme
    # start et end en keV
    if par == 'p':
        particle = 'gamma'
    elif par == 'b':
        particle = 'beta-'
    elif par == 'bp':
        particle = 'beta+'


    if end <= 200:      # delta_Ei = 1
        delta_Ei = 1
        #d_end = int(end/201*1000)
        i_st = int(start-1)
        i_end = int(end-1)
        #x = i_end - i_st +1
        end_point = 2 + 5*(end+1)

    elif end <= 2000:    # delta_Ei = 1
        delta_Ei = 2
        #d_end = int(end/2001*1000)
        i_st = int((start-200)/delta_Ei)
        i_end = int((end-200)/delta_Ei)
        #x = i_end - i_st + 1
        end_point = int(end/2+1)

    else:
        delta_Ei = 10
        #d_end = int(end/(1e4+1)*1000)
        i_st = int((start-2000)/10)
        i_end = int((end-2000)/10)
        #x = i_end - i_st +1
        end_point = int(end/10+1)

    #zz = matrice_p[1:d_end+1,i_st:i_end+1]  # matrice 
    zz = matrice_p[1:,i_st:i_end+1]
    #zz = matrice_p[1:end_point,i_st:i_end+1]
    n = int((end-start)/delta_Ei)
    x = np.linspace(start,end,n+1)
    y = e
    #y = e[0:end_point-1]  
    #x_size = np.size(x)
    #y_size = np.size(y)

    '''
    for i in range(x):
        xx[:,i] = i * delta_Ei + start
        yy[:,i] = e[0:d_end]
    '''
    '''
    if end>20:
        zz = zz[2:,:]
        xx = xx[2:,:]
        yy = yy[2:,:]
    
    ''' 
    if end > 20 and par=='p':
        zz = zz[2:,:]
        y = y[2:]
    
    #print(xx.shape,yy.shape,zz.shape)
    #print(zz[0,0],yy[0,0])
    xx,yy = np.meshgrid(x,y)
    print(xx[0,0],yy.shape,zz[0,0])
    h = plt.pcolormesh(xx,yy,zz,cmap = plt.cm.hot)
    cb = plt.colorbar(h)
    cb.set_label("probabilité en log")
    #plt.xticks(xs,s)
    #plt.yticks(np.linspace(0,end,10))
    #x_maj = MultipleLocator(1)
    #ax = plt.gca()
    #ax.xaxis.set_major_locator(x_maj)
    #ax.yaxis.set_major_locator(LogLocator(base=10))
    #plt.yscale('log')
    plt.ylim(y[0],y[end_point])
    #plt.xlim(start,end)
    plt.xlabel("énergie incidente/keV")
    plt.ylabel("énergie déposée/MeV")
    title = "probabilité d'énergie déposée par " + particle + " de " + str(start) + "-" + str(end) + "keV"
    plt.title(title)
    name = "matrice/matrice_fig_log6-6" + particle +"_" + str(start) + "_" + str(end) + "k.png"
    plt.savefig(name)
    return 0
#'''




def ecrit_matrice(matrice,niveau,par):
    if par == 'p':
        name1 = 'gamma_'
    elif par == 'b':
        name1 = 'beta-_'
    elif par == 'bp':
        name1 = 'beta+_'

    if niveau == 0:
        end_energy = 200
        start_energy = 1
        taille_x = 200      # 1k-200k où delta=1k

    elif niveau == 1:
        end_energy = 2000
        start_energy = 200
        taille_x = 901     # 200k-2000k où delta=1k

    else:
        end_energy = 10000
        start_energy = 2000
        taille_x = 801      #2M-10M où delta = 0.1M

    taille_y = 1003
    name = 'matrice/matrice_' + name1 + str(start_energy) + '_' + str(end_energy) + 'k.txt'
    with open(name,'w') as file:
    #file.write('# matrice energy\n')
        for i in range(taille_y):
            for j in range(taille_x):
                file.write("%e"%matrice[i][j])
                file.write('         ')
            file.write('\n')
        #file.write('\n')
    
        '''
        file.write('\n')
        for i in range(taille_y):
            for j in range(taille_x):
                file.write("%e"%matrice_cdf[i][j])
                file.write('         ')
            file.write('\n')
        '''


def find_info(energy,niveau,par,npas=1000,mode='N'):
    if par == 'p':
        name1 = 'p/'
    else:
        name1 = 'b/'

    if niveau == 0:
        output = 'output1'
    elif niveau == 1:
        output = 'output2'
    else:
        output = 'output3'

    name_doc = output + name1 + 'output_' + str(int(energy))+'keV.o'
    e = []
    p = []
    f = open(name_doc)
    data = f.readlines()
    f.close()

    for i in range(len(data)):
        m = i
        if data[i].find('bremsstrahlung') == 1: break 
    return m

position = find_info(206,1,'b')

#================ tracer la matrice ========================================
#e,matrice_p = creat_matrice(1,par='b')
#print(e[-1],matrice_p.shape,np.size(e))
#ecri = ecrit_matrice(matrice_p,1,par='b') 
#fig1 = matrice_fig(matrice_p,41,80,e)
#for i in range(1003):
    #for j in range(901):
        #matrice_p[i][j] = np.log(matrice_p[i][j]+6e-6)
#fig2 = matrice_fig(matrice_p,200,900,e,'b')
#print(fig2)

#============= tracer la proba à Ei=Ed ====================================== 
'''
max_valeur = []
for i in range(901):
    a = matrice_p[1:,i]
    max_valeur.append(max(a))
e_inci = matrice_p[0,:]
#print(e_inci[0:10],max_valeur[0:10])
plt.plot(e_inci,max_valeur)
plt.xticks(np.arange(200,2002,200))
plt.yticks(np.arange(0.4,1.05,0.05))
plt.xlabel('énergie incidente/keV')
plt.ylabel('probabilité')
plt.title("probabilté d'énergie déposée de particule béta à Ei=E_déposée")
plt.savefig('beta.png')
'''













'''
name = 'matrice/matrice_' + str(start_energy) + '_' + str(end_energy) + 'k.txt'
with open(name,'w') as file:
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