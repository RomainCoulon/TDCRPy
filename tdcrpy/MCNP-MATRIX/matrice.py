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

    if par == 'p':
        name1 = 'p/'
    elif par=='b':
        name1 = 'b/'
    else:
        name1='bp/'

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

#e1, p1 = readMCNP(1000,1,'b')
#e2, p2 = readMCNP(1000,1,'bp')  
#plt.plot(e1,p1,label='beta-')
#plt.plot(e2,p2,label='beta+')
#plt.yscale('log')
#plt.legend(fontsize=10,loc='best')
#plt.savefig('matrice/spectre comparison.png')
#print(e[-10:])
# ---------------------- créer la matrice de pdf et cdf -----------------------

def creat_matrice(niveau,par,mode='pdf',npas=1000):
    '''
    **************************************************
    pour former une matrice de proba d'énergie déposée
    **************************************************
    ---------
    PARAMETRE
    ---------
    niveau ------- type: int ---- l'ordre de matrice 
        0: 1ère matrice -- 1-200 keV   
        1: 2ème matrice -- 200-2000 keV
        2: 3ème matrice -- 2000-10000 keV
    
    par ---------- type: str ---- type de particule
        p: photon
        b: béta- (électron)
        bp: béta+ (positron)

    mode --------- type: str ---- type de distribution
        pdf: distribution de probabilite
        cdf: distribution de probabilité cumulée

    npas --------- type: int ---- nombre de bins 

    ------
    RETURN
    ------
    e ---------- type: array ------ le vecteur d'énergie déposée
    matrice_p -- type: array n*1000 dimension --- la matrice de probabilité dont l'entête est l'énergie incidente

    '''

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
        #print(energy)
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



def matrice_fig(matrice_p,start,end,e,par,v):
    '''
    *************************************************
    pour tracer la matrice de proba d'énergie déposée
    *************************************************
    ---------
    PARAMETRE
    ---------
    matrice_p ---- type: array n*1000 dimension ---- la matrice de probabilité dont l'entête est l'énergie incidente
    start -------- type: int ------ définir la gamme à tracer
    end ---------- type: int ------ définir la gamme à tracer
    e ------------ type: array 1*1000 dimension ------ le vecteur de chaque bin d'énergie
    par ---------- type: str ------ type de particule
        p: photon gamma
        b: béta-
        bp: béta+
    v ------------ type: int ------- volume de solution
    ''' 
    # matrice_p : la matrice complète de chaque gamme
    # vecteur de l'énergie déposée pour chaque gamme
    # start et end en keV
    if par == 'p':
        particle = 'photon'
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
        if end < 200:
            end_point = 2 + 5*(end+1)
        elif end ==200:
            end_point = -1    

    elif end <= 2000:    # delta_Ei = 1
        delta_Ei = 2
        #d_end = int(end/2001*1000)
        i_st = int((start-200)/delta_Ei)
        i_end = int((end-200)/delta_Ei)
        #x = i_end - i_st + 1
        if end < 2000:
            end_point = int(end/2+1)
        elif end == 2000:
            end_point = -1

    else:
        delta_Ei = 10
        #d_end = int(end/(1e4+1)*1000)
        i_st = int((start-2000)/10)
        i_end = int((end-2000)/10)
        #x = i_end - i_st +1
        if end < 10000:
            end_point = int(end/10+1)
        elif end == 10000:
            end_point = -1
        

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
    if end > 50 and par=='p':
        zz = zz[2:,:]
        y = y[2:]
    #print(xx.shape,yy.shape,zz.shape)
    #print(zz[0,0],yy[0,0])
    xx,yy = np.meshgrid(x,y)
    #print(xx[0,0],yy.shape,zz[0,0])
    h = plt.pcolormesh(xx,yy,zz,cmap = plt.cm.hot)
    cb = plt.colorbar(h)
    cb.set_label("probabilité")
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
    name = "matrice/matrice_" +str(v)+ 'ml-' + particle +"_" + str(start) + "_" + str(end) + "k.png"
    plt.savefig(name)
    return 0
#'''




def ecrit_matrice(matrice,niveau,par,v):
    '''
    ***************************************************************************
    pour écrire une matrice complète de proba d'énergie déposée dans un fichier
    ***************************************************************************
    ---------
    PARAMETRE
    ---------
    matrice ----- type: array n*1000 dimension ---- matrice complète à écrire dans un fichier
    niveau ------ type: int ---- l'ordre de matrice 
        0: 1ère matrice -- 1-200 keV   
        1: 2ème matrice -- 200-2000 keV
        2: 3ème matrice -- 2000-10000 keV
    par --------- type: str ---- type de particule
        p: photon
        b: béta- (électron)
        bp: béta+ (positron)
    '''
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
    name = 'matrice/matrice_' + str(v) + name1 + str(start_energy) + '_' + str(end_energy) + 'k.txt'
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


def find_info(niveau,par,info,npas=1000,mode='N'):
    '''
    ********************************************************
    pour trouver les valeurs particulières que l'on s'intéresse dans un fichier de MCNP
    ********************************************************
    ---------
    PARAMETRE
    ---------
    niveau ------- type: int ---- l'ordre de matrice 
        0: 1ère matrice -- 1-200 keV   
        1: 2ème matrice -- 200-2000 keV
        2: 3ème matrice -- 2000-10000 keV
    
    par ---------- type: str ---- type de particule
        p: photon
        b: béta- (électron)
        bp: béta+ (positron)

    info -------- type: str ----- mot clé que l'on s'intéresse
    npas -------- type: int ----- nombre de bins
    mode

    ------
    RETURN
    ------
    nombre ------ type: list ----- l'ensemble de valeurs intéressantes
    
    '''
    if par == 'p':
        name1 = 'p/'
    else:
        name1 = 'b/'

    if niveau == 0:
        output = 'output1'
        energy_start = 1
        energy_end = 200
        nb = 200

    elif niveau == 1:
        output = 'output2'
        energy_start = 200
        energy_end = 2000
        nb = 901
    else:
        output = 'output3'
        energy_start = 2000
        energy_end = 10000
        nb = 801
    
    energy = np.linspace(energy_start,energy_end,nb)
    nombre = []
    for i in energy:
        name_doc = output + name1 + 'output_' + str(int(i))+'keV.o'
        f = open(name_doc)
        data = f.readlines()
        f.close()
        #m = data[info].split()
        #      for i in range(len(data)):
        #       m = i
        #      if data[i].find(info) == 1: break

        # lis = data[m].split() 

        #m = []
        for i in range(len(data)):
            data[i] = data[i].split()
        for i, p in enumerate(data):
            if info in p:
                m=i
                #nombre.append(float(data[m][1]))
                break
        if data[m+7][0] == '2':
            nombre.append(float(data[m+7][2]))
        #else:
            #print('bad position')
    return nombre

#================== tracer nb breamsstrahlung en fct de Ei ========================
#nb_brem = []
#for i in np.linspace(200,2000,901):
    #nb_brem.append(find_info(i,1,'b')/1e6)
#print(nb_brem[0:10])
#print(nb_brem[525:550])
#plt.plot(np.linspace(200,2000,901),nb_brem)
#plt.yticks(np.arange(12710,238780,20000))
#plt.xticks(np.arange(200,2002,200))
#plt.xlabel('énergie incidente/keV')
#plt.ylabel('nombre de photon produits par Bremsstrahlung')
#plt.title('nombre de photons produits par Bremsstrahlung en fonction de Ei')
#plt.savefig('evolution nb brem.png')
#print(value)

#=================tracer le nb electron qui sortent ========================
#a = np.array(find_info(1,'b','1electron'))
#print(a[899:901])
#plt.plot(np.linspace(200,2000,901),a/1e6)
#plt.xticks(np.arange(200,2002,200))
#plt.yticks(np.linspace(0,0.6,13))
#plt.xlim(200,2000)
#plt.xlabel('énergie initiale/keV')
#plt.ylabel("fraction d'electrons atteignant le flacon")
#plt.title("fraction d'electrons atteignant le flacon")
#plt.savefig('beta_atteind_flacon.png')

#===================================
#a = np.array(find_info(1,'b','bremsstrahlung'))
#print(a[899:901])
#plt.plot(np.linspace(200,2000,901),a/1e6)
#plt.xticks(np.arange(200,2002,200))
#plt.yticks(np.linspace(0,0.25,11))
#plt.xlim(200,2000)
#plt.xlabel('énergie incidente/keV')
#plt.ylabel("fraction d'photon Bremsstralung")
#plt.title("fraction d'photon Bremsstralung en fonction Ei")
#plt.savefig('nombre photon brem.png')




#================ tracer la matrice ========================================
#e,matrice_p = creat_matrice(0,par='p')
#print(len(e))
#print(matrice_p[76:90,25])
#print(matrice_p.shape,matrice_p[0,541:543])
#ecri = ecrit_matrice(matrice_p,1,'b',16) 
#fig1 = matrice_fig(matrice_p,50,100,e,'p',10)
#for i in range(1003):
  #  for j in range(901):
 #       matrice_p[i][j] = np.log(matrice_p[i][j]+1e-8)
#fig2 = matrice_fig(matrice_p,200,1000,e,'b',v=16)
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
plt.savefig('beta1.png')
'''


#===============================ecrire un fichier===================================
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
