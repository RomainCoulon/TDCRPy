import numpy as np
import linecache as lc
import matplotlib.pyplot as plt

def readMCNP(energy,np1):
    e = []
    p = []
    f = open('output/output_'+str(int(energy))+'keV.o')
    data = f.readlines()
    f.close()

    for i in range(len(data)):
        m = i
        if data[i].find('cell  5') == 1: break


    end_point = m + 2 + np1 + 2
    for j in range(m+2,end_point):
        data[j] = data[j].split()
        e.append(float(data[j][0]))
        p.append(float(data[j][1]))
        
    p /= sum(np.asarray(p)) # normaliser p
    return e,p

#e,p = readMCNP(10,1000)
#print(sum(p))

'''
plt.plot(e,p)
plt.xlim(-0.001,0.011)
plt.xlabel('Energy/MeV')
plt.ylabel('Probability')
plt.yscale('log')
plt.title('spectre gamma at E=100keV')
plt.savefig('proba = f(E)_test.png')
'''

#f  =open('matrice/matrice_0_20k.txt')
#data = f.readlines()

e,p = readMCNP(1986,1000)
pc = np.cumsum(p[2:])
for j in range(100):
    trial = float(np.random.rand(1))
    for i, p in enumerate(pc):
        if p> trial: break
    print(trial,i,e[i+2])
'''
me = np.zeros((1002,3))
x = [4,204,2000]
for i in range(3):
    E = x[i]
    e,p = readMCNP(E,1000)
    me[:,i] = e
#print(me[0:5,:])

with open("matrice/E_depose.txt","w") as file:
    for i in range(1002):
        for j in range(3):
            file.write("%e"%me[i][j])
            file.write('         ')
        file.write('\n')
    file.write('\n')
'''
