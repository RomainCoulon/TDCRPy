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

e,p = readMCNP(4,1000)
print(sum(p))

'''
plt.plot(e,p)
plt.xlabel('Energy/MeV')
plt.ylabel('Probability')
plt.yscale('log')
plt.title('spectre gamma at E=1MeV')
plt.savefig('proba = f(E)_1MeV_1e6.png')
'''