import numpy as np
import linecache as lc

def readMCNP(energy,np):
    e = []
    p = []
    f = open('output/output_'+str(int(energy))+'keV.o')
    data = f.readlines()
    f.close()

    for i in range(len(data)):
        m = i
        if data[i].find('cell  5') == 1: break


    end_point = m + 2 + np + 2
    for j in range(m+2,end_point):
        data[j] = data[j].split()
        e.append(data[j][0])
        p.append(data[j][1])
    return e,p

e,p = readMCNP(8,100)
print(p)