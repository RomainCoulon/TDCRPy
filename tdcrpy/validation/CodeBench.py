# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 10:02:04 2024

@author: romain.coulon
"""
# import tdcrpy as td
import numpy as np
import matplotlib.pyplot as plt

rad ="Fe-55"
# kB = [6, 10, 15]
kB = [10]
CIEMATCODE = "NUR"
# CIEMATCODE = "EFFY"

BIPMCODE = "TDCRPy"
# BIPMCODE = "TDCRPyA"

TDCR17 = True

# Fe-55
# ExpTD = [5.015440131719763883e-01, 3.033921221898202569e-01, 2.337696889928395638e-01, 1.963266220117187155e-01, 1.525520789152907120e-01, 1.290251872314254478e-01]
# ExpTD = [5.036195861139387553e-01, 3.068002450941523795e-01, 2.347349238197918575e-01, 1.868167474073787748e-01, 1.492620432358786609e-01, 1.111062699027455108e-01]
# ExpTD = [5.129091648648070878e-01,3.151410221276935708e-01,2.473618561405817318e-01,2.038446684986770507e-01,1.555745292587844830e-01,1.225838288948343929e-01]
Dtdcr17 = [0.7954, .7925, .7075, .5999, .7760, .7757]

xx = 100*(1-2*5.129091648648070878e-01)


# Co-60
# ExpTD = [9.818128981512568298e-01, 9.714731781760674867e-01, 9.645707367051945536e-01, 9.554547633467871393e-01, 9.488064225575219002e-01, 9.397242455346862533e-01]

# Cd-109
# ExpTD = [8.522594641728009623e-01, 8.995376644905349606e-01, 9.175877883182167460e-01, 9.242570160461722750e-01, 9.220671007346888937e-01, 8.994349939851205011e-01]

# C-14
# ExpTD = [0.96512047, 0.94196669, 0.92484093, 0.90833458, 0.88807086, 0.850513  ]

# H-3
ExpTD = [0.77784194, 0.77443954, 0.67543926, 0.55646075, 0.75506601, 0.7546677]
Dtdcr17 = [0.7954, .7925, .7075, .5999, .7760, .7757]

if CIEMATCODE == "EFFY": endHeader=72
elif CIEMATCODE == "NUR": endHeader=4

nd = ['d = 0', 'd = 0.1', 'd = 0.2', 'd = 0.3', 'd = 0.4', 'd = 0.5' ]

for kBn in kB:
    """
    NUR DATA
    """
    folderNUR = rad+CIEMATCODE+"/"
    
    with open(folderNUR+rad+'_'+str(1)+'_'+str(kBn)+'.TAB', 'r') as file:
        # Read all lines into a list
        lines = file.readlines()
    
    # Print the content of the file
    lNUR = []; tNUR = []; dNUR = []
    for iline, line in enumerate(lines):
        if iline == 2:
            header=line.strip().split()
        if line == "\n" and iline>endHeader:
            break
        if "--------" in line and iline>endHeader:
            break
        if iline >= endHeader:
            A=line.split()
            lNUR.append(1/float(A[0]))
            if CIEMATCODE=="EFFY":
                tNUR.append(float(A[4])/100)
                dNUR.append(float(A[5])/100)
            else:
                tNUR.append(float(A[4]))
                dNUR.append(float(A[5]))                
    tNUR = np.asarray(tNUR); dNUR = np.asarray(dNUR)
            
            
            
    """
    TDCRPy DATA
    """
    folderTDCRPy = rad+BIPMCODE+"/"
    
    with open(folderTDCRPy+"result_A_"+rad+'_'+str(kBn)+'.txt', 'r') as file:
        # Read all lines into a list
        lines = file.readlines()        
    
    # Print the content of the file
    lTDCRPy = []; tTDCRPy = []; dTDCRPy = []; tTDCRPys = []; dTDCRPys = []
    for iline, line in enumerate(lines):
        A=line.split()
        lTDCRPy.append(float(A[0]))
        dTDCRPy.append(float(A[3]))
        tTDCRPy.append(float(A[5]))
        dTDCRPys.append(float(A[4]))
        tTDCRPys.append(float(A[6]))    
    
    dTDCRPy = np.asarray(dTDCRPy)
    tTDCRPy = np.asarray(tTDCRPy)
    
    """
    DISPLAY
    """
    
    
    # plt.figure("T vs L")
    # plt.clf()
    # plt.plot(lNUR,tNUR,label="NUR")
    # plt.errorbar(lTDCRPy,tTDCRPy,yerr=tTDCRPys,fmt="r", label="TDCRPy")
    # plt.legend()
    # plt.xlabel(r"$l$ /keV-1")
    
    tdNUR  = tNUR/dNUR
    tdTDCRPy = tTDCRPy/dTDCRPy
       
    dev = []
    tdvec= []
    indm = 0
    for u, utd in enumerate(tdTDCRPy):
        ind = abs(tdNUR-utd).argmin()
        if indm != ind:
            tdvec.append(utd)
            dev.append(100*(dTDCRPy[u]/dNUR[ind]-1))
            indm=ind
    
    tdexp =[]
    dexp =[]
    count = 0
    for v in ExpTD:
        if rad == "Cd-109":
            index=abs(np.asarray(tdTDCRPy)-v).argmin()
            if count == 0:
                index=abs(np.asarray(tdTDCRPy)-v).argmin()
                index0=index
                count+=1
            else:
                count+=1
                index=abs(np.asarray(tdTDCRPy[:index0])-v).argmin()
                index0=index
            print(index0,tdTDCRPy[index0])
        else:
            index=abs(np.asarray(tdTDCRPy)-v).argmin()
        tdexp.append(tdTDCRPy[index])
        dexp.append(dTDCRPy[index])
    
    plt.figure("D vs TDCR", figsize=(8,6))
    plt.clf()
    if rad == "Cd-109":
        plt.plot(tdNUR[:-900],dNUR[:-900],'-b',label="NUR")
    elif rad == "Co-60":
        plt.plot(tdNUR[:-900],dNUR[:-900],'-b',label="NUR")
    else:
        plt.plot(tdNUR,dNUR,'-b',label=CIEMATCODE)
    plt.errorbar(tdTDCRPy,dTDCRPy,yerr=dTDCRPys,fmt="-*k", label=BIPMCODE)
    plt.plot(tdexp,dexp,"or", label="Exp")
    if TDCR17:
        plt.plot(tdexp,Dtdcr17,"og", label="TDCR17")
    plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    for itxt, xtxt in enumerate(tdexp):
        plt.text(xtxt+0.005,dexp[itxt]-0.00,nd[itxt],fontstyle='italic')
    plt.ylabel(r"$\epsilon_{D}$",fontsize=16)
    plt.xlabel(r"$\epsilon_{T}$/$\epsilon_{D}$",fontsize=16)
    plt.show()
    plt.savefig(rad+'.png', dpi=100)
    
        
    plt.figure('deviation')
    plt.clf()
    plt.plot(tdvec,dev,'-b')
    plt.plot(tdTDCRPy,100*np.asarray(dTDCRPys)/dTDCRPy,'-r')
    # plt.plot(tdexp,dexp,'or')
    plt.ylabel(r"$d$ (%)")
    plt.xlabel(r"TDCR")    
            
        
