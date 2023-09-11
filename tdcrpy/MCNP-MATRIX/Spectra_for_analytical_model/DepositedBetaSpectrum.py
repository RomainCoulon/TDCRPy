import tdcrpy
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import savgol_filter
import numpy as np
import random

radv=["H-3", "C-14", "S-35", "Ca-45", "Ni-63", "Sr-89", "Sr-90", "Tc-99", "Pm-147", "Pu-241"]

for rad in radv:
    e, p = tdcrpy.TDCR_model_lib.readBetaShape(rad, 'beta-', 'tot')
    n=len(e)
    print(n)
    N=1000000
    # N=1000
    
    ei=[]
    ed=[]
    for i in range(N):
        ind = random.choices(range(len(p)), weights=p)[0]
        ei.append(e[ind])
        ed.append(tdcrpy.TDCR_model_lib.energie_dep_beta(ei[-1]))
    
    kde = gaussian_kde(ei, bw_method=0.02)
    kde2 = gaussian_kde(ed, bw_method=0.02)
    e2 = np.linspace(min(ei), max(ei), 1000)
    
    
    p2 = kde(e2)
    p2 /= sum(p2)
    
    p3 =kde2(e2)
    p3 /= sum(p3)
    
    p4 = savgol_filter(p3, 201, 2)
    p4 = [x if x >= 0 else 0 for x in p4]
    p4 /=sum(p4)
    
    p5 = savgol_filter(p2, 201, 2)
    p5 = [x if x >= 0 else 0 for x in p5]
    p5 /=sum(p5)
    
    file_path = "dep_spectrum_"+str(rad)+".txt"
    with open(file_path, "w") as file:
        for item1, item2 in zip(e2, p4):
            row = f"{item1}\t{item2}\n"
            file.write(row)
    
    
    # plt.figure("Spectra")
    # plt.clf()
    # plt.plot(e,p,"--b",label="BetaShape'")
    # plt.plot(e2,p2,"-g",label="re-sample")
    # plt.plot(e2,p3,"-r",label="MCNP6")
    # plt.plot(e2,p4,"-k",label="MCNP6 smoothed")
    # plt.legend()
    # plt.xlabel("Energy /keV")
    # plt.ylabel("Probability")
    # plt.show()
    
    plt.figure("Spectra")
    plt.clf()
    plt.plot(e2,p5,"-k",label="Emitted spectrum")
    plt.plot(e2,p4,"--r",label="Deposited spectrum")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10),ncol=2,fontsize=16)
    plt.xlabel("Energy /keV", fontsize=20)
    plt.ylabel(r"Probability /keV$^{-1}$", fontsize=20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig("Spectra_"+rad+".png",format='png',bbox_inches='tight')
    plt.show()

