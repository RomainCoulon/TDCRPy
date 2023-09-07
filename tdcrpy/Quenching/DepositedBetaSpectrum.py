import tdcrpy
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import random

# rad=[H-3, C-14, S-35, Ca-45, Ni-63, Sr-89, Sr-90, Tc-99, Pm-147, Pu-241]
rad="Cs-137"


e, p0 = tdcrpy.TDCR_model_lib.readBetaShape(rad, 'beta-', 'tot')
n=len(e)
p=[]

for k, ek in enumerate(e):
    if k==0:
        p.append(p0[k] * (e[k+1])-ek)
    else:
        p.append(p0[k] * (ek-e[k-1]))
    # print(k,n)
    # if k==n-1:
    #     p.append(p0[k] * (ek-e[k-1]))
    # else:
    #     p.append(p0[k] * (e[k+1]-ek))
    
p/=sum(p)

N=10000

ei=[]
ed=[]
for i in range(N):
    # ind = tdcrpy.TDCR_model_lib.sampling(p)
    ind = random.choices(range(len(p)), weights=p)[0]
    ei.append(e[ind])
    ed.append(tdcrpy.TDCR_model_lib.energie_dep_beta(ei[-1]))
    # print(ind,ed[-1])
    
# kde = gaussian_kde(ed, bw_method='scott')
# kde = gaussian_kde(ed, bw_method='silverman')
kde = gaussian_kde(ei, bw_method=0.02)
kde2 = gaussian_kde(ed, bw_method=0.02)
# e2 = e
e2 = np.linspace(min(ei), max(ei), n)


p2 = kde(e2)
p2 /= sum(p2)




p3 =kde2(e2)
p3 /= sum(p3)
# p4=[]
# indexes = [index for index, _ in sorted(enumerate(ed), key=lambda x: x[1])]
# for i in indexes:
#     p4.append(p3[i])
# ed.sort() 




# ed2 = []
# p2 = []
# for y in range(N):
#     ed = []
#     for i, ei in enumerate(e):
#         # calcul of the deposited energy
#         ed.append(tdcrpy.TDCR_model_lib.energie_dep_beta(ei))
    
#     indexes = [index for index, _ in sorted(enumerate(ed), key=lambda x: x[1])]
#     # print(indexes)

#     for i in indexes:
#         p2.append(p[i])
#     ed.sort()    
#     ed2 += ed

# p3 = savgol_filter(p2, 11, 2)


plt.figure("Spectra")
plt.clf()
plt.plot(e,p0,"-k",label="BetaShape")
plt.plot(e,p,"--b",label="BetaShape'")
plt.plot(e2,p2,"-g",label="re-sample")
# plt.plot(e2,p3,"-r",label="MCNP6")
plt.legend()
plt.xlabel("Energy /keV")
plt.ylabel("Probability")
plt.show()
