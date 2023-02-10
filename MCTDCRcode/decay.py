# pip install radioactivedecay
import radioactivedecay as rd
import matplotlib.pyplot as plt

Ac225_t0 = rd.Inventory({'Ac-225': 1}, 'Bq')
Ac225_t1 = Ac225_t0.decay(2, 'd')
Ac225_t1.activities('Bq')
print(Ac225_t1)

nuc = rd.Nuclide('Ac-225')
decayChain = nuc.plot()
plt.savefig('decayChain.png')