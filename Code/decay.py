"""
This program aims to calculate the relative composition of radionuclides
decaying in chain using the Bateman model implemented in the radioactivedecay module form (MIT)

Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures
"""

# pip install radioactivedecay
import radioactivedecay as rd
import matplotlib.pyplot as plt

rad = "Ra-223" # Radionuclide
A0 = 1         # Inital activity
unit = 'Bq'    # unit for the activity
decaytime = 2  # decay time
unitt = 'd'    # unit for the time 

rad_t0 = rd.Inventory({rad: A0}, unit)
rad_t1 = rad_t0.decay(decaytime, unitt)
rad_t1.activities(unit)
print(rad_t1)

decayGraph = rad_t0.plot(decaytime, unitt, yunits=unit)
plt.savefig('decayGraph.png')

nuc = rd.Nuclide(rad)
decayChain = nuc.plot()
plt.savefig('decayChain.png')
