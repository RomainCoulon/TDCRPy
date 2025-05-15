# -*- coding: utf-8 -*-
"""
Created on Wed May 14 16:17:16 2025

@author: romain.coulon
"""
import numpy as np

m = 0.1
N = 1000000

ps0 = 1-np.exp(-m)

psi=np.random.poisson(m,N)
ps1 = sum(psi>0)/N
ups1 = np.sqrt(sum(psi>0))/N

print(ps0,ps1,ups1)
print(ps0-ps1,ups1)
print(abs(ps0-ps1)<2*ups1)