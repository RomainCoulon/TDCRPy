# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:45:10 2024

@author: romain.coulon
"""

import tdcrpy as td
import numpy as np
import matplotlib.pyplot as plt
import cProfile

"""
Test efficiency
"""

# def testEff():
td.TDCRPy.TDCRPy(1, 0.97, 0.97, 0.97, 0.97, "Co-60", "1", 10000, 1e-5, 16, "eff", "sym", Smodel=False)
    
# cProfile(testEff())