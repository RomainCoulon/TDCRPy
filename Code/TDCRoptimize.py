# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 10:04:53 2023

@author: romain.coulon
"""

import numpy as np
import TDCRPy as td
import scipy.optimize as opt
import sys, time
sys.path.insert(1, 'G:\Python_modules\BIPM_RI_PyModules')
import TDCRcalculation as tc 



def eff(TDCR_measure, TAB, TBC, TAC, Rad, kB):
    
    N=50000
    L=1
    # r=opt.minimize_scalar(td.TDCRPy, args=(TDCR_measure, TAB, TBC, TAC, Rad, "1", N, kB, 0.98, 1000, "res", "sym"), method='golden')#,options={'xatol': 1e-7, 'disp': True, 'maxiter':100})
    # r=opt.minimize_scalar(td.TDCRPy, args=(TDCR_measure, TAB, TBC, TAC, Rad, "1", N, kB, 0.98, 1000, "res", "sym"), method='brent')
    r=opt.minimize_scalar(td.TDCRPy, args=(TDCR_measure, TAB, TBC, TAC, Rad, "1", N, kB, 0.98, 1000, "res", "sym"), method='bounded', bounds=[0.5, 2])
    
    # r=opt.minimize(td.TDCRPy, L, args=(TDCR_measure, TAB, TBC, TAC, Rad, "1", N, kB, 0.98, 1000, "res", "sym"), method='nelder-mead',options={'xatol': 1e-7, 'disp': True, 'maxiter':100})
    L=r.x
    print(r)
    #L=(L*0.995, L*1.021, L*0.988)           # Free paramete in keV-1
    # L = (L, L, L)
    # r=opt.minimize(td.TDCRPy, L, args=(TDCR_measure, TAB, TBC, TAC, Rad, "1", N, kB, 0.98, 1000, "res", "asym"), method='nelder-mead',options={'xatol': 1e-7, 'disp': True, 'maxiter':100})
    # L=r.x
    # print(r)
    out=td.TDCRPy(L,TDCR_measure, TAB, TBC, TAC, Rad, "1", N, kB, 0.98, 1000, "eff", "sym")
    return np.mean(L), L, out[2], out[2], out[3]
    
def TicTocGenerator():
    """
    Generator that returns time differences
    """
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    """
    Prints the time difference yielded by generator instance TicToc
    """
    
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    """
    Records a time in TicToc, marks the beginning of a time interval
    """
    toc(False)


"""
TEST
"""

Rad="Co-60"    # list of radionuclides (Na-24)
kB =1.0e-5       # Birks constant in cm/keV  
TDCR_measure = 0.977667386529166        # Measured TDCR value
TAB = 0.992232838598821
TBC = 0.992343419459002
TAC = 0.99275350064608



tic()
F1, FFF, eff1, eff2 = tc.I2calc(TDCR_measure, TAB, TBC, TAC, Rad, kB)
toc()

print("/nanalytical model")
print("free parameter = ", F1)
print("free parameters = ", FFF)
print("Double count rate efficiency (sym) = ", eff1)
print("Double count rate efficiency (asym) = ", eff2)


tic()
F1, FFF, eff1, eff2, u = eff(TDCR_measure, TAB, TBC, TAC, Rad, kB)
toc()

print("/nstochastic model")
print("free parameter = ", F1)
print("free parameters = ", FFF)
print("Double count rate efficiency (sym) = ", eff1)
print("Double count rate efficiency (asym) = ", eff2)
print("u(eff) = ", u)