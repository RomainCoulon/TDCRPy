# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 10:04:53 2023

@author: romain.coulon
"""

import numpy as np
import tdcrpy.TDCRPy as td
import scipy.optimize as opt
# import sys, time
# sys.path.insert(1, 'G:\Python_modules\BIPM_RI_PyModules')
# import TDCRcalculation as tc 

def eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, mode2, N=1000, RHO=0.98, nE=1000, L=1):
    """
    Caclulation of the efficiency of a TDCR system based on the model TDCRPy

    Parameters
    ----------
    L : TYPE
        DESCRIPTION.
    TD : TYPE
        DESCRIPTION.
    TAB : TYPE
        DESCRIPTION.
    TBC : TYPE
        DESCRIPTION.
    TAC : TYPE
        DESCRIPTION.
    Rad : TYPE
        DESCRIPTION.
    pmf_1 : TYPE
        DESCRIPTION.
    kB : TYPE
        DESCRIPTION.
    mode2 : TYPE
        DESCRIPTION.
    N : TYPE, optional
        DESCRIPTION. The default is 1000.
    RHO : TYPE, optional
        DESCRIPTION. The default is 098.
    nE : TYPE, optional
        DESCRIPTION. The default is 1000.
    L : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    """
    # Estimation of the free parameter that minimize the residuals
    r=opt.minimize_scalar(td.TDCRPy, args=(TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, RHO, nE, "res", "sym"), method='bounded', bounds=[0.5, 2])
    L=r.x
    print(r)
    
    if mode2 == "asym":
        L=(L, L, L)           # Free paramete in keV-1
        r=opt.minimize(td.TDCRPy, L, args=(TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, RHO, nE, "res", "asym"), method='nelder-mead',options={'xatol': 1e-7, 'disp': True, 'maxiter':100})
        L=r.x
        print(r)
        out=td.TDCRPy(L,TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, RHO, nE, "eff", "asym")
    else:
        out=td.TDCRPy(L,TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, RHO, nE, "eff", "sym")
    return np.mean(L), L, out[2], out[2], out[3]
    