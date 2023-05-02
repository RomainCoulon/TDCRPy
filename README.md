# TDCRPy

TDCRPy is a Python code to calculate detection efficiency of a liquide scintillation counter using 3-photomultiplier tubes.
The calculation is based on the photo-physical model called of the Triple-to-Double-Coincidence-Ratio method (TDCR) [[1]](#1) and a Monte-Carlo sampling allowing to adress complexe decay schemes and radionuclide mixtures.

The code directly reads decay data from the Decay Data Evaluation Project (DDEP) web interface [[2]](#2) that is recommanded to be used by the radionuclide metrology community. The PenNuc format [[3]](#3) is used to simulate decays and the spectra from the BetaShape code are used.  

<a id="1">[1]</a> https://doi.org/10.1016/0883-2889(88)90161-X

<b id="2">[2]</b> http://www.lnhb.fr/ddep_wg/

<c id="2">[2]</c> https://doi.org/10.1016/j.cpc.2019.08.002
