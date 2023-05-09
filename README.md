# TDCRPy

TDCRPy is a Python code to calculate detection efficiency of a liquide scintillation counter using 3-photomultiplier tubes.
The calculation is based on the photo-physical model called of the Triple-to-Double-Coincidence-Ratio method (TDCR) [[1]](#1) and a Monte-Carlo sampling allowing to adress complexe decay schemes and radionuclide mixtures.

The code directly reads decay data from the Decay Data Evaluation Project (DDEP) web interface [[2]](#2) that is recommanded to be used by the radionuclide metrology community. The PenNuc format [[3]](#3) is used to simulate decays and the $\beta$ spectra from the BetaShape code <c id="4">[4]</c> <c id="5">[5]</c> are used. The BetaShape code estimates accurate $\beta$ spectra by taking the atomic exchange effect. It has been demonstrated to offer significant improvement in the context of liquid scintillation counting [[6]](#6).

The stopping power of electrons between 20 keV and 1000 keV is a mixture of a radiative loss model and a collision model that has been validated agaisnt the NIST model ESTAR [[7]](#7) recommanded by the ICRU Report 37. At low energy - between 10 eV and 20 keV - the model from Tan and Xia [[9]](#9) is implemented.

The stopping power of $\alpha$ particles of energy comprises between 1 keV and 8 MeV comes from the NIST code ASTAR [[7]](#7) recommanded in the ICRU Report 49 [[10]](#10). For energy below 1 keV, an extrapolation is made.


<a id="1">[1]</a> https://doi.org/10.1016/0883-2889(88)90161-X

<b id="2">[2]</b> http://www.lnhb.fr/ddep_wg/

<c id="3">[3]</c> https://doi.org/10.1016/j.cpc.2019.08.002

<c id="4">[4]</c> https://doi.org/10.1103/PhysRevC.92.059902

<c id="5">[5]</c> https://doi.org/10.1016/j.apradiso.2019.108884

<c id="6">[6]</c> https://doi.org/10.1016/j.apradiso.2020.109478

<c id="7">[7]</c> https://dx.doi.org/10.18434/T4NC7P

<c id="8">[8]</c> ICRU Report 37, Stopping Powers for Electrons and Positrons

<c id="9">[9]</c> https://doi.org/10.1016/j.apradiso.2011.08.012

<c id="10">[10]</c> https://www.icru.org/report/stopping-power-and-ranges-for-protons-and-alpha-particles-report-49/
