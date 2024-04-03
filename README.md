# TDCRPy

`TDCRPy` is a Python code to calculate detection efficiency of a liquid scintillation counter using 3-photomultiplier tubes.
The calculation is based on the photo-physical model called of the Triple-to-Double-Coincidence-Ratio method (TDCR) and a Monte-Carlo sampling allowing to adress complexe decay schemes and radionuclide mixtures.

The code is developped and maintained by the BIPM (MIT license).

## Installation

`TDCRPy` requires that the following packages are installed in your `Python` environement.

```shell
pip install importlib.resources configparser numpy tqdm setuptools scipy
```
or in `conda` environement:

```shell
conda install importlib.resources configparser numpy tqdm setuptools scipy
```

Then, `TDCRPy` can be installed.

```shell
pip install TDCRPy
```

To obtain the last version.

```shell
pip install TDCRPy --upgrade
```

The module can be imported in your Python code such as.

```python
import tdcrpy
```

## Test

To run the unit tests of the package:

```shell
python -m unittest tdcrpy.test.test_tdcrpy
```

or (using the coverage package)

```shell
coverage run -m unittest tdcrpy.test.test_tdcrpy
coverage report -m
```

## User guide

At first, the TDCRPy module must be imported.

```python
import tdcrpy
```
### TDCRPy()

The main function of the TDCRPy module is `TDCRPy()` in which the Monte-Carlo Triple-to-Double Coincidence Ratio model is implemented. The computation is made for a given solution containing a radionuclide (or a mixture of radionuclides) `rad`, a given volume of scintillator `V` and a given Birks constant `kB`.

It can operates in two modes:
*  In `mode="eff"`, it calculates the efficiency of the TDCR system as a function of a value (or a triplet) of the free parameter(s) `L`. The measurement data is not used is this mode.
*  In `mode="res"`, it calculates the residual of the TDCR model parametrized by a value (or triplet) of the free parameter(s) `L` and the measurement data `TD`, `TAB`, `TBC`, `TAC`. This residual is used in the optimisation procedure `tdcrpy.TDCRoptimize.eff()` that aims to estimate the dection efficiency of the TDCR system (see above).

Also, two configurations can be set:

* In `mode2="sym"`, the symmetry is considered between the 3 photomultiplier tubes (PMTs). Here, the free parameter `L` is a scalar and only the global TDCR value `TD` is used as measurement data.
* In `mode2="asym"`, a possible asymmetry between the 3 photomultiplier tubes is considered in the model. Here, the free parameter `L` is a triplet corresponding to each PMTs and only the specific TDCR values `TAB`, `TBC` and `TAC` are used as measurement data.

The parmeter `N` sets the number of Monte-Carlo trails used for the estimation. Each MC trial corresponds to a simulated radiactive decay.

TDCRPy() is parametrised by the following parameters:
*  `L` is the free parameter(s) in keV^-1^. It represents the number of photoelectrons produced by the PMTs per unit of ionisation energy released in the scintillator during a decay without considering the scintillation quenching. It is of type `float` when `mode2="sym"` and of type `tuple` or `list` when `mode2="asym"`.
*  `TD` (type `float`) is the quotient of the triple coincidence count rate divided by the logic sum of the double coincidence count rate.
*  `TAB` (type `float`) is the quotient of the triple coincidence count rate divided by the double coincidence count rate between the channel A and B.
*  `TBC` (type `float`) is the quotient of the triple coincidence count rate divided by the double coincidence count rate between the channel B and C.
*  `TAC` (type `float`) is the quotient of the triple coincidence count rate divided by the double coincidence count rate between the channel A and C.
*  `rad` (type `string`) is the list of radionuclides contained in the solution, e.g. `"H-3"`, `"H-3, Co-60"`.
*  `pmf_1` (type `string`) is the probability of each radionuclide contained in `rad`. If `rad = "H-3, Co-60"` and `pmf_1="0.8, 0.2"`, the solution contains 80 % of tritium and 20 % of cobalt 60.
*  `N` (type `integer`) is the number of Monte-Carlo trials corresponding each to a simulated nuclear decay in the scintillator. Monte-Carlo calculation is not used in the case of pure beta emitting radionuclides where the analytical model is implemented.
*  `kB` (type `float`) is the Birks constant caracterizing the scintillator. It is espressed in cm/keV.
*  `V` (type `float`) is the volume in ml of the scintillator.
*  `mode` sets the type of output: `mode="res"` to return the residual, `mode="eff"` to return efficiencies.
*  `mode2` sets whether the TDCR system has to be considered as symetrical by the model `mode2="sym"` or not symetrical `mode2="asym"`.
*  `Display` (type `boolean`) is an optional parameter set by default, `Display=False`. If `Display=True`, detailed on the simulated decays are displayed.
*  `barp` (type `boolean`) is an optional parameter set by default, `barp=False`. If `barp=True`, the progression bar of the calculation is displayed.


In  `mode="eff"`, `TDCRPy()` returns a `tuple` composed of:
*  the estimation of the counting efficiency of single events *S*,
*  the standard uncertainty of the estimation of *S*,
*  the estimation of the counting efficiency of double coincidences *D*,
*  the standard uncertainty of the estimation of *D*,
*  the estimation of the counting efficiency of triple coincidences *T*,
*  the standard uncertainty of the estimation of *T*,

Example: 

```python
import tdcrpy

L = 1.5
TD = 0.977667386529166
TAB = 0.992232838598821
TBC = 0.992343419459002
TAC = 0.99275350064608
Rad="Co-60"
pmf_1="1"
N = 100
kB =1.0e-5
V = 10
mode = "eff"
mode2 = "sym"

result = tdcrpy.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V, mode, mode2)

print("efficiency S = ", round(result[0],4), "+/-", round(result[1],4))
print("efficiency D = ", round(result[2],4), "+/-", round(result[3],4))
print("efficiency T = ", round(result[4],4), "+/-", round(result[5],4))
```

```Console
efficiency S =  0.9875 +/- 0.007
efficiency D =  0.9903 +/- 0.0072
efficiency T =  0.9749 +/- 0.0121
```

In  `mode="res"`, TDCRPy() returns the residuals *R* (type `float`) of the TDCR model for the given value of the free parameter `L`.

Example: 

```python
import tdcrpy

L = 1.5
TD = 0.977667386529166
TAB = 0.992232838598821
TBC = 0.992343419459002
TAC = 0.99275350064608
Rad="Co-60"
pmf_1="1"
N = 100
kB =1.0e-5
V = 10
mode = "res"
mode2 = "sym"

result = tdcrpy.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V, mode, mode2)

print("R = ",result)
```

```Console
R =  9.72238771394384e-06
```

When a triplet of free parameters is considered, `L` is of type `tuple` and `mode2 = "asym"`. It must be noted that the calculation time is three times higher than in the symetrical model. 

Example: 

```python
import tdcrpy

L = (1.5, 1.2, 1.4)
TD = 0.977667386529166
TAB = 0.992232838598821
TBC = 0.992343419459002
TAC = 0.99275350064608
Rad="Co-60"
pmf_1="1"
N = 100
kB =1.0e-5
V = 10
mode = "res"
mode2 = "asym"

result = tdcrpy.TDCRPy.TDCRPy(L, TD, TAB, TBC, TAC, Rad, pmf_1, N, kB, V, mode, mode2)

print("R = ",result)
```

```Console
R =  1.4353680366705997e-05
```

By default, the Monte-Carlo model is not applied for the following list radionuclides when used solely (not in a mixture): ^3^H, ^14^C, ^35^S, ^45^Ca, ^63^Ni, ^89^Sr, ^90^Sr, ^99^Tc, ^147^Pm, ^241^Pu. In the latter cases, the analytical model is applied and the calculation time does not dependent to the parameter `N`.

### Advance settings

An advanced setting can be configured in the `config.toml` file for functions `TDCRPy.TDCRPy()` and `TDCRoptimize.eff()`. In this file:
*  By default `Y = True` so that the analytical model is applied for solution containing only pure beta emitting radionuclides. If you would like to apply the MC calculation also for these nuclides, set `Y = False`.
*  The list of radionuclides for which the analytical model is applied is defined by default such as `radListPureBeta = H-3, C-14, S-35, Ca-45, Ni-63, Sr-89, Sr-90, Tc-99, Pm-147, Pu-241`.
*  The number of bins to discretize the linear energy space for quenching calculation for the radionuclides listed above has been set to induce an error from numerical approximation below 10^-4^. Thus the parameter `nE = 7000, 1000, 1000, 500, 2000, 50, 200, 500, 1000, 7000`.
*  In the case of Monte-Carlo calculation, the number of bins to discretize the linear energy space for quenching calculation can be adjusted. `nE_electron` and `nE_alpha` parameters for respectively electrons and alpha particles are respectiveley set by default such as `nE_electron = 1000`
and `nE_alpha = 1000`. These values ensure an error on quenched energy estimation from numerical approximation below 10^-3^.
* By default the calculation is set for Ultima-Gold cocktail mixed with a small amount of aqueous solution. You can adapt for a specific scintillator by changing the `density` (default `density=0.96`), the mean charge number `Z` (default `Z=5.2`) and the mean mass number `A` (default `A=11.04`) of the scintillator.
* To optimize the speed of the Monte-Carlo calculation, a spline interpolation on precalculated quenched energy is applied. The parameter `depthSpline` (default `depthSpline = 5`) sets the number of bins on each side of the energy point on which the interpolation is applied. The parameter `Einterp` (default `Einterp = 1`) set the energy (in keV) above which the interpolation is applied.

### TDCRoptimize.eff()

An optimization procedure `TDCRoptimize.eff()` is proposed to estimate the counting efficiency of a TDCR system. It approximates the free parameter(s) `L` in order to minimize the residual estimated by `TDCRPy()` for a given set of measurement data. In `mode2 = "sym"`, the *golden section search technique* is applied with a maximum number of iterations of 20. In `mode2 = "asym"`, the *Nelder-Mead technique* is applied with a maximum number of iterations of 20.

`TDCRoptimize.eff()` is parametrised by the following parameters:
*  `TD` (type `float`) is the quotient of the triple coincidence count rate divided by the logic sum of the double coincidence count rate.
*  `TAB` (type `float`) is the quotient of the triple coincidence count rate divided by the double coincidence count rate between the channel A and B.
*  `TBC` (type `float`) is the quotient of the triple coincidence count rate divided by the double coincidence count rate between the channel B and C.
*  `TAC` (type `float`) is the quotient of the triple coincidence count rate divided by the double coincidence count rate between the channel A and C.
*  `rad` (type `string`) is the list of radionuclides contained in the solution, e.g. `"H-3"`, `"H-3, Co-60"`.
*  `pmf_1` (type `string`) is the probability of each radionuclide contained in `rad`. If `rad = "H-3, Co-60"` and `pmf_1="0.8, 0.2"`, the solution contains 80 % of tritium and 20 % of cobalt 60.
*  `kB` (type `float`) is the Birks constant caracterizing the scintillator. It is espressed in cm/keV.
*  `V` (type `float`) is the volume in ml of the scintillator.
*  `mode2` sets whether the TDCR system has to be considered symetrical `mode2="sym"` or not symetrical `mode2="asym"`.
*  `N` (type `integer`) is the number of Monte-Carlo trials corresponding each, to a simulated nuclear decay in the scintillator. Monte-Carlo calculation is not used in the case of pure beta emitting radionuclides where the analytical model is implemented. It is an optional parameter set by default at `N = 10000`.
*  `L` (type `float`) initial guess of the free parameter(s) in keV^-1^. It represents the number of photoelectrons produced by the PMTs per unit of ionisation energy released in the scintillator during a decay and in the abscence of scintillation quenching. It is set by default at `L = 1`.

`TDCRoptimize.eff()` returns:
*  the global free parameter *L*,
*  the free parameters related to each channel (*L<sub>A</sub>*, *L<sub>B</sub>*, *L<sub>C</sub>*),
*  the estimation of the counting efficiency of single events *S*,
*  the standard uncertainty of the estimation of *S*,
*  the estimation of the counting efficiency of double coincidences *D*,
*  the standard uncertainty of the estimation of *D*,
*  the estimation of the counting efficiency of triple coincidences *T*,
*  the standard uncertainty of the estimation of *T*.

Example with `mode2 = "sym"`: 

```python
import tdcrpy

TD = 0.977667386529166
TAB = 0.992232838598821
TBC = 0.992343419459002
TAC = 0.99275350064608
Rad="Co-60"
pmf_1="1"
N = 250
kB =1.0e-5
V = 10
mode2 = "sym"

result = tdcrpy.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, V, mode2, N=N)

print("Global free parameter = \t", round(result[0],4), " keV-1")
print("Free parameter (PMT A) = \t", round(result[1][0],4) , " keV-1")
print("Free parameter (PMT B) = \t", round(result[1][1],4) , " keV-1")
print("Free parameter (PMT C) = \t", round(result[1][2],4) , " keV-1")
print("efficiency S = \t", round(result[2],4), "+/-", round(result[3],4))
print("efficiency D = \t", round(result[4],4), "+/-", round(result[5],4))
print("efficiency T = \t", round(result[6],4), "+/-", round(result[7],4))
```

```Console
Global free parameter = 	 1.3061  keV-1
Free parameter (PMT A) = 	 1.3061  keV-1
Free parameter (PMT B) = 	 1.3061  keV-1
Free parameter (PMT C) = 	 1.3061  keV-1
efficiency S = 	 0.979 +/- 0.0067
efficiency D = 	 0.9833 +/- 0.0064
efficiency T = 	 0.963 +/- 0.0096
```

Example with `mode2 = "asym"`: 

```python
import tdcrpy

TD = 0.977667386529166
TAB = 0.992232838598821
TBC = 0.992343419459002
TAC = 0.99275350064608
Rad="Co-60"
pmf_1="1"
N = 250
kB =1.0e-5
V = 10
mode2 = "asym"

result = tdcrpy.TDCRoptimize.eff(TD, TAB, TBC, TAC, Rad, pmf_1, kB, V, mode2, N=N)

print("Global free parameter = \t", round(result[0],4), " keV-1")
print("Free parameter (PMT A) = \t", round(result[1][0],4) , " keV-1")
print("Free parameter (PMT B) = \t", round(result[1][1],4) , " keV-1")
print("Free parameter (PMT C) = \t", round(result[1][2],4) , " keV-1")
print("efficiency S = \t", round(result[2],4), "+/-", round(result[3],4))
print("efficiency D = \t", round(result[4],4), "+/-", round(result[5],4))
print("efficiency T = \t", round(result[6],4), "+/-", round(result[7],4))
```

```Console
Global free parameter = 	 1.3061  keV-1
Free parameter (PMT A) = 	 1.3061  keV-1
Free parameter (PMT B) = 	 1.3061  keV-1
Free parameter (PMT C) = 	 1.3061  keV-1
efficiency S = 	 0.979 +/- 0.0067
efficiency D = 	 0.9833 +/- 0.0064
efficiency T = 	 0.963 +/- 0.0096
```
       
