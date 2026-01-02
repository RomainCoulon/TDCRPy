# TDCRPy

<div align="center">

<img src="docs/logo1.png" alt="TDCRPy Logo" width="200"/>

**A Photo-Physical Stochastic Model for Liquid Scintillation Counting**

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![Status](https://img.shields.io/badge/status-stable-green)
![BIPM](https://img.shields.io/badge/maintained%20by-BIPM-005696)

</div>

---

## 📖 Overview

**TDCRPy** is a Python package developed and maintained by the **BIPM** (Bureau International des Poids et Mesures). It estimates detection efficiencies of liquid scintillation counters using the **TDCR** (Triple to Double Coincidence Ratio) or **CIEMAT/NIST** methods.

The calculation is based on a photo-physical stochastic model, allowing users to address:
* Complex decay schemes (Beta spectra via BetaShape).
* Radionuclide mixtures.
* Ionization quenching (Birks model).
* Micelle effects in scintillator cocktails.
* Dynamic efficiency evolution over time.

Technical details can be found in:
* [Coulon et al., Applied Radiation and Isotopes (2024)](https://doi.org/10.1016/j.apradiso.2024.111518)
* [Coulon et al., BIPM Technical Report](http://dx.doi.org/10.13140/RG.2.2.15682.80321)

---

## 📦 Installation

TDCRPy requires a standard Python scientific environment.

### 1. Install Dependencies
You can install the required libraries via `pip` or `conda`.

```shell
# Using pip
pip install numpy scipy configparser tqdm importlib-resources

# Optional (for visualization and dynamic decay features)
pip install opencv-python radioactivedecay matplotlib
````

### 2\. Install TDCRPy

```shell
pip install TDCRPy
```

To upgrade to the latest version:

```shell
pip install TDCRPy --upgrade
```

### 3\. Run Tests

Verify the installation by running the unit tests:

```shell
python -m unittest tdcrpy.test.test_tdcrpy
```

-----

## ⚡ Quick Start

Here is a basic example to estimate detection efficiencies for **Co-60** using a symmetric PMT configuration.

```python
import tdcrpy

# --- 1. Define Parameters ---
mode = "eff"           # Calculation mode
L = 1.2                # Free parameter in keV^-1
Rad = "Co-60"          # Radionuclide
pmf_1 = "1"            # Relative fraction (100%)
N = 1000               # Number of Monte Carlo trials
kB = 1.0e-5            # Birks constant in cm keV^-1
V = 10                 # Volume of scintillator in mL

# --- 2. Run Calculation ---
result = tdcrpy.TDCRPy.TDCRPy(L, Rad, pmf_1, N, kB, V, mode)

# --- 3. Display Results ---
print(f"Efficiency S (Single): {result[0]:.4f} +/- {result[1]:.4f}")
print(f"Efficiency D (Double): {result[2]:.4f} +/- {result[3]:.4f}")
print(f"Efficiency T (Triple): {result[4]:.4f} +/- {result[5]:.4f}")
print(f"Efficiency D (CIEMAT/NIST): {result[12]:.4f} +/- {result[13]:.4f}")
```

### Calculate Free Parameter from Measured TDCR

If you have an experimental TDCR value ($R_T/R_D$), you can reverse-calculate the free parameter $L$:

```python
TD = 0.9776  # Measured TDCR parameter
result = tdcrpy.TDCRPy.eff(TD, Rad, pmf_1, kB, V)

print(f"Global free parameter L = {result[0]:.4f} keV^-1")
```

-----

## 🛠 Advanced Features

### Asymmetric PMTs

TDCRPy supports calculations where the quantum efficiency differs between PMTs. Pass a tuple for the free parameter $L$:

```python
# L for (PMT A, PMT B, PMT C)
L = (1.1, 1.3, 1.2) 
result = tdcrpy.TDCRPy.TDCRPy(L, "Co-60", "1", 1000, 1.0e-5, 10)

# Efficiencies for specific pairs (AB, BC, AC) are available in the result tuple
print(f"Efficiency AB: {result[6]:.4f}")
```

### Radionuclide Mixtures

Simulate mixtures by providing comma-separated nuclides and their relative fractions.

```python
# Example: 80% Co-60 and 20% H-3
Rad = "Co-60, H-3"
pmf_1 = "0.8, 0.2" 

result = tdcrpy.TDCRPy.TDCRPy(L, Rad, pmf_1, N, kB, V)
```

### Dynamic Decay & Efficiency

Combine `TDCRPy` with `radioactivedecay` to simulate how efficiency changes as a sample decays (e.g., Mo-99/Tc-99m).

```python
import radioactivedecay as rd
import tdcrpy as td
import numpy as np

# Define inventory
rad_t0 = rd.Inventory({'Mo-99': 1}, 'Bq')

# Decay for 30 hours
rad_t1 = rad_t0.decay(30.0, 'h')

# Calculate current composition for TDCRPy
A_t1 = rad_t1.activities('Bq')
total_activity = sum(A_t1.values())

# Format strings for TDCRPy
nuclides = ", ".join([k for k, v in A_t1.items() if v > 0])
fractions = ", ".join([str(v/total_activity) for k, v in A_t1.items() if v > 0])

# Run simulation
result = td.TDCRPy.TDCRPy(1.0, nuclides, fractions, 1000, 1e-5, 10, "eff")
```

-----

## ⚙️ Configuration & Physics

You can customize the underlying physics model using `tdcrpy.TDCR_model_lib`.

To view current settings:

```python
import tdcrpy as td
td.TDCR_model_lib.readParameters(disp=True)
```

### Common Configurations

| Parameter | Method | Description |
| :--- | :--- | :--- |
| **Electron Bins** | `modifynE_electron(n)` | Integration bins for electrons (default: 1000) |
| **Alpha Bins** | `modifynE_alpha(n)` | Integration bins for alpha particles (default: 1000) |
| **Density** | `modifyDensity(rho)` | Scintillator density in g/cm³ (default: 0.96) |
| **Mean Z / A** | `modifyZ(z)`, `modifyA(a)` | Mean atomic/mass numbers of the cocktail |
| **Micelle Effect** | `modifyMicCorr(bool)` | Activate reverse micelle correction (default: False) |
| **Micelle Size** | `modifyDiam_micelle(d)` | Diameter in nm (default: 2.0) |
| **Dead Time** | `modifyDeadTime(t)` | Extended dead time in µs (default: 10) |
| **Coincidence Time** | `modifyTau(ns)` | Resolving time in ns (default: 50) |

## Tutorials

* [Use the stochastic model for all nuclides](https://github.com/RomainCoulon/TDCRPy/blob/main/notebooks/tuturial.ipynb)
* [Use the analytical model for pure beta emitters](https://github.com/RomainCoulon/TDCRPy/blob/main/notebooks/analyticalModel.ipynb)
* [Modify the parameters of the model](https://github.com/RomainCoulon/TDCRPy/blob/main/notebooks/changeParameters.ipynb)

-----

## 📚 Citations

If you use **TDCRPy** in your work, please cite the following:

> **TDCRPy: A python package for TDCR measurements** \> R. Coulon, J. Hu  
> *Applied Radiation and Isotopes* (2024)  
> DOI: [10.1016/j.apradiso.2024.111518](https://doi.org/10.1016/j.apradiso.2024.111518)

-----

## ⚖️ License

This project is licensed under the **MIT License**.

Copyright © BIPM (Bureau International des Poids et Mesures).

