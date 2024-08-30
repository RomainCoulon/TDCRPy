<a href="https://colab.research.google.com/github/RomainCoulon/TDCRPy/blob/main/Tuturial.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# Quick start with the TDCRPy code


```python
#pip install TDCRPy --upgrade
```


```python
import tdcrpy
```

## 1. Symmetric PMTs

### 1.1 Estimation of detection efficiencies given the free parameter


```python
mode = "eff"               # ask for efficiency calculation
L = 1.2                    # free parameter in keV-1
Rad="Co-60"                # radionuclide
pmf_1="1"                  # relative fraction of the radionulide
N = 1000                   # number of Monte Carlo trials
kB =1.0e-5                 # Birks constant in cm keV-1
V = 10                     # volume of scintillator in mL
```


```python
result = tdcrpy.TDCRPy.TDCRPy(L, Rad, pmf_1, N, kB, V, mode)
```


```python
print(f"efficiency S = {round(result[0],4)} +/- {round(result[1],4)}")
print(f"efficiency D = {round(result[2],4)} +/- {round(result[3],4)}")
print(f"efficiency T = {round(result[4],4)} +/- {round(result[5],4)}")
print(f"efficiency D (C/N system) = {round(result[12],4)} +/- {round(result[13],4)}")
```

    efficiency S = 0.9859 +/- 0.0031
    efficiency D = 0.9738 +/- 0.0043
    efficiency T = 0.9527 +/- 0.0057
    efficiency D (C/N system) = 0.9702 +/- 0.0045
    

### 1.2 Estimation of detection efficiencies given the measured TDCR parameter


```python
TD = 0.977667386529166     # TDCR parameter
```


```python
result = tdcrpy.TDCRPy.eff(TD, Rad, pmf_1, kB, V)
```

    global free parameter = 1.1991314352332345 keV-1
    


```python
print(f"free parameter = {round(result[0],4)} keV-1")
print(f"efficiency S = {round(result[2],4)} +/- {round(result[3],4)}")
print(f"efficiency D = {round(result[4],4)} +/- {round(result[5],4)}")
print(f"efficiency T = {round(result[6],4)} +/- {round(result[7],4)}")
print(f"efficiency D (C/N system) = {round(result[14],4)} +/- {round(result[15],4)}")
```

    free parameter = 1.1991 keV-1
    efficiency S = 0.9858 +/- 0.001
    efficiency D = 0.973 +/- 0.0014
    efficiency T = 0.9512 +/- 0.0018
    efficiency D (C/N system) = 0.9692 +/- 0.0014
    

## 2. Asymmetric PMTs

### 2.1 Estimation of detection efficiencies given the free parameters


```python
mode = "eff"               # ask for efficiency calculation
L = (1.1, 1.3, 1.2)        # free parameter in keV-1
Rad="Co-60"                # radionuclide
pmf_1="1"                  # relative fraction of the radionulide
N = 1000                   # number of Monte Carlo trials
kB =1.0e-5                 # Birks constant in cm keV-1
V = 10                     # volume of scintillator in mL
```


```python
result = tdcrpy.TDCRPy.TDCRPy(L, Rad, pmf_1, N, kB, V)
```


```python
print(f"efficiency S = {round(result[0],4)} +/- {round(result[1],4)}")
print(f"efficiency D = {round(result[2],4)} +/- {round(result[3],4)}")
print(f"efficiency T = {round(result[4],4)} +/- {round(result[5],4)}")
print(f"efficiency AB = {round(result[6],4)} +/- {round(result[7],4)}")
print(f"efficiency BC = {round(result[8],4)} +/- {round(result[9],4)}")
print(f"efficiency AC = {round(result[10],4)} +/- {round(result[11],4)}")
print(f"efficiency D (C/N system) = {round(result[12],4)} +/- {round(result[13],4)}")
```

    efficiency S = 0.9904 +/- 0.0024
    efficiency D = 0.9803 +/- 0.0037
    efficiency T = 0.9625 +/- 0.005
    efficiency AB = 0.9684 +/- 0.0045
    efficiency BC = 0.9696 +/- 0.0044
    efficiency AC = 0.9674 +/- 0.0046
    efficiency D (C/N system) = 0.9772 +/- 0.0039
    

### 2.2 Estimation of detection efficiencies given the measured TDCR parameters


```python
TD = [0.977667386529166, 0.992232838598821, 0.992343419459002, 0.99275350064608]
```


```python
result = tdcrpy.TDCRPy.eff(TD, Rad, pmf_1, kB, V)
```

    global free parameter = 1.2369501044018718 keV-1
    free parameters = [1.2369501 1.2369501 1.2369501] keV-1
    


```python
print(f"Global free parameter = {round(result[0],4)} keV-1")
print(f"free parameter of PMT A = {round(result[1][0],4)} keV-1")
print(f"free parameter of PMT B = {round(result[1][1],4)} keV-1")
print(f"free parameter of PMT C = {round(result[1][2],4)} keV-1")
print(f"efficiency S = {round(result[2],4)} +/- {round(result[3],4)}")
print(f"efficiency D = {round(result[4],4)} +/- {round(result[5],4)}")
print(f"efficiency T = {round(result[6],4)} +/- {round(result[7],4)}")
print(f"efficiency AB = {round(result[8],4)} +/- {round(result[9],4)}")
print(f"efficiency BC = {round(result[10],4)} +/- {round(result[11],4)}")
print(f"efficiency AC = {round(result[12],4)} +/- {round(result[13],4)}")
print(f"efficiency D (C/N system) = {round(result[14],4)} +/- {round(result[15],4)}")
```

    Global free parameter = 1.237 keV-1
    free parameter of PMT A = 1.237 keV-1
    free parameter of PMT B = 1.237 keV-1
    free parameter of PMT C = 1.237 keV-1
    efficiency S = 0.9872 +/- 0.0009
    efficiency D = 0.9751 +/- 0.0013
    efficiency T = 0.9533 +/- 0.0018
    efficiency AB = 0.9606 +/- 0.0016
    efficiency BC = 0.9606 +/- 0.0016
    efficiency AC = 0.9606 +/- 0.0016
    efficiency D (C/N system) = 0.9714 +/- 0.0014
    

## 3. Radionuclide mixture

### 3.1 Estimation of detection efficiencies given the free parameter


```python
mode = "eff"               # ask for efficiency calculation
mode2 = "sym"              # specify that symmetric PMTs is considered
L = 1.2                    # free parameter in keV-1
Rad="Co-60, H-3"           # radionuclides
pmf_1="0.8, 0.2"                  # relatives fractions of the radionulides
N = 1000                   # number of Monte Carlo trials
kB =1.0e-5                 # Birks constant in cm keV-1
V = 10                     # volume of scintillator in mL
```


```python
result = tdcrpy.TDCRPy.TDCRPy(L, Rad, pmf_1, N, kB, V)
```


```python
print(f"efficiency S = {round(result[0],4)} +/- {round(result[1],4)}")
print(f"efficiency D = {round(result[2],4)} +/- {round(result[3],4)}")
print(f"efficiency T = {round(result[4],4)} +/- {round(result[5],4)}")
print(f"efficiency D (C/N system) = {round(result[12],4)} +/- {round(result[13],4)}")
```

    efficiency S = 0.9622 +/- 0.0045
    efficiency D = 0.9163 +/- 0.007
    efficiency T = 0.8482 +/- 0.0096
    efficiency D (C/N system) = 0.9037 +/- 0.0074
    

### 3.2 Estimation of detection efficiencies given the measured TDCR parameter


```python
TD = 0.977667386529166     # TDCR parameter
```


```python
result = tdcrpy.TDCRPy.eff(TD, Rad, pmf_1, kB, V)
```

    global free parameter = 4.999573818598962 keV-1
    


```python
print(f"free parameter = {round(result[0],4)} keV-1")
print(f"efficiency S = {round(result[2],4)} +/- {round(result[3],4)}")
print(f"efficiency D = {round(result[4],4)} +/- {round(result[5],4)}")
print(f"efficiency T = {round(result[6],4)} +/- {round(result[7],4)}")
print(f"efficiency D (C/N system) = {round(result[14],4)} +/- {round(result[15],4)}")
```

    free parameter = 4.9996 keV-1
    efficiency S = 0.9835 +/- 0.001
    efficiency D = 0.9677 +/- 0.0015
    efficiency T = 0.9397 +/- 0.002
    efficiency D (C/N system) = 0.9629 +/- 0.0016
    


```python

```
