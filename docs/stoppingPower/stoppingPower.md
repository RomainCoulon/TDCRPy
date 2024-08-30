# Estimation of stopping powers in liquid scintillator


```python
# pip install TDCRPy --upgrade
```


```python
import tdcrpy as td
import matplotlib.pyplot as plt
import numpy as np
```

## Stopping power at a given energy


```python
density = 0.96 # in g/cm3
energy = 100 # keV
result_alpha = td.TDCR_model_lib.stoppingpowerA(energy,rho=density)
result_electron = td.TDCR_model_lib.stoppingpower(energy*1e3,rho=density)
print(f"dE/dx = {result_alpha:.7g} keV/cm for alpha particles" )
print(f"dE/dx = {result_electron*1e3:.4g} keV/cm for electrons" )
```

    dE/dx = 1461120 keV/cm for alpha particles
    dE/dx = 3419 keV/cm for electrons
    

## Plot stopping power curves


```python
energy_vec = np.logspace(-2,4,1000) # in keV
w_a, w_e = [], []
for e in energy_vec:
        w_a.append(td.TDCR_model_lib.stoppingpowerA(e,rho=density))
        w_e.append(1e3*td.TDCR_model_lib.stoppingpower(e*1e3,rho=density))
plt.figure("Stopping power")
plt.clf()
plt.plot(energy_vec,w_a,"-k",label="alpha particles")
plt.plot(energy_vec,w_e,"-r",label="electrons")
plt.xscale('log')
plt.xlabel(r'$E$ /keV', fontsize=14)
plt.ylabel(r'd$E$/d$x$ /(keV/cm)', fontsize=14)
plt.legend()
```




    <matplotlib.legend.Legend at 0x25d905a4c90>




    
![png](output_6_1.png)
    

