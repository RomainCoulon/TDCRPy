# Modify the parameters of the model (advanced settings)


```python
#pip install TDCRPy --upgrade
```


```python
import tdcrpy as td
```

## Read the parameters


```python
print("\nparameters in a list: ", td.TDCR_model_lib.readParameters(disp=True))
```

    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    
    parameters in a list:  (1000, 1000, 0.96, 5.2, 11.04, 5, 100.0, 1.5, 2.0, 0.1, 50, 10, 20, False)
    

### Change energy binning the quenching function of electrons 


```python
td.TDCR_model_lib.modifynE_electron(100)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifynE_electron(1000)
```

    New configuration:
    number of integration bins for electrons = 100
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change energy binning the quenching function of alpha particles


```python
td.TDCR_model_lib.modifynE_alpha(100)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifynE_alpha(1000)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 100
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change the density (in g/cm3) of the LS source


```python
td.TDCR_model_lib.modifyDensity(1.02)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyDensity(0.96)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 1.02 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change the mean charge number of the LS source


```python
td.TDCR_model_lib.modifyZ(5.7)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyZ(5.2)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.7
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change the mean atomic mass number of the LS source


```python
td.TDCR_model_lib.modifyA(12.04)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyA(11.04)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 12.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change the depht paramerter of the spline interpolation of the quenching function


```python
td.TDCR_model_lib.modifyDepthSpline(7)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyDepthSpline(5)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 7
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change the energy threshold (in keV) above which the interpolation is applied (for alpha particles)


```python
td.TDCR_model_lib.modifyEinterp_a(200)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyEinterp_a(100)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 200.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change the energy threshold (in keV) above which the interpolation is applied (for electrons)


```python
td.TDCR_model_lib.modifyEinterp_e(2.0)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyEinterp_e(1.5)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 2.0 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Activate/desactivate the micelles correction


```python
td.TDCR_model_lib.modifyMicCorr(False)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyMicCorr(True)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change the diameter of reverse micelles (in nm)


```python
td.TDCR_model_lib.modifyDiam_micelle(4.0)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyDiam_micelle(2.0)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 4.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change the acqueous fraction of the LS source


```python
td.TDCR_model_lib.modifyfAq(0.2)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyfAq(0.1)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.2
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change the coincidence resolving time of the TDCR counter (in ns)


```python
td.TDCR_model_lib.modifyTau(100)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyTau(50)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 100 ns
    extended dead time = 10 µs
    measurement time = 20 min
    

### Change the extended dead time of the TDCR counter (in µs)¶


```python
td.TDCR_model_lib.modifyDeadTime(100)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyDeadTime(10)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 100 µs
    measurement time = 20 min
    

### Change the measurement time (in min)¶


```python
td.TDCR_model_lib.modifyMeasTime(60)
print("New configuration:")
td.TDCR_model_lib.readParameters(disp=True)
# back to the default value
td.TDCR_model_lib.modifyMeasTime(20)
```

    New configuration:
    number of integration bins for electrons = 1000
    number of integration bins for alpha = 1000
    density = 0.96 g/cm3
    Z = 5.2
    A = 11.04
    depth of spline interp. = 5
    energy above which interp. in implemented (for alpha) = 100.0 keV
    energy above which interp. in implemented (for electron) = 1.5 keV
    activation of the micelle correction = False
    diameter of micelle = 2.0 nm
    acqueous fraction = 0.1
    coincidence resolving time = 50 ns
    extended dead time = 10 µs
    measurement time = 60 min
    
