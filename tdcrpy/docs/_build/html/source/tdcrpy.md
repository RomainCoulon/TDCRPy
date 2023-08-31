::: {.document}
::: {.documentwrapper}
::: {.bodywrapper}
::: {.body role="main"}
::: {#tdcrpy-package .section}
# tdcrpy package[¶](#tdcrpy-package "Permalink to this heading"){.headerlink}

This section provides documentation for the modules and classes in the
project.

::: {#submodules .section}
## Submodules[¶](#submodules "Permalink to this heading"){.headerlink}
:::

::: {#module-tdcrpy.TDCRPy .section}
[]{#tdcrpy-tdcrpy-module}

## tdcrpy.TDCRPy module[¶](#module-tdcrpy.TDCRPy "Permalink to this heading"){.headerlink}

Created on Mon Jan 23 16:01:49 2023

A Monte-Carlo code to calculate detection efficiency in TDCR
measurements

\@author: Romain Coulon, Jialin Hu Bureau International des Poids et
Mesures

[[tdcrpy.TDCRPy.]{.pre}]{.sig-prename .descclassname}[[TDCRPy]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[L]{.pre}]{.n}*, *[[TD]{.pre}]{.n}*, *[[TAB]{.pre}]{.n}*, *[[TBC]{.pre}]{.n}*, *[[TAC]{.pre}]{.n}*, *[[Rad]{.pre}]{.n}*, *[[pmf_1]{.pre}]{.n}*, *[[N]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[V]{.pre}]{.n}*, *[[mode]{.pre}]{.n}*, *[[mode2]{.pre}]{.n}*, *[[Display]{.pre}]{.n}[[=]{.pre}]{.o}[[False]{.pre}]{.default_value}*, *[[barp]{.pre}]{.n}[[=]{.pre}]{.o}[[False]{.pre}]{.default_value}*[)]{.sig-paren}[¶](#tdcrpy.TDCRPy.TDCRPy "Permalink to this definition"){.headerlink}

:   This is the main function of the TDCRPy package running the
    Monte-Carlo Triple-to-Double Coincidence Ratio model. The
    computation is made for a given solution containing a radionuclide
    (or a mixture of radionuclides), a given volume of scintillator V
    and a given Birks constant kB.

    It can operates in two modes:

    > <div>
    >
    > --\> In mode="eff", it calculates the efficiency of the TDCR
    > system as a function of a value (triplet) of free parameter(s) L,
    > the measurement data is not used;
    >
    > ---\> In mode="res", it calculates the residual of the TDCR model
    > parametrized by a value (or triplet) of free parameter(s) L and
    > the measurement data TD, TAB, TBC, TAC.
    >
    > </div>

    also, two configuration can be set:

    > <div>
    >
    > --\> mode2="sym", where symmetry is considered between the 3
    > photomultiplier tubes - here L is a scalar and only the global
    > TDCR value TD is used as measurement data.
    >
    > ---\> mode2="asym", where an asymmetry between the 3
    > photomultiplier tubes is possible - here L is a triplet and only
    > the specific TDCR values TAB, TBC, TAC are used as measurement
    > data.
    >
    > </div>

    The parmeter N sets the number of Monte-Carlo trails used for the
    estimation. Each MC trial corresponds to a simulated radiactive
    decay. TDCRPY() used a set of fonctions from the
    tdcrpy.TDCR_model_lib module.

    Advanced settings can be configured in the config.toml file.

    > <div>
    >
    > --\> By default Y = True so that the analytical model is applied
    > for solution containing only pure beta emitting radionuclides. If
    > you would like to apply the MC calculation also for these
    > nuclides, set Y = False.
    >
    > ---\> If you would like to change the number of bins nE to
    > discretize the linear energy space for quenching calculation, you
    > can change nE_electron and nE_alpha parameters for respectively
    > electrons and alpha particles.
    >
    > </div>

    > <div>
    >
    > --\> By default the calculation is set for Ultima-Gold cocktail
    > mixed with a small amount of aqueous solution. You can adapt for a
    > specific scintillator by changing the density, the mean charge
    > number Z and the mean mass number A of the scintillator.
    >
    > </div>

    L[Float (if mode2="sym") or a tuple (if mode2="asym")]{.classifier}

    :   Free parameter in keV-1.

    TD[float]{.classifier}

    :   triple-to-double coincidence ratio. Not consider if
        mode2="asym". Not consider if mode2="asym".

    TAB[float]{.classifier}

    :   triple-to-double coincidence ratio (coincidences between channel
        A and B). Not consider if mode2="sym".

    TBC[float]{.classifier}

    :   triple-to-double coincidence ratio (coincidences between channel
        B and C). Not consider if mode2="sym".

    TAC[float]{.classifier}

    :   triple-to-double coincidence ratio (coincidences between channel
        A and C). Not consider if mode2="sym".

    Rad[string]{.classifier}

    :   List of radionuclides (eg. "H-3, Co-60").

    pmf_1[string]{.classifier}

    :   list of probability of each radionuclide (eg. "0.8, 0.2").

    N[integer]{.classifier}

    :   Number of Monte-Carlo trials. recommanded N\>10000 (see JCGM
        101). Not applied in the case of pure beta emitting
        radionuclides.

    kB[float]{.classifier}

    :   Birks constant in cm/keV.

    V[float]{.classifier}

    :   volume of the scintillator in ml. run only for 10 ml

    mode[string]{.classifier}

    :   "res" to return the residual, "eff" to return efficiencies.

    mode2[string]{.classifier}

    :   "sym" for symetrical model, "asym" for symetrical model.

    Display[Boolean, optional]{.classifier}

    :   "True" to display details on the decay sampling. The default is
        False.

    barp[Boolean, optional]{.classifier}

    :   "True" to display the calculation progress. The default is True.

    res[float]{.classifier}

    :   Residuals of the model compared the measurement data for (a)
        given free parmeters L. (only in mode="res")

    mean_efficiency_S[float]{.classifier}

    :   Estimation of the efficiency of single counting events. (only in
        mode="eff")

    std_efficiency_S[float]{.classifier}

    :   Standard uncertainty from calculation associated with the
        estimation of the efficiency of single counting events. (only in
        mode="eff")

    mean_efficiency_D[float]{.classifier}

    :   Estimation of the efficiency of logic sum of double
        coincidences. (only in mode="eff")

    std_efficiency_D[float]{.classifier}

    :   Standard uncertainty from calculation associated with the
        estimation of the efficiency of logic sum of double
        coincidences. (only in mode="eff")

    mean_efficiency_T[float]{.classifier}

    :   Estimation of the efficiency of triple coincidences. (only in
        mode="eff")

    std_efficiency_T[float]{.classifier}

    :   Standard uncertainty from calculation associated with the
        estimation of the efficiency of triple coincidences. (only in
        mode="eff")
:::

::: {#module-tdcrpy.TDCR_model_lib .section}
[]{#tdcrpy-tdcr-model-lib-module}

## tdcrpy.TDCR_model_lib module[¶](#module-tdcrpy.TDCR_model_lib "Permalink to this heading"){.headerlink}

Created on Mon Jan 23 16:04:46 2023

Library of function of the TDCRpy code

\@author: Romain Coulon, Jialin Hu Bureau International des Poids et
Mesures

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[E_quench_a]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[nE]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.E_quench_a "Permalink to this definition"){.headerlink}

:   This function calculate the quenched energy alpha particles
    according to the Birks model of scintillation quenching

    e[float]{.classifier}

    :   energy of the alpha particle in keV.

    kB[float]{.classifier}

    :   Birks constant in cm/keV.

    float

    :   Quenched energy in keV.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[E_quench_e]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[nE]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.E_quench_e "Permalink to this definition"){.headerlink}

:   This function calculate the quenched energy of electrons according
    to the Birks model of scintillation quenching

    e[float]{.classifier}

    :   energy of the electron in eV.

    kB[float]{.classifier}

    :   Birks constant in cm/MeV.

    float

    :   Quenched energy in eV.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[TicTocGenerator]{.pre}]{.sig-name .descname}[(]{.sig-paren}[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.TicTocGenerator "Permalink to this definition"){.headerlink}

:   Generator that returns time differences

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[clear_terminal]{.pre}]{.sig-name .descname}[(]{.sig-paren}[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.clear_terminal "Permalink to this definition"){.headerlink}

:   Function to clear the terminal screen

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[display_header]{.pre}]{.sig-name .descname}[(]{.sig-paren}[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.display_header "Permalink to this definition"){.headerlink}

:   Function to display the header.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[energie_dep_beta]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e_inci]{.pre}]{.n}*, *[[\*]{.pre}]{.o}*, *[[matrice10_1]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[1.0,]{.pre} [2.0,]{.pre} [3.0,]{.pre} [\...,]{.pre} [198.0,]{.pre} [199.0,]{.pre} [200.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [5.1e-05,]{.pre} [4.7e-05,]{.pre} [5.7e-05\],]{.pre} [\[0.0,]{.pre} [2e-06,]{.pre} [2e-06,]{.pre} [\...,]{.pre} [7.7e-05,]{.pre} [7.4e-05,]{.pre} [7e-05\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice10_2]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[200.0,]{.pre} [202.0,]{.pre} [204.0,]{.pre} [\...,]{.pre} [1996.0,]{.pre} [1998.0,]{.pre} [2000.0\],]{.pre} [\[5.7e-05,]{.pre} [4.9e-05,]{.pre} [5.6e-05,]{.pre} [\...,]{.pre} [0.001795,]{.pre} [0.001826,]{.pre} [0.00181\],]{.pre} [\[0.000494,]{.pre} [0.000489,]{.pre} [0.000492,]{.pre} [\...,]{.pre} [0.000912,]{.pre} [0.000912,]{.pre} [0.000913\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.48222,]{.pre} [6.8e-05,]{.pre} [7e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [5e-06,]{.pre} [0.481935,]{.pre} [6.6e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3e-06,]{.pre} [6e-06,]{.pre} [0.481512\]\])]{.pre}]{.default_value}*, *[[matrice10_3]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[2000.0,]{.pre} [2010.0,]{.pre} [2020.0,]{.pre} [\...,]{.pre} [9980.0,]{.pre} [9990.0,]{.pre} [10000.0\],]{.pre} [\[0.00181,]{.pre} [0.001776,]{.pre} [0.001775,]{.pre} [\...,]{.pre} [0.012191,]{.pre} [0.012222,]{.pre} [0.012252\],]{.pre} [\[0.004217,]{.pre} [0.004213,]{.pre} [0.004184,]{.pre} [\...,]{.pre} [0.004783,]{.pre} [0.004785,]{.pre} [0.004751\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.000646,]{.pre} [1.4e-05,]{.pre} [1.5e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [1e-06,]{.pre} [0.000645,]{.pre} [1.5e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.00065\]\])]{.pre}]{.default_value}*, *[[ed]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0002008,]{.pre} [0.001999,]{.pre} [0.009991\],]{.pre} [\[0.0004016,]{.pre} [0.003998,]{.pre} [0.019982\],]{.pre} [\...,]{.pre} [\[0.2006,]{.pre} [1.997,]{.pre} [9.981\],]{.pre} [\[0.2008,]{.pre} [1.999,]{.pre} [9.991\],]{.pre} [\[0.201,]{.pre} [2.001,]{.pre} [10.001\]\])]{.pre}]{.default_value}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.energie_dep_beta "Permalink to this definition"){.headerlink}

:   This function samples the energy deposited by an electron in the
    scintillator using response calculated by the Monte-Carlo code
    MCNP6.

    e_inci[float]{.classifier}

    :   energy of the electron in keV.

    matrice10_1[list\[list\], optional]{.classifier}

    :   response matrix for electrons in the range \[1-200\] keV and for
        a scintillator volume of 10 ml.

    matrice10_2[list\[list\], optional]{.classifier}

    :   response matrix for electrons in the range \[200-2000\] keV and
        for a scintillator volume of 10 ml.

    matrice10_3[list\[list\], optional]{.classifier}

    :   response matrix for electrons in the range \[2000-10000\] keV
        and for a scintillator volume of 10 ml.

    ed[list\[list\], optional]{.classifier}

    :   matrix of input energies. column 0: \[1-200\] keV; column 1:
        \[200-2000\] keV; column 2: \[2000-10000\] keV

    result[float]{.classifier}

    :   deposited energy in keV.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[energie_dep_beta2]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e_inci]{.pre}]{.n}*, *[[\*]{.pre}]{.o}*, *[[matrice10_1]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[1.0,]{.pre} [2.0,]{.pre} [3.0,]{.pre} [\...,]{.pre} [198.0,]{.pre} [199.0,]{.pre} [200.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [5.1e-05,]{.pre} [4.7e-05,]{.pre} [5.7e-05\],]{.pre} [\[0.0,]{.pre} [2e-06,]{.pre} [2e-06,]{.pre} [\...,]{.pre} [7.7e-05,]{.pre} [7.4e-05,]{.pre} [7e-05\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice10_2]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[200.0,]{.pre} [202.0,]{.pre} [204.0,]{.pre} [\...,]{.pre} [1996.0,]{.pre} [1998.0,]{.pre} [2000.0\],]{.pre} [\[5.7e-05,]{.pre} [4.9e-05,]{.pre} [5.6e-05,]{.pre} [\...,]{.pre} [0.001795,]{.pre} [0.001826,]{.pre} [0.00181\],]{.pre} [\[0.000494,]{.pre} [0.000489,]{.pre} [0.000492,]{.pre} [\...,]{.pre} [0.000912,]{.pre} [0.000912,]{.pre} [0.000913\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.48222,]{.pre} [6.8e-05,]{.pre} [7e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [5e-06,]{.pre} [0.481935,]{.pre} [6.6e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3e-06,]{.pre} [6e-06,]{.pre} [0.481512\]\])]{.pre}]{.default_value}*, *[[matrice10_3]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[2000.0,]{.pre} [2010.0,]{.pre} [2020.0,]{.pre} [\...,]{.pre} [9980.0,]{.pre} [9990.0,]{.pre} [10000.0\],]{.pre} [\[0.00181,]{.pre} [0.001776,]{.pre} [0.001775,]{.pre} [\...,]{.pre} [0.012191,]{.pre} [0.012222,]{.pre} [0.012252\],]{.pre} [\[0.004217,]{.pre} [0.004213,]{.pre} [0.004184,]{.pre} [\...,]{.pre} [0.004783,]{.pre} [0.004785,]{.pre} [0.004751\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.000646,]{.pre} [1.4e-05,]{.pre} [1.5e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [1e-06,]{.pre} [0.000645,]{.pre} [1.5e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.00065\]\])]{.pre}]{.default_value}*, *[[ed]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0002008,]{.pre} [0.001999,]{.pre} [0.009991\],]{.pre} [\[0.0004016,]{.pre} [0.003998,]{.pre} [0.019982\],]{.pre} [\...,]{.pre} [\[0.2006,]{.pre} [1.997,]{.pre} [9.981\],]{.pre} [\[0.2008,]{.pre} [1.999,]{.pre} [9.991\],]{.pre} [\[0.201,]{.pre} [2.001,]{.pre} [10.001\]\])]{.pre}]{.default_value}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.energie_dep_beta2 "Permalink to this definition"){.headerlink}

:   This function samples the energy deposited by an electron in the
    scintillator using response calculated by the Monte-Carlo code
    MCNP6.

    e_inci[float]{.classifier}

    :   energy of the electron in keV.

    matrice10_1[list\[list\], optional]{.classifier}

    :   response matrix for electrons in the range \[1-200\] keV and for
        a scintillator volume of 10 ml.

    matrice10_2[list\[list\], optional]{.classifier}

    :   response matrix for electrons in the range \[200-2000\] keV and
        for a scintillator volume of 10 ml.

    matrice10_3[list\[list\], optional]{.classifier}

    :   response matrix for electrons in the range \[2000-10000\] keV
        and for a scintillator volume of 10 ml.

    ed[list\[list\], optional]{.classifier}

    :   matrix of input energies. column 0: \[1-200\] keV; column 1:
        \[200-2000\] keV; column 2: \[2000-10000\] keV

    result[float]{.classifier}

    :   deposited energy in keV.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[energie_dep_gamma]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e_inci]{.pre}]{.n}*, *[[v]{.pre}]{.n}*, *[[matrice10_1]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[1.0,]{.pre} [2.0,]{.pre} [3.0,]{.pre} [\...,]{.pre} [198.0,]{.pre} [199.0,]{.pre} [200.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [1e-06,]{.pre} [2e-06,]{.pre} [2e-06\],]{.pre} [\[0.0,]{.pre} [0.001795,]{.pre} [0.005885,]{.pre} [\...,]{.pre} [0.873706,]{.pre} [0.87388,]{.pre} [0.87406\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice10_2]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[200.0,]{.pre} [202.0,]{.pre} [204.0,]{.pre} [\...,]{.pre} [1996.0,]{.pre} [1998.0,]{.pre} [2000.0\],]{.pre} [\[2e-06,]{.pre} [1e-06,]{.pre} [0.0,]{.pre} [\...,]{.pre} [4.5e-05,]{.pre} [4.2e-05,]{.pre} [4.8e-05\],]{.pre} [\[0.877508,]{.pre} [0.877865,]{.pre} [0.878157,]{.pre} [\...,]{.pre} [0.951362,]{.pre} [0.951381,]{.pre} [0.951434\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3e-06,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [4e-06,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [4e-06\]\])]{.pre}]{.default_value}*, *[[matrice10_3]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[2000.0,]{.pre} [2010.0,]{.pre} [2020.0,]{.pre} [\...,]{.pre} [9980.0,]{.pre} [9990.0,]{.pre} [10000.0\],]{.pre} [\[6e-05,]{.pre} [3e-05,]{.pre} [8e-05,]{.pre} [\...,]{.pre} [0.00017,]{.pre} [0.00019,]{.pre} [0.00015\],]{.pre} [\[0.95212,]{.pre} [0.95221,]{.pre} [0.95236,]{.pre} [\...,]{.pre} [0.97913,]{.pre} [0.97918,]{.pre} [0.97918\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice16_1]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[1.0,]{.pre} [2.0,]{.pre} [3.0,]{.pre} [\...,]{.pre} [198.0,]{.pre} [199.0,]{.pre} [200.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [2e-06,]{.pre} [1e-06,]{.pre} [3e-06\],]{.pre} [\[0.0,]{.pre} [0.001533,]{.pre} [0.005069,]{.pre} [\...,]{.pre} [0.855719,]{.pre} [0.855893,]{.pre} [0.856118\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice16_2]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[200.0,]{.pre} [202.0,]{.pre} [204.0,]{.pre} [\...,]{.pre} [1996.0,]{.pre} [1998.0,]{.pre} [2000.0\],]{.pre} [\[3e-06,]{.pre} [1e-06,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3.8e-05,]{.pre} [4.8e-05,]{.pre} [4.1e-05\],]{.pre} [\[0.859989,]{.pre} [0.860337,]{.pre} [0.860676,]{.pre} [\...,]{.pre} [0.943906,]{.pre} [0.943951,]{.pre} [0.944006\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3e-06,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [4e-06,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [4e-06\]\])]{.pre}]{.default_value}*, *[[matrice16_3]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[2000.0,]{.pre} [2010.0,]{.pre} [2020.0,]{.pre} [\...,]{.pre} [9980.0,]{.pre} [9990.0,]{.pre} [10000.0\],]{.pre} [\[4.1e-05,]{.pre} [3.1e-05,]{.pre} [4.4e-05,]{.pre} [\...,]{.pre} [0.000194,]{.pre} [0.00019,]{.pre} [0.000205\],]{.pre} [\[0.944572,]{.pre} [0.944712,]{.pre} [0.944892,]{.pre} [\...,]{.pre} [0.975521,]{.pre} [0.975553,]{.pre} [0.975526\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[ed]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0002008,]{.pre} [0.001999,]{.pre} [0.009991\],]{.pre} [\[0.0004016,]{.pre} [0.003998,]{.pre} [0.019982\],]{.pre} [\...,]{.pre} [\[0.2006,]{.pre} [1.997,]{.pre} [9.981\],]{.pre} [\[0.2008,]{.pre} [1.999,]{.pre} [9.991\],]{.pre} [\[0.201,]{.pre} [2.001,]{.pre} [10.001\]\])]{.pre}]{.default_value}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.energie_dep_gamma "Permalink to this definition"){.headerlink}

:   This function samples the energy deposited by a x or gamma rays in
    the scintillator using response calculated by the Monte-Carlo code
    MCNP6.

    e_inci[float]{.classifier}

    :   energy of the photon in keV.

    v[float]{.classifier}

    :   volume of the scintillator in ml.

    matrice10_1[list\[list\], optional]{.classifier}

    :   response matrix for photons in the range \[1-200\] keV and for a
        scintillator volume of 10 ml.

    matrice10_2[list\[list\], optional]{.classifier}

    :   response matrix for photons in the range \[200-2000\] keV and
        for a scintillator volume of 10 ml.

    matrice10_3[list\[list\], optional]{.classifier}

    :   response matrix for photons in the range \[2000-10000\] keV and
        for a scintillator volume of 10 ml.

    ed[list\[list\], optional]{.classifier}

    :   matrix of input energies. column 0: \[1-200\] keV; column 1:
        \[200-2000\] keV; column 2: \[2000-10000\] keV

    result[float]{.classifier}

    :   deposited energy in keV.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[energie_dep_gamma2]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e_inci]{.pre}]{.n}*, *[[v]{.pre}]{.n}*, *[[matrice10_1]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[1.0,]{.pre} [2.0,]{.pre} [3.0,]{.pre} [\...,]{.pre} [198.0,]{.pre} [199.0,]{.pre} [200.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [1e-06,]{.pre} [2e-06,]{.pre} [2e-06\],]{.pre} [\[0.0,]{.pre} [0.001795,]{.pre} [0.005885,]{.pre} [\...,]{.pre} [0.873706,]{.pre} [0.87388,]{.pre} [0.87406\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice10_2]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[200.0,]{.pre} [202.0,]{.pre} [204.0,]{.pre} [\...,]{.pre} [1996.0,]{.pre} [1998.0,]{.pre} [2000.0\],]{.pre} [\[2e-06,]{.pre} [1e-06,]{.pre} [0.0,]{.pre} [\...,]{.pre} [4.5e-05,]{.pre} [4.2e-05,]{.pre} [4.8e-05\],]{.pre} [\[0.877508,]{.pre} [0.877865,]{.pre} [0.878157,]{.pre} [\...,]{.pre} [0.951362,]{.pre} [0.951381,]{.pre} [0.951434\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3e-06,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [4e-06,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [4e-06\]\])]{.pre}]{.default_value}*, *[[matrice10_3]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[2000.0,]{.pre} [2010.0,]{.pre} [2020.0,]{.pre} [\...,]{.pre} [9980.0,]{.pre} [9990.0,]{.pre} [10000.0\],]{.pre} [\[6e-05,]{.pre} [3e-05,]{.pre} [8e-05,]{.pre} [\...,]{.pre} [0.00017,]{.pre} [0.00019,]{.pre} [0.00015\],]{.pre} [\[0.95212,]{.pre} [0.95221,]{.pre} [0.95236,]{.pre} [\...,]{.pre} [0.97913,]{.pre} [0.97918,]{.pre} [0.97918\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice16_1]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[1.0,]{.pre} [2.0,]{.pre} [3.0,]{.pre} [\...,]{.pre} [198.0,]{.pre} [199.0,]{.pre} [200.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [2e-06,]{.pre} [1e-06,]{.pre} [3e-06\],]{.pre} [\[0.0,]{.pre} [0.001533,]{.pre} [0.005069,]{.pre} [\...,]{.pre} [0.855719,]{.pre} [0.855893,]{.pre} [0.856118\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice16_2]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[200.0,]{.pre} [202.0,]{.pre} [204.0,]{.pre} [\...,]{.pre} [1996.0,]{.pre} [1998.0,]{.pre} [2000.0\],]{.pre} [\[3e-06,]{.pre} [1e-06,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3.8e-05,]{.pre} [4.8e-05,]{.pre} [4.1e-05\],]{.pre} [\[0.859989,]{.pre} [0.860337,]{.pre} [0.860676,]{.pre} [\...,]{.pre} [0.943906,]{.pre} [0.943951,]{.pre} [0.944006\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3e-06,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [4e-06,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [4e-06\]\])]{.pre}]{.default_value}*, *[[matrice16_3]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[2000.0,]{.pre} [2010.0,]{.pre} [2020.0,]{.pre} [\...,]{.pre} [9980.0,]{.pre} [9990.0,]{.pre} [10000.0\],]{.pre} [\[4.1e-05,]{.pre} [3.1e-05,]{.pre} [4.4e-05,]{.pre} [\...,]{.pre} [0.000194,]{.pre} [0.00019,]{.pre} [0.000205\],]{.pre} [\[0.944572,]{.pre} [0.944712,]{.pre} [0.944892,]{.pre} [\...,]{.pre} [0.975521,]{.pre} [0.975553,]{.pre} [0.975526\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[ed]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0002008,]{.pre} [0.001999,]{.pre} [0.009991\],]{.pre} [\[0.0004016,]{.pre} [0.003998,]{.pre} [0.019982\],]{.pre} [\...,]{.pre} [\[0.2006,]{.pre} [1.997,]{.pre} [9.981\],]{.pre} [\[0.2008,]{.pre} [1.999,]{.pre} [9.991\],]{.pre} [\[0.201,]{.pre} [2.001,]{.pre} [10.001\]\])]{.pre}]{.default_value}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.energie_dep_gamma2 "Permalink to this definition"){.headerlink}

:   This function samples the energy deposited by a x or gamma rays in
    the scintillator using response calculated by the Monte-Carlo code
    MCNP6.

    e_inci[float]{.classifier}

    :   energy of the photon in keV.

    v[float]{.classifier}

    :   volume of the scintillator in ml.

    matrice10_1[list\[list\], optional]{.classifier}

    :   response matrix for photons in the range \[1-200\] keV and for a
        scintillator volume of 10 ml.

    matrice10_2[list\[list\], optional]{.classifier}

    :   response matrix for photons in the range \[200-2000\] keV and
        for a scintillator volume of 10 ml.

    matrice10_3[list\[list\], optional]{.classifier}

    :   response matrix for photons in the range \[2000-10000\] keV and
        for a scintillator volume of 10 ml.

    matrice16_1[list\[list\], optional]{.classifier}

    :   response matrix for photons in the range \[1-200\] keV and for a
        scintillator volume of 16 ml.

    matrice16_2[list\[list\], optional]{.classifier}

    :   response matrix for photons in the range \[200-2000\] keV and
        for a scintillator volume of 16 ml.

    matrice16_3[list\[list\], optional]{.classifier}

    :   response matrix for photons in the range \[2000-10000\] keV and
        for a scintillator volume of 16 ml.

    ed[list\[list\], optional]{.classifier}

    :   matrix of input energies. column 0: \[1-200\] keV; column 1:
        \[200-2000\] keV; column 2: \[2000-10000\] keV

    result[float]{.classifier}

    :   deposited energy in keV.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[modelAnalytical]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[L]{.pre}]{.n}*, *[[TD]{.pre}]{.n}*, *[[TAB]{.pre}]{.n}*, *[[TBC]{.pre}]{.n}*, *[[TAC]{.pre}]{.n}*, *[[rad]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[V]{.pre}]{.n}*, *[[mode]{.pre}]{.n}*, *[[mode2]{.pre}]{.n}*, *[[ne]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.modelAnalytical "Permalink to this definition"){.headerlink}

:   TDCR analytical model that is used for pure beta emitting
    radionuclides

    L[float or tuple]{.classifier}

    :   free parameter(s).

    TD[float]{.classifier}

    :   triple-to-double coincidence ratio that was measured (logic
        sum).

    TAB[float]{.classifier}

    :   triple-to-double coincidence ratio that was measured (channels A
        and B).

    TBC[flat]{.classifier}

    :   triple-to-double coincidence ratio that was measured (channels B
        and C).

    TAC[float]{.classifier}

    :   triple-to-double coincidence ratio that was measured (channels A
        and C).

    rad[string]{.classifier}

    :   radionuclide (eg. "Na-22").

    kB[float]{.classifier}

    :   Birks constant in cm/keV.

    V[float]{.classifier}

    :   volume of the scintillator in ml. run only for 10 ml

    mode[string]{.classifier}

    :   "res" to return the residual, "eff" to return efficiencies.

    mode2[string]{.classifier}

    :   "sym" for symetrical model, "asym" for symetrical model.

    nE[integer]{.classifier}

    :   Number of bins for the quenching function.

    res[float]{.classifier}

    :   Residuals of the model compared the measurement data for (a)
        given free parmeters L. (only in mode="res")

    mean_efficiency_S[float]{.classifier}

    :   Estimation of the efficiency of single counting events. (only in
        mode="eff")

    mean_efficiency_D[float]{.classifier}

    :   Estimation of the efficiency of logic sum of double
        coincidences. (only in mode="eff")

    mean_efficiency_T[float]{.classifier}

    :   Estimation of the efficiency of triple coincidences. (only in
        mode="eff")

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[normalise]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[p_x]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.normalise "Permalink to this definition"){.headerlink}

:   This function is used to ensure that the sum of probability is equal
    to 1.

    p_x[list]{.classifier}

    :   vector of probabilities.

    p[list]{.classifier}

    :   normalized probability vector.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[readBetaShape]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[rad]{.pre}]{.n}*, *[[mode]{.pre}]{.n}*, *[[level]{.pre}]{.n}*, *[[z=\<zipfile.ZipFile]{.pre} [filename=\'C:\\\\Users\\\\romain.coulon\\\\AppData\\\\Local\\\\Continuum\\\\anaconda3\\\\lib\\\\site-packages\\\\tdcrpy\\\\decayData\\\\All-nuclides_BetaShape.zip\']{.pre} [mode=\'r\'\>]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.readBetaShape "Permalink to this definition"){.headerlink}

:   This funcion reads the beta spectra calculated by the code BetaShape
    and published in the DDEP web page.

    refs:

    > <div>
    >
    > <https://doi.org/10.1103/PhysRevC.92.059902>
    >
    > <http://www.lnhb.fr/ddep_wg/>
    >
    > </div>

    rad[string]{.classifier}

    :   identifier of the radionuclide. e.g. 'Na-22'

    mode[string]{.classifier}

    :   identifier of the decay mode. 'beta-' or 'beta+'

    level[int or string]{.classifier}

    :   level of the daughter after decay. 0,1,2,3 .... or 'tot' in case
        of pure beta emitting radionuclides

    e[list]{.classifier}

    :   the energy vector in keV.

    dNdx[list]{.classifier}

    :   the probability density in keV-1.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[readEShape]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[rad]{.pre}]{.n}*, *[[\*]{.pre}]{.n}*, *[[z=\<zipfile.ZipFile]{.pre} [filename=\'C:\\\\Users\\\\romain.coulon\\\\AppData\\\\Local\\\\Continuum\\\\anaconda3\\\\lib\\\\site-packages\\\\tdcrpy\\\\decayData\\\\All-nuclides_Ensdf.zip\']{.pre} [mode=\'r\'\>]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.readEShape "Permalink to this definition"){.headerlink}

:   This function reads the ENSDF zip files and format the data to be
    processed by TDCRPy.

    rad[string]{.classifier}

    :   name of the radionuclide such as 'Ag-108'.

    z[ZipFile object]{.classifier}

    :   zip ENSDF file.

    daug_name[list]{.classifier}

    :   daughter nucleus of the decay

    Energy[list]{.classifier}

    :   comprise all transition energies of the daughter nucleus.

    Prob[list]{.classifier}

    :   comprise all transtion probabilities of the daughter nucleus.

    Type[list]{.classifier}

    :   comprise all type of transition of the daughter nucleus.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[readPenNuc2]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[rad]{.pre}]{.n}*, *[[z1=\<zipfile.ZipFile]{.pre} [filename=\'C:\\\\Users\\\\romain.coulon\\\\AppData\\\\Local\\\\Continuum\\\\anaconda3\\\\lib\\\\site-packages\\\\tdcrpy\\\\decayData\\\\All-nuclides_PenNuc.zip\']{.pre} [mode=\'r\'\>]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.readPenNuc2 "Permalink to this definition"){.headerlink}

:   This function is used to read PenNuc files to format the decay data
    in lists readable by TDCRPy.

    rad[string]{.classifier}

    :   name of the radionculide (for example: "Am-241").

    daughter[list]{.classifier}

    :   list of the daughter nucleus -- indice 0.

    prob_daug[list]{.classifier}

    :   list of probabilities to produce daugter nuclei -- indice 1.

    energy_Q[list]{.classifier}

    :   list of Q value for each transition to a given daughter nucleus
        -- indice 2.

    desin_type_tot[list\[list\]]{.classifier}

    :   list of type of decay branch / emitted particules -- indice 3.
        It contains a sub-list for all possible branches of a given
        daughter nucleus and a sub-sub list related to possible decay
        mode of each branch.

    desin_energy_tot[list\[list\]]{.classifier}

    :   list of the energies of decay transition or the emitted
        particles -- indice 4. It contains a sub-list for all possible
        branches of a given daughter nucleus and a sub-sub list related
        to possible decay mode of each branch.

    desin_prob_tot[list\[list\]]{.classifier}

    :   list of the prabability of decay transition or the emitted
        particles -- indice 5. It contains a sub-list for all possible
        branches of a given daughter nucleus and a sub-sub list related
        to possible decay mode of each branch.

    desin_level_tot[list\[list\]]{.classifier}

    :   list of energy level that the daughter nucleus can have just
        after the decay of the mother nucleus -- indice 6. It contains a
        sub-list for all possible branches of a given daughter nucleus
        and a sub-sub list related to possible decay mode of each
        branch.

    prob_branch_tot[list]{.classifier}

    :   list of branch probabilities -- indice 7. It contains a sub-list
        for all possible branches of a given daughter nucleus.

    tran_type_tot[list\[list\]]{.classifier}

    :   list of all possible transitions -- indice 8. It contains a
        sub-list for all possible branches of a given daughter nucleus
        and a sub-sub list related to possible decay mode of each
        branch.

    tran_energy_tot[list\[list\]]{.classifier}

    :   list of energy associated with transitions -- indice 9. It
        contains a sub-list for all possible branches of a given
        daughter nucleus and a sub-sub list related to possible decay
        mode of each branch.

    tran_prob_tot[list\[list\]]{.classifier}

    :   list of probability associated with transitions -- indice 10. It
        contains a sub-list for all possible branches of a given
        daughter nucleus and a sub-sub list related to possible decay
        mode of each branch.

    tran_level_tot[list\[list\]]{.classifier}

    :   list of corresponding branch levels -- indice 11. It contains a
        sub-list for all possible branches of a given daughter nucleus
        and a sub-sub list related to the level before the transition.

    tran_level_end_tot[list\[list\]]{.classifier}

    :   list of level following given transitions -- indice 12. It
        contains a sub-list for all possible branches of a given
        daughter nucleus and a sub-sub list related to the level after
        the transition.

    level_energy_tot[list\[list\]]{.classifier}

    :   list of energy levels -- indice 13. It contains a sub-list for
        all possible branches of a given daughter nucleus and a sub-sub
        list related to possible decay mode of each branch.

    prob_tran_tot[list\[list\]]{.classifier}

    :   list of sum of transition of each branches -- indice 14. It
        contains a sub-list for all possible branches of a given
        daughter nucleus and a sub-sub list related to possible decay
        mode of each branch.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[read_matrice]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[path]{.pre}]{.n}*, *[[niveau]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.read_matrice "Permalink to this definition"){.headerlink}

:   This function read the response matrix calculated by MCNP6
    simulation.

    path[string]{.classifier}

    :   path to the response matrix.

    niveau[integer or string]{.classifier}

    :   energy range of the response matrix. 0: \[1-200\] keV; 1:
        \[200-2000\] keV; 2: \[2000-10000\] keV. "e" for the input
        energy matrix.

    matrice[list\[list\]]{.classifier}

    :   formatted response matrix.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[relaxation_atom]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[daugther]{.pre}]{.n}*, *[[rad]{.pre}]{.n}*, *[[lacune]{.pre}]{.n}[[=]{.pre}]{.o}[[\'defaut\']{.pre}]{.default_value}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.relaxation_atom "Permalink to this definition"){.headerlink}

:   This function simulates the atomic rearangement following a missing
    electron an inner shell of the daughter atom.

    daugther[string]{.classifier}

    :   The daughter nucleus (for example NB95,PD110 etc.)

    rad[string]{.classifier}

    :   The mother nucleus (for exemple Am-241, C-11 etc.)

    lacune[string]{.classifier}

    :   The shell where the electron is missing (for example
        'Atom_K','Atom_L' etc.)

    Type : type of transition Auger L or K, or X Ray. Energy :
    corresponding energy in keV.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[sampling]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[p_x]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.sampling "Permalink to this definition"){.headerlink}

:   This function aims to sample in a pdf or a pmf

    p_x[float vector]{.classifier}

    :   Probability Density (or mass) Function (PDF or PMF) of the
        random variable x.

    i[integer]{.classifier}

    :   index in x pointing the sampled value of the random variable X.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[stoppingpower]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e]{.pre}]{.n}*, *[[rho]{.pre}]{.n}[[=]{.pre}]{.o}[[0.96]{.pre}]{.default_value}*, *[[Z]{.pre}]{.n}[[=]{.pre}]{.o}[[5.2]{.pre}]{.default_value}*, *[[A]{.pre}]{.n}[[=]{.pre}]{.o}[[11.04]{.pre}]{.default_value}*, *[[emin]{.pre}]{.n}[[=]{.pre}]{.o}[[0]{.pre}]{.default_value}*, *[[file]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[0.0,]{.pre} [0.7,]{.pre} [1.4,]{.pre} [\...,]{.pre} [12.491336,]{.pre} [12.490668,]{.pre} [8.55292038e-72\])]{.pre}]{.default_value}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.stoppingpower "Permalink to this definition"){.headerlink}

:   The stopping power of electrons between 20 keV and 1000 keV is a
    mixture of a radiative loss model \[1\], and a collision model \[2\]
    that has been validated agaisnt the NIST model ESTAR \[3\]
    recommanded by the ICRU Report 37 \[4\]. At low energy - between 10
    eV and 20 keV - the model from Tan and Xia \[5\] is implemented.

    Refs:

    > <div>
    >
    > \[1\] <https://doi.org/10.1016/0020-708x(82)90244-7>
    >
    > \[2\]
    > <https://www.ijstr.org/final-print/jan2017/Calculations-Of-Stopping-Power-And-Range-Of-Electrons-Interaction-With-Different-Material-And-Human-Body-Parts.pdf>
    >
    > \[3\] <https://dx.doi.org/10.18434/T4NC7P>
    >
    > \[4\] ICRU Report 37, Stopping Powers for Electrons and Positrons
    >
    > \[5\] <https://doi.org/10.1016/j.apradiso.2011.08.012>
    >
    > </div>

    e[float]{.classifier}

    :   Energy of the electron in eV.

    rho[float, optional]{.classifier}

    :   density of the source in g.cm-3. The default is 0.96.

    Z[float, optional]{.classifier}

    :   mean charge number of the source. The default is 5.2.

    A[float, optional]{.classifier}

    :   mean mass number of the source. The default is 11.04.

    emin[float, optional]{.classifier}

    :   the minimal energy to consider. The default is 0.

    file[list, optional]{.classifier}

    :   tabulated data form the Tan and Xia model. The default is
        data_TanXia_f.

    dEdx[float]{.classifier}

    :   Calculated stopping power in MeV.cm-1.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[stoppingpowerA]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e]{.pre}]{.n}*, *[[rho]{.pre}]{.n}[[=]{.pre}]{.o}[[0.96]{.pre}]{.default_value}*, *[[energy_alpha]{.pre}]{.n}[[=]{.pre}]{.o}[[\[1.0,]{.pre} [1.5,]{.pre} [2.0,]{.pre} [2.5,]{.pre} [3.0,]{.pre} [4.0,]{.pre} [5.0,]{.pre} [6.0,]{.pre} [7.0,]{.pre} [8.0,]{.pre} [9.0,]{.pre} [10.0,]{.pre} [12.5,]{.pre} [15.0,]{.pre} [17.5,]{.pre} [20.0,]{.pre} [22.5,]{.pre} [25.0,]{.pre} [27.5,]{.pre} [30.0,]{.pre} [35.0,]{.pre} [40.0,]{.pre} [45.0,]{.pre} [50.0,]{.pre} [55.0,]{.pre} [60.0,]{.pre} [65.0,]{.pre} [70.0,]{.pre} [75.0,]{.pre} [80.0,]{.pre} [85.0,]{.pre} [90.0,]{.pre} [95.0,]{.pre} [100.0,]{.pre} [125.0,]{.pre} [150.0,]{.pre} [175.0,]{.pre} [200.0,]{.pre} [225.0,]{.pre} [250.0,]{.pre} [275.0,]{.pre} [300.0,]{.pre} [350.0,]{.pre} [400.0,]{.pre} [450.0,]{.pre} [500.0,]{.pre} [550.0,]{.pre} [600.0,]{.pre} [650.0,]{.pre} [700.0,]{.pre} [750.0,]{.pre} [800.0,]{.pre} [850.0,]{.pre} [900.0,]{.pre} [950.0,]{.pre} [1000.0,]{.pre} [1250.0,]{.pre} [1500.0,]{.pre} [1750.0,]{.pre} [2000.0,]{.pre} [2250.0,]{.pre} [2500.0,]{.pre} [2750.0,]{.pre} [3000.0,]{.pre} [3500.0,]{.pre} [4000.0,]{.pre} [4500.0,]{.pre} [5000.0,]{.pre} [5500.0,]{.pre} [6000.0,]{.pre} [6500.0,]{.pre} [7000.0,]{.pre} [7500.0,]{.pre} [8000.0\]]{.pre}]{.default_value}*, *[[dEdx_alpha]{.pre}]{.n}[[=]{.pre}]{.o}[[\[426600.0,]{.pre} [442500.0,]{.pre} [456000.0,]{.pre} [468700.0,]{.pre} [480800.0,]{.pre} [504300.0,]{.pre} [526800.0,]{.pre} [548500.0,]{.pre} [569400.0,]{.pre} [589600.0,]{.pre} [609100.0,]{.pre} [628100.0,]{.pre} [673000.0,]{.pre} [715000.0,]{.pre} [754500.0,]{.pre} [791900.0,]{.pre} [827400.0,]{.pre} [861200.0,]{.pre} [893600.0,]{.pre} [924700.0,]{.pre} [983500.0,]{.pre} [1038000.0,]{.pre} [1090000.0,]{.pre} [1139000.0,]{.pre} [1185000.0,]{.pre} [1229000.0,]{.pre} [1271000.0,]{.pre} [1311000.0,]{.pre} [1350000.0,]{.pre} [1387000.0,]{.pre} [1423000.0,]{.pre} [1457000.0,]{.pre} [1490000.0,]{.pre} [1522000.0,]{.pre} [1668000.0,]{.pre} [1792000.0,]{.pre} [1900000.0,]{.pre} [1993000.0,]{.pre} [2074000.0,]{.pre} [2145000.0,]{.pre} [2206000.0,]{.pre} [2258000.0,]{.pre} [2340000.0,]{.pre} [2398000.0,]{.pre} [2436000.0,]{.pre} [2457000.0,]{.pre} [2465000.0,]{.pre} [2463000.0,]{.pre} [2453000.0,]{.pre} [2436000.0,]{.pre} [2415000.0,]{.pre} [2389000.0,]{.pre} [2361000.0,]{.pre} [2331000.0,]{.pre} [2299000.0,]{.pre} [2267000.0,]{.pre} [2111000.0,]{.pre} [1962000.0,]{.pre} [1819000.0,]{.pre} [1684000.0,]{.pre} [1561000.0,]{.pre} [1457000.0,]{.pre} [1369000.0,]{.pre} [1292000.0,]{.pre} [1164000.0,]{.pre} [1061000.0,]{.pre} [977400.0,]{.pre} [907300.0,]{.pre} [847700.0,]{.pre} [796300.0,]{.pre} [751400.0,]{.pre} [711800.0,]{.pre} [676600.0,]{.pre} [645000.0\]]{.pre}]{.default_value}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.stoppingpowerA "Permalink to this definition"){.headerlink}

:   Estimation of the stopping power of alpha particles using tabulated
    values form the ASTAR code

    ref:

    > <div>
    >
    > <https://dx.doi.org/10.18434/T4NC7P>
    >
    > </div>

    e[float]{.classifier}

    :   energy of the alpha particle in keV.

    rho[float, optional]{.classifier}

    :   density of the source in g.cm-3. The default is 0.96.

    energy_alpha[list, optional]{.classifier}

    :   the list of energy (in keV) for which the stopping power was
        calculated with ASTAR. The default is energy_alph.

    dEdx_alpha[list, optional]{.classifier}

    :   the list of stopping powers (in keV.cm2/g) associated with the
        energy vector. The default is dEdx_alph.

    float

    :   Interpolated ASTAR estimation of the stopping power.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[tic]{.pre}]{.sig-name .descname}[(]{.sig-paren}[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.tic "Permalink to this definition"){.headerlink}

:   Records a time in TicToc, marks the beginning of a time interval

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[toc]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[tempBool]{.pre}]{.n}[[=]{.pre}]{.o}[[True]{.pre}]{.default_value}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.toc "Permalink to this definition"){.headerlink}

:   Prints the time difference yielded by generator instance TicToc

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[transf_name]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[rad]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.transf_name "Permalink to this definition"){.headerlink}

:   This function format the name of the nuclide to match with the
    PenNuc format.

    rad[string]{.classifier}

    :   name of the radionculdie such as '108AG'.

    RAD[string]{.classifier}

    :   name of the radionuclide such as 'AG108' that match with PenNuc
        format.

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[writeEffcurves]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[x]{.pre}]{.n}*, *[[y]{.pre}]{.n}*, *[[uy]{.pre}]{.n}*, *[[rad]{.pre}]{.n}*, *[[p]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[SDT]{.pre}]{.n}*[)]{.sig-paren}[¶](#tdcrpy.TDCR_model_lib.writeEffcurves "Permalink to this definition"){.headerlink}

:   This function writes efficiency curves

    x[TYPE]{.classifier}

    :   DESCRIPTION.

    y[TYPE]{.classifier}

    :   DESCRIPTION.

    uy[TYPE]{.classifier}

    :   DESCRIPTION.

    rad[TYPE]{.classifier}

    :   DESCRIPTION.

    p[TYPE]{.classifier}

    :   DESCRIPTION.

    kB[TYPE]{.classifier}

    :   DESCRIPTION.

    SDT[TYPE]{.classifier}

    :   DESCRIPTION.

    None.
:::

::: {#module-tdcrpy.TDCRoptimize .section}
[]{#tdcrpy-tdcroptimize-module}

## tdcrpy.TDCRoptimize module[¶](#module-tdcrpy.TDCRoptimize "Permalink to this heading"){.headerlink}

Created on Wed Jul 5 10:04:53 2023

\@author: romain.coulon, jialin.hu

[[tdcrpy.TDCRoptimize.]{.pre}]{.sig-prename .descclassname}[[eff]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[TD]{.pre}]{.n}*, *[[TAB]{.pre}]{.n}*, *[[TBC]{.pre}]{.n}*, *[[TAC]{.pre}]{.n}*, *[[Rad]{.pre}]{.n}*, *[[pmf_1]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[V]{.pre}]{.n}*, *[[mode2]{.pre}]{.n}*, *[[N]{.pre}]{.n}[[=]{.pre}]{.o}[[1000]{.pre}]{.default_value}*, *[[L]{.pre}]{.n}[[=]{.pre}]{.o}[[1]{.pre}]{.default_value}*[)]{.sig-paren}[¶](#tdcrpy.TDCRoptimize.eff "Permalink to this definition"){.headerlink}

:   Caclulation of the efficiency of a TDCR system based on the model
    TDCRPy.

    TD[float]{.classifier}

    :   triple-to-double coincidence ratio. Not consider if
        mode2="asym". Not consider if mode2="asym".

    TAB[float]{.classifier}

    :   triple-to-double coincidence ratio (coincidences between channel
        A and B). Not consider if mode2="sym".

    TBC[float]{.classifier}

    :   triple-to-double coincidence ratio (coincidences between channel
        B and C). Not consider if mode2="sym".

    TAC[float]{.classifier}

    :   triple-to-double coincidence ratio (coincidences between channel
        A and C). Not consider if mode2="sym".

    Rad[string]{.classifier}

    :   List of radionuclides.

    pmf_1[string]{.classifier}

    :   list of probability of each radionuclide..

    kB[float]{.classifier}

    :   Birks constant.

    V[float]{.classifier}

    :   volume of the scintillator in ml. run only for 10 ml

    mode2[string]{.classifier}

    :   "sym" for symetrical model, "asym" for symetrical model.

    N[interger, optional]{.classifier}

    :   number of Monte-Carlo trials. The default is 1000.

    L[float, optional]{.classifier}

    :   free parameter(s) as initial guess. The default is 1.

    L0[float]{.classifier}

    :   global free parameter.

    L[tuple]{.classifier}

    :   free parameters (relevant for the asymetric model).

    eff_S[float]{.classifier}

    :   counting efficiency of single events.

    u_eff_S[float]{.classifier}

    :   standard uncertainty of eff_S.

    eff_D[float]{.classifier}

    :   counting efficiency of double coincidences.

    u_eff_D[float]{.classifier}

    :   standard uncertainty of eff_D.

    eff_T[float]{.classifier}

    :   counting efficiency of triple coincidences.

    u_eff_T[float]{.classifier}

    :   standard uncertainty of eff_T.
:::

::: {#module-tdcrpy .section}
[]{#module-contents}

## Module contents[¶](#module-tdcrpy "Permalink to this heading"){.headerlink}
:::
:::
:::
:::
:::

::: {.sphinxsidebar role="navigation" aria-label="main navigation"}
::: {.sphinxsidebarwrapper}
# [TDCRPy](../index.html) {#tdcrpy .logo}

### Navigation

::: {.relations}
### Related Topics

-   [Documentation overview](../index.html)
:::

::: {#searchbox style="display: none" role="search"}
### Quick search {#searchlabel}

::: {.searchformwrapper}
:::
:::
:::
:::

::: {.clearer}
:::
:::

::: {.footer}
©2023, Romain Coulon, Jialin Hu. \| Powered by [Sphinx
5.0.2](http://sphinx-doc.org/) & [Alabaster
0.7.12](https://github.com/bitprophet/alabaster) \| [Page
source](../_sources/source/tdcrpy.rst.txt)
:::
