::: {.wy-grid-for-nav}
::: {.wy-side-scroll}
::: {.wy-side-nav-search}
[TDCRPy](../index.html){.icon .icon-home}

::: {role="search"}
:::
:::

::: {.wy-menu .wy-menu-vertical spy="affix" role="navigation" aria-label="Navigation menu"}
::: {.local-toc}
-   [tdcrpy package](#){.reference .internal}
    -   [Submodules](#submodules){.reference .internal}
    -   [tdcrpy.Activity_TDCR
        module](#tdcrpy-activity-tdcr-module){.reference .internal}
    -   [tdcrpy.EfficiencyProfils
        module](#tdcrpy-efficiencyprofils-module){.reference .internal}
    -   [tdcrpy.TDCRPy module](#module-tdcrpy.TDCRPy){.reference
        .internal}
        -   [`TDCRPy()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCRPy.TDCRPy){.reference .internal}
    -   [tdcrpy.TDCR_model_lib
        module](#module-tdcrpy.TDCR_model_lib){.reference .internal}
        -   [`E_quench_a()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.E_quench_a){.reference
            .internal}
        -   [`E_quench_e()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.E_quench_e){.reference
            .internal}
        -   [`TicTocGenerator()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.TicTocGenerator){.reference
            .internal}
        -   [`clear_terminal()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.clear_terminal){.reference
            .internal}
        -   [`display_header()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.display_header){.reference
            .internal}
        -   [`energie_dep_beta()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.energie_dep_beta){.reference
            .internal}
        -   [`energie_dep_gamma()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.energie_dep_gamma){.reference
            .internal}
        -   [`modelAnalytical()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.modelAnalytical){.reference
            .internal}
        -   [`normalise()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.normalise){.reference
            .internal}
        -   [`readBetaShape()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.readBetaShape){.reference
            .internal}
        -   [`readEShape()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.readEShape){.reference
            .internal}
        -   [`readPenNuc2()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.readPenNuc2){.reference
            .internal}
        -   [`read_matrice()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.read_matrice){.reference
            .internal}
        -   [`relaxation_atom()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.relaxation_atom){.reference
            .internal}
        -   [`sampling()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.sampling){.reference
            .internal}
        -   [`stoppingpower()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.stoppingpower){.reference
            .internal}
        -   [`stoppingpowerA()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.stoppingpowerA){.reference
            .internal}
        -   [`tic()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.tic){.reference
            .internal}
        -   [`toc()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.toc){.reference
            .internal}
        -   [`transf_name()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.transf_name){.reference
            .internal}
        -   [`writeEffcurves()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCR_model_lib.writeEffcurves){.reference
            .internal}
    -   [tdcrpy.TDCRoptimize
        module](#module-tdcrpy.TDCRoptimize){.reference .internal}
        -   [`eff()`{.docutils .literal
            .notranslate}](#tdcrpy.TDCRoptimize.eff){.reference
            .internal}
    -   [tdcrpy.decay module](#tdcrpy-decay-module){.reference
        .internal}
    -   [Module contents](#module-tdcrpy){.reference .internal}
:::
:::
:::

::: {.section .wy-nav-content-wrap toggle="wy-nav-shift"}
[TDCRPy](../index.html)

::: {.wy-nav-content}
::: {.rst-content}
::: {role="navigation" aria-label="Page navigation"}
-   [](../index.html){.icon .icon-home}
-   tdcrpy package
-   [View page source](../_sources/source/tdcrpy.rst.txt)

------------------------------------------------------------------------
:::

::: {.document role="main" itemscope="itemscope" itemtype="http://schema.org/Article"}
::: {itemprop="articleBody"}
::: {#tdcrpy-package .section}
# tdcrpy package[](#tdcrpy-package "Permalink to this heading"){.headerlink}

::: {#submodules .section}
## Submodules[](#submodules "Permalink to this heading"){.headerlink}
:::

::: {#tdcrpy-activity-tdcr-module .section}
## tdcrpy.Activity_TDCR module[](#tdcrpy-activity-tdcr-module "Permalink to this heading"){.headerlink}
:::

::: {#tdcrpy-efficiencyprofils-module .section}
## tdcrpy.EfficiencyProfils module[](#tdcrpy-efficiencyprofils-module "Permalink to this heading"){.headerlink}
:::

::: {#module-tdcrpy.TDCRPy .section}
[]{#tdcrpy-tdcrpy-module}

## tdcrpy.TDCRPy module[](#module-tdcrpy.TDCRPy "Permalink to this heading"){.headerlink}

Created on Mon Jan 23 16:01:49 2023

A Monte-Carlo code to calculate detection efficiency in TDCR
measurements

\@author: Romain Coulon, Jialin Hu Bureau International des Poids et
Mesures

[[tdcrpy.TDCRPy.]{.pre}]{.sig-prename .descclassname}[[TDCRPy]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[L]{.pre}]{.n}*, *[[TD]{.pre}]{.n}*, *[[TAB]{.pre}]{.n}*, *[[TBC]{.pre}]{.n}*, *[[TAC]{.pre}]{.n}*, *[[Rad]{.pre}]{.n}*, *[[pmf_1]{.pre}]{.n}*, *[[N]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[mode]{.pre}]{.n}*, *[[mode2]{.pre}]{.n}*, *[[Display]{.pre}]{.n}[[=]{.pre}]{.o}[[False]{.pre}]{.default_value}*, *[[barp]{.pre}]{.n}[[=]{.pre}]{.o}[[True]{.pre}]{.default_value}*[)]{.sig-paren}[](#tdcrpy.TDCRPy.TDCRPy "Permalink to this definition"){.headerlink}

:   This is a Monte-Carlo TDCR model

    ::: {#parameters .section}
    ### Parameters[](#parameters "Permalink to this heading"){.headerlink}

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

    :   Number of Monte-Carlo trials. recommanded N\>10000. Not applied
        in the case of pure beta emitting radionuclides.

    kB[float]{.classifier}

    :   Birks constant in cm/keV.

    mode[string]{.classifier}

    :   "res" to return the residual, "eff" to return efficiencies.

    mode2[string]{.classifier}

    :   "sym" for symetrical model, "asym" for symetrical model.

    Display[Boolean, optional]{.classifier}

    :   "True" to display details on the decay sampling. The default is
        False.

    barp[Boolean, optional]{.classifier}

    :   "True" to display the calculation progress. The default is True.
    :::

    ::: {#returns .section}
    ### Returns[](#returns "Permalink to this heading"){.headerlink}

    Tuple

    :   if mode=="res", the residual (float). if mode=="eff", the
        efficiencies (list)
    :::
:::

::: {#module-tdcrpy.TDCR_model_lib .section}
[]{#tdcrpy-tdcr-model-lib-module}

## tdcrpy.TDCR_model_lib module[](#module-tdcrpy.TDCR_model_lib "Permalink to this heading"){.headerlink}

Created on Mon Jan 23 16:04:46 2023

Library of function of the TDCRpy code

\@author: Romain Coulon, Jialin Hu Bureau International des Poids et
Mesures

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[E_quench_a]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[nE]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.E_quench_a "Permalink to this definition"){.headerlink}

:   This function calculate the quenched energy alpha particles
    according to the Birks model of scintillation quenching

    ::: {#id1 .section}
    ### Parameters[](#id1 "Permalink to this heading"){.headerlink}

    e[float]{.classifier}

    :   energy of the alpha particle in keV.

    kB[float]{.classifier}

    :   Birks constant in cm/keV.
    :::

    ::: {#id2 .section}
    ### Returns[](#id2 "Permalink to this heading"){.headerlink}

    float

    :   Quenched energy in keV.
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[E_quench_e]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[nE]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.E_quench_e "Permalink to this definition"){.headerlink}

:   This function calculate the quenched energy of electrons according
    to the Birks model of scintillation quenching

    ::: {#id3 .section}
    ### Parameters[](#id3 "Permalink to this heading"){.headerlink}

    e[float]{.classifier}

    :   energy of the electron in eV.

    kB[float]{.classifier}

    :   Birks constant in cm/MeV.
    :::

    ::: {#id4 .section}
    ### Returns[](#id4 "Permalink to this heading"){.headerlink}

    float

    :   Quenched energy in eV.
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[TicTocGenerator]{.pre}]{.sig-name .descname}[(]{.sig-paren}[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.TicTocGenerator "Permalink to this definition"){.headerlink}

:   Generator that returns time differences

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[clear_terminal]{.pre}]{.sig-name .descname}[(]{.sig-paren}[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.clear_terminal "Permalink to this definition"){.headerlink}

:   

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[display_header]{.pre}]{.sig-name .descname}[(]{.sig-paren}[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.display_header "Permalink to this definition"){.headerlink}

:   

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[energie_dep_beta]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e_inci]{.pre}]{.n}*, *[[\*]{.pre}]{.o}*, *[[matrice10_1]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[1.0,]{.pre} [2.0,]{.pre} [3.0,]{.pre} [\...,]{.pre} [198.0,]{.pre} [199.0,]{.pre} [200.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [5.1e-05,]{.pre} [4.7e-05,]{.pre} [5.7e-05\],]{.pre} [\[0.0,]{.pre} [2e-06,]{.pre} [2e-06,]{.pre} [\...,]{.pre} [7.7e-05,]{.pre} [7.4e-05,]{.pre} [7e-05\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice10_2]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[200.0,]{.pre} [202.0,]{.pre} [204.0,]{.pre} [\...,]{.pre} [1996.0,]{.pre} [1998.0,]{.pre} [2000.0\],]{.pre} [\[5.7e-05,]{.pre} [4.9e-05,]{.pre} [5.6e-05,]{.pre} [\...,]{.pre} [0.001795,]{.pre} [0.001826,]{.pre} [0.00181\],]{.pre} [\[0.000494,]{.pre} [0.000489,]{.pre} [0.000492,]{.pre} [\...,]{.pre} [0.000912,]{.pre} [0.000912,]{.pre} [0.000913\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.48222,]{.pre} [6.8e-05,]{.pre} [7e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [5e-06,]{.pre} [0.481935,]{.pre} [6.6e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3e-06,]{.pre} [6e-06,]{.pre} [0.481512\]\])]{.pre}]{.default_value}*, *[[matrice10_3]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[2000.0,]{.pre} [2010.0,]{.pre} [2020.0,]{.pre} [\...,]{.pre} [9980.0,]{.pre} [9990.0,]{.pre} [10000.0\],]{.pre} [\[0.00181,]{.pre} [0.001776,]{.pre} [0.001775,]{.pre} [\...,]{.pre} [0.012191,]{.pre} [0.012222,]{.pre} [0.012252\],]{.pre} [\[0.004217,]{.pre} [0.004213,]{.pre} [0.004184,]{.pre} [\...,]{.pre} [0.004783,]{.pre} [0.004785,]{.pre} [0.004751\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.000646,]{.pre} [1.4e-05,]{.pre} [1.5e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [1e-06,]{.pre} [0.000645,]{.pre} [1.5e-05\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.00065\]\])]{.pre}]{.default_value}*, *[[ed]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0002008,]{.pre} [0.001999,]{.pre} [0.009991\],]{.pre} [\[0.0004016,]{.pre} [0.003998,]{.pre} [0.019982\],]{.pre} [\...,]{.pre} [\[0.2006,]{.pre} [1.997,]{.pre} [9.981\],]{.pre} [\[0.2008,]{.pre} [1.999,]{.pre} [9.991\],]{.pre} [\[0.201,]{.pre} [2.001,]{.pre} [10.001\]\])]{.pre}]{.default_value}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.energie_dep_beta "Permalink to this definition"){.headerlink}

:   ::: {#id5 .section}
    ### Parameters[](#id5 "Permalink to this heading"){.headerlink}

    e_inci[float]{.classifier}

    :   l'énergie incidente de particule.

    -   

        : TYPE

        :   DESCRIPTION.

    matrice10_1[matrix]{.classifier}

    :   matrice d'électrons de 1-200keV de solution 10ml.

    matrice2[TYPE, optional]{.classifier}

    :   matrice d'électrons de 200-2000keV de solution 10ml.

    matrice3[TYPE, optional]{.classifier}

    :   matrice d'électrons de 2000-10000keV de solution 10ml.

    ed[TYPE, optional]{.classifier}

    :   matrice de bins d'énergie. colonne 0: 1-200keV; colonne 1:
        200-2000keV

    ::: {#id6 .section}
    #### Returns[](#id6 "Permalink to this heading"){.headerlink}

    result[float]{.classifier}

    :   l'énergie déposée.
    :::
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[energie_dep_gamma]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e_inci]{.pre}]{.n}*, *[[v]{.pre}]{.n}*, *[[matrice10_1]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[1.0,]{.pre} [2.0,]{.pre} [3.0,]{.pre} [\...,]{.pre} [198.0,]{.pre} [199.0,]{.pre} [200.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [1e-06,]{.pre} [2e-06,]{.pre} [2e-06\],]{.pre} [\[0.0,]{.pre} [0.001795,]{.pre} [0.005885,]{.pre} [\...,]{.pre} [0.873706,]{.pre} [0.87388,]{.pre} [0.87406\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice10_2]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[200.0,]{.pre} [202.0,]{.pre} [204.0,]{.pre} [\...,]{.pre} [1996.0,]{.pre} [1998.0,]{.pre} [2000.0\],]{.pre} [\[2e-06,]{.pre} [1e-06,]{.pre} [0.0,]{.pre} [\...,]{.pre} [4.5e-05,]{.pre} [4.2e-05,]{.pre} [4.8e-05\],]{.pre} [\[0.877508,]{.pre} [0.877865,]{.pre} [0.878157,]{.pre} [\...,]{.pre} [0.951362,]{.pre} [0.951381,]{.pre} [0.951434\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3e-06,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [4e-06,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [4e-06\]\])]{.pre}]{.default_value}*, *[[matrice10_3]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[2000.0,]{.pre} [2010.0,]{.pre} [2020.0,]{.pre} [\...,]{.pre} [9980.0,]{.pre} [9990.0,]{.pre} [10000.0\],]{.pre} [\[6e-05,]{.pre} [3e-05,]{.pre} [8e-05,]{.pre} [\...,]{.pre} [0.00017,]{.pre} [0.00019,]{.pre} [0.00015\],]{.pre} [\[0.95212,]{.pre} [0.95221,]{.pre} [0.95236,]{.pre} [\...,]{.pre} [0.97913,]{.pre} [0.97918,]{.pre} [0.97918\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice16_1]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[1.0,]{.pre} [2.0,]{.pre} [3.0,]{.pre} [\...,]{.pre} [198.0,]{.pre} [199.0,]{.pre} [200.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [2e-06,]{.pre} [1e-06,]{.pre} [3e-06\],]{.pre} [\[0.0,]{.pre} [0.001533,]{.pre} [0.005069,]{.pre} [\...,]{.pre} [0.855719,]{.pre} [0.855893,]{.pre} [0.856118\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[matrice16_2]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[200.0,]{.pre} [202.0,]{.pre} [204.0,]{.pre} [\...,]{.pre} [1996.0,]{.pre} [1998.0,]{.pre} [2000.0\],]{.pre} [\[3e-06,]{.pre} [1e-06,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3.8e-05,]{.pre} [4.8e-05,]{.pre} [4.1e-05\],]{.pre} [\[0.859989,]{.pre} [0.860337,]{.pre} [0.860676,]{.pre} [\...,]{.pre} [0.943906,]{.pre} [0.943951,]{.pre} [0.944006\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [3e-06,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [4e-06,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [4e-06\]\])]{.pre}]{.default_value}*, *[[matrice16_3]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[2000.0,]{.pre} [2010.0,]{.pre} [2020.0,]{.pre} [\...,]{.pre} [9980.0,]{.pre} [9990.0,]{.pre} [10000.0\],]{.pre} [\[4.1e-05,]{.pre} [3.1e-05,]{.pre} [4.4e-05,]{.pre} [\...,]{.pre} [0.000194,]{.pre} [0.00019,]{.pre} [0.000205\],]{.pre} [\[0.944572,]{.pre} [0.944712,]{.pre} [0.944892,]{.pre} [\...,]{.pre} [0.975521,]{.pre} [0.975553,]{.pre} [0.975526\],]{.pre} [\...,]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [\...,]{.pre} [0.0,]{.pre} [0.0,]{.pre} [0.0\]\])]{.pre}]{.default_value}*, *[[ed]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[\[0.0,]{.pre} [0.0,]{.pre} [0.0\],]{.pre} [\[0.0002008,]{.pre} [0.001999,]{.pre} [0.009991\],]{.pre} [\[0.0004016,]{.pre} [0.003998,]{.pre} [0.019982\],]{.pre} [\...,]{.pre} [\[0.2006,]{.pre} [1.997,]{.pre} [9.981\],]{.pre} [\[0.2008,]{.pre} [1.999,]{.pre} [9.991\],]{.pre} [\[0.201,]{.pre} [2.001,]{.pre} [10.001\]\])]{.pre}]{.default_value}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.energie_dep_gamma "Permalink to this definition"){.headerlink}

:   ::: {#id7 .section}
    ### Parameters[](#id7 "Permalink to this heading"){.headerlink}

    e_inci[float]{.classifier}

    :   l'énergie incidente de particule.

    -   

        : TYPE

        :   DESCRIPTION.

    matrice10_1[matrix]{.classifier}

    :   matrice de photon de 1-200keV de solution 10ml.

    matrice2[TYPE, optional]{.classifier}

    :   matrice de photon de 200-2000keV de solution 10ml.

    matrice3[TYPE, optional]{.classifier}

    :   matrice de photon de 2000-10000keV de solution 10ml.

    ed[TYPE, optional]{.classifier}

    :   matrice de bins d'énergie. colonne 0: 1-200keV; colonne 1:
        200-2000keV

    ::: {#id8 .section}
    #### Returns[](#id8 "Permalink to this heading"){.headerlink}

    result[float]{.classifier}

    :   l'énergie déposée.
    :::
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[modelAnalytical]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[L]{.pre}]{.n}*, *[[TD]{.pre}]{.n}*, *[[TAB]{.pre}]{.n}*, *[[TBC]{.pre}]{.n}*, *[[TAC]{.pre}]{.n}*, *[[rad]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[mode]{.pre}]{.n}*, *[[mode2]{.pre}]{.n}*, *[[ne]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.modelAnalytical "Permalink to this definition"){.headerlink}

:   TDCR analytical model that can be used for pure beta emitting
    radionuclides

    ::: {#id9 .section}
    ### Parameters[](#id9 "Permalink to this heading"){.headerlink}

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

    mode[string]{.classifier}

    :   "res" to return the residual, "eff" to return efficiencies.

    mode2[string]{.classifier}

    :   "sym" for symetrical model, "asym" for symetrical model.

    nE[integer]{.classifier}

    :   Number of bins for the quenching function.
    :::

    ::: {#id10 .section}
    ### Returns[](#id10 "Permalink to this heading"){.headerlink}

    Tuple

    :   if mode=="res", the residual (float). if mode=="eff", the
        efficiencies (list)
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[normalise]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[p_x]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.normalise "Permalink to this definition"){.headerlink}

:   This function is used to ensure that the sum of probability is equal
    to 1.

    ::: {#id11 .section}
    ### Parameters[](#id11 "Permalink to this heading"){.headerlink}

    p_x[list]{.classifier}

    :   vector of probabilities.
    :::

    ::: {#id12 .section}
    ### Returns[](#id12 "Permalink to this heading"){.headerlink}

    p[list]{.classifier}

    :   normalized probability vector.
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[readBetaShape]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[rad]{.pre}]{.n}*, *[[mode]{.pre}]{.n}*, *[[level]{.pre}]{.n}*, *[[z=\<zipfile.ZipFile]{.pre} [filename=\'C:\\\\Users\\\\romain.coulon\\\\Anaconda3\\\\lib\\\\site-packages\\\\tdcrpy\\\\decayData\\\\All-nuclides_BetaShape.zip\']{.pre} [mode=\'r\'\>]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.readBetaShape "Permalink to this definition"){.headerlink}

:   This funcion reads the beta spectra calculated by the code BetaShape
    and published in the DDEP web page. refs:

    > <div>
    >
    > <https://doi.org/10.1103/PhysRevC.92.059902>
    > <http://www.lnhb.fr/ddep_wg/>
    >
    > </div>

    ::: {#id13 .section}
    ### Parameters[](#id13 "Permalink to this heading"){.headerlink}

    rad[string]{.classifier}

    :   identifier of the radionuclide. e.g. 'Na-22'

    mode[string]{.classifier}

    :   identifier of the decay mode. 'beta-' or 'beta+'

    level[int]{.classifier}

    :   level of the daughter after decay. 0,1,2,3 ....
    :::

    ::: {#id14 .section}
    ### Returns[](#id14 "Permalink to this heading"){.headerlink}

    e[list]{.classifier}

    :   the energy vector in keV.

    dNdx[list]{.classifier}

    :   the probability density in keV-1.
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[readEShape]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[rad]{.pre}]{.n}*, *[[\*]{.pre}]{.n}*, *[[z=\<zipfile.ZipFile]{.pre} [filename=\'C:\\\\Users\\\\romain.coulon\\\\Anaconda3\\\\lib\\\\site-packages\\\\tdcrpy\\\\decayData\\\\All-nuclides_Ensdf.zip\']{.pre} [mode=\'r\'\>]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.readEShape "Permalink to this definition"){.headerlink}

:   ::: {#pour-lire-les-fichiers-dans-all-nuclides-ensdf-zip .section}
    ### pour lire les fichiers dans All-nuclides_Ensdf.zip[](#pour-lire-les-fichiers-dans-all-nuclides-ensdf-zip "Permalink to this heading"){.headerlink}
    :::

    ::: {#parametre .section}
    ### PARAMETRE[](#parametre "Permalink to this heading"){.headerlink}

    rad -- type: str par exemple 'Ag-108' z -- ENSDF files
    :::

    ::: {#return .section}
    ### RETURN[](#return "Permalink to this heading"){.headerlink}

    daug_name -- type: list -- les filles de désintégration Energy -----
    type: list -- chaque élément comprend toutes les énergies de
    transition de la fille de même indice Prob ------- type: list --
    chaque élément comprend toutes les proba de transition de la fille
    de même indice Type ------- type: list -- chaque élément comprend
    touts les types de transition de la fille de même indice
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[readPenNuc2]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[rad]{.pre}]{.n}*, *[[z1=\<zipfile.ZipFile]{.pre} [filename=\'C:\\\\Users\\\\romain.coulon\\\\Anaconda3\\\\lib\\\\site-packages\\\\tdcrpy\\\\decayData\\\\All-nuclides_PenNuc.zip\']{.pre} [mode=\'r\'\>]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.readPenNuc2 "Permalink to this definition"){.headerlink}

:   This function is used to read PenNuc files to format the decay data
    in lists.

    > <div>
    >
    > rad -- type: str (par exemple: "Am-241") -- radionucléide
    >
    > daughter -- indice 0 -- des noyaux fils -- len = nb de noyaux fils
    > prob_daug -- indice 1 -- des probabilités de noyaux fils -- len =
    > nb de noyaux fils energy_Q -- indice 2 -- des énergies de
    > désintégrations -- len = nb de noyaux fils
    >
    > desin_type_tot -- indice 3 -- des types de désintégrations/particules émis
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des désintégrations possibles de
    >     chaque branch -- len de sous-list de sous-list = nb de type de
    >     désintégrations de chaque branch
    >
    > desin_energy_tot -- indice 4 -- des énergies de désintégrations/énergies de patricules émis
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des énergies de désintégrations
    >     possibles de chaque branch -- len de sous-list de sous-list =
    >     nb de type de désintégrations de chaque branch
    >
    > desin_prob_tot -- indice 5 -- des probabilités de désintégrations
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des probabilités de désintégrations
    >     possibles de chaque branch -- len de sous-list de sous-list =
    >     nb de type de désintégrations de chaque branch
    >
    > desin_level_tot -- indice 6 -- des niveaux atteints après des désintégrations
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des niveaux après des
    >     désintégrations de chaque branch -- len de sous-list de
    >     sous-list = nb de type de désintégrations de chaque branch
    >
    > prob_branch_tot -- indice 7 -- probabilités de chaque branch
    >
    > :   len = nb de noyaux fils sous-list -- des probabilités de
    >     branchs de noyau fil -- len de sous-list = nb de branch de
    >     chaque fil
    >
    > tran_type_tot -- indice 8 -- transitions possibles
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des transitions possibles de chaque
    >     branch -- len de sous-list de sous-list = nb de type de
    >     transitions de chaque branch
    >
    > tran_energy_tot -- indice 9 -- énergies de transitions
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des énergies de transitions
    >     possibles de chaque branch -- len de sous-list de sous-list =
    >     nb de type de transitions de chaque branch
    >
    > tran_prob_tot -- indice 10 -- probabilités de transitions
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des probabilités de transitions
    >     possibles de chaque branch -- len de sous-list de sous-list =
    >     nb de type de transitions de chaque branch
    >
    > tran_level_tot -- indice 11 -- niveaux de branch correspondants
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des niveaux de chaque branch avant
    >     des transitions -- len de sous-list de sous-list = 1
    >
    > tran_level_end_tot -- indice 12 -- niveaux après des transitions
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des niveaux après des transitions de
    >     chaque branch -- len de sous-list de sous-list = nb de type de
    >     transitions de chaque branch
    >
    > level_energy_tot -- indice 13 -- énergies de niveaux
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des énergies de niveaux de chaque
    >     branch -- len de sous-list de sous-list = 1
    >
    > prob_tran_tot -- indice 14 -- la somme de transition de chaque branch
    >
    > :   len = nb de noyaux fils sous-list -- des branchs possibles de
    >     noyau fil -- len de sous-list = nb de branch de chaque fil
    >     sous-list de sous-list -- des énergies de niveaux de chaque
    >     branch -- len de sous-list de sous-list = 1
    >
    > </div>

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[read_matrice]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[path]{.pre}]{.n}*, *[[niveau]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.read_matrice "Permalink to this definition"){.headerlink}

:   

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[relaxation_atom]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[daugther]{.pre}]{.n}*, *[[rad]{.pre}]{.n}*, *[[lacune]{.pre}]{.n}[[=]{.pre}]{.o}[[\'defaut\']{.pre}]{.default_value}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.relaxation_atom "Permalink to this definition"){.headerlink}

:   ::: {#id15 .section}
    ### PARAMETRE[](#id15 "Permalink to this heading"){.headerlink}

    daugther -- type: str -- la fille tirée dans cette itération (par
    exemple NB95,PD110 etc.) rad ------- type: str -- le radionucléide
    étudié (par exemple Am-241, C-11 etc.) lacuen ---- type: str -- la
    lacune atomique (par exemple 'Atom_K','Atom_L' etc.)
    :::

    ::: {#id16 .section}
    ### RETURN[](#id16 "Permalink to this heading"){.headerlink}

    Type ---- type de transition: Auger L ou K ou Rayon X Energy --
    énergie correspondante
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[sampling]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[p_x]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.sampling "Permalink to this definition"){.headerlink}

:   This function aims to sample in a pdf or a pmf

    ::: {#id17 .section}
    ### Parameters[](#id17 "Permalink to this heading"){.headerlink}

    p_x[float vector]{.classifier}

    :   Probability Density (or mass) Function (PDF or PMF) of the
        random variable x.
    :::

    ::: {#id18 .section}
    ### Returns[](#id18 "Permalink to this heading"){.headerlink}

    i[integer]{.classifier}

    :   index in x pointing the sampled value of the random variable X.
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[stoppingpower]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e]{.pre}]{.n}*, *[[rho]{.pre}]{.n}[[=]{.pre}]{.o}[[0.96]{.pre}]{.default_value}*, *[[Z]{.pre}]{.n}[[=]{.pre}]{.o}[[5.2]{.pre}]{.default_value}*, *[[A]{.pre}]{.n}[[=]{.pre}]{.o}[[11.04]{.pre}]{.default_value}*, *[[emin]{.pre}]{.n}[[=]{.pre}]{.o}[[0]{.pre}]{.default_value}*, *[[file]{.pre}]{.n}[[=]{.pre}]{.o}[[array(\[0.0,]{.pre} [0.7,]{.pre} [1.4,]{.pre} [\...,]{.pre} [12.491336,]{.pre} [12.490668,]{.pre} [8.98619049e-43\])]{.pre}]{.default_value}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.stoppingpower "Permalink to this definition"){.headerlink}

:   The stopping power of electrons between 20 keV and 1000 keV is a
    mixture of a radiative loss model \[1\], and a collision model \[2\]
    that has been validated agaisnt the NIST model ESTAR \[3\]
    recommanded by the ICRU Report 37 \[4\]. At low energy - between 10
    eV and 20 keV - the model from Tan and Xia \[5\] is implemented.
    Refs:

    > <div>
    >
    > \[1\] <https://doi.org/10.1016/0020-708x(82)90244-7> \[2\]
    > <https://www.ijstr.org/final-print/jan2017/Calculations-Of-Stopping-Power-And-Range-Of-Electrons-Interaction-With-Different-Material-And-Human-Body-Parts.pdf>
    > \[3\] <https://dx.doi.org/10.18434/T4NC7P> \[4\] ICRU Report 37,
    > Stopping Powers for Electrons and Positrons \[5\]
    > <https://doi.org/10.1016/j.apradiso.2011.08.012>
    >
    > </div>

    ::: {#id19 .section}
    ### Parameters[](#id19 "Permalink to this heading"){.headerlink}

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
    :::

    ::: {#id20 .section}
    ### Returns[](#id20 "Permalink to this heading"){.headerlink}

    dEdx[float]{.classifier}

    :   Calculated stopping power in MeV.cm-1.
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[stoppingpowerA]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[e]{.pre}]{.n}*, *[[rho]{.pre}]{.n}[[=]{.pre}]{.o}[[0.96]{.pre}]{.default_value}*, *[[energy_alpha]{.pre}]{.n}[[=]{.pre}]{.o}[[\[1.0,]{.pre} [1.5,]{.pre} [2.0,]{.pre} [2.5,]{.pre} [3.0,]{.pre} [4.0,]{.pre} [5.0,]{.pre} [6.0,]{.pre} [7.0,]{.pre} [8.0,]{.pre} [9.0,]{.pre} [10.0,]{.pre} [12.5,]{.pre} [15.0,]{.pre} [17.5,]{.pre} [20.0,]{.pre} [22.5,]{.pre} [25.0,]{.pre} [27.5,]{.pre} [30.0,]{.pre} [35.0,]{.pre} [40.0,]{.pre} [45.0,]{.pre} [50.0,]{.pre} [55.0,]{.pre} [60.0,]{.pre} [65.0,]{.pre} [70.0,]{.pre} [75.0,]{.pre} [80.0,]{.pre} [85.0,]{.pre} [90.0,]{.pre} [95.0,]{.pre} [100.0,]{.pre} [125.0,]{.pre} [150.0,]{.pre} [175.0,]{.pre} [200.0,]{.pre} [225.0,]{.pre} [250.0,]{.pre} [275.0,]{.pre} [300.0,]{.pre} [350.0,]{.pre} [400.0,]{.pre} [450.0,]{.pre} [500.0,]{.pre} [550.0,]{.pre} [600.0,]{.pre} [650.0,]{.pre} [700.0,]{.pre} [750.0,]{.pre} [800.0,]{.pre} [850.0,]{.pre} [900.0,]{.pre} [950.0,]{.pre} [1000.0,]{.pre} [1250.0,]{.pre} [1500.0,]{.pre} [1750.0,]{.pre} [2000.0,]{.pre} [2250.0,]{.pre} [2500.0,]{.pre} [2750.0,]{.pre} [3000.0,]{.pre} [3500.0,]{.pre} [4000.0,]{.pre} [4500.0,]{.pre} [5000.0,]{.pre} [5500.0,]{.pre} [6000.0,]{.pre} [6500.0,]{.pre} [7000.0,]{.pre} [7500.0,]{.pre} [8000.0\]]{.pre}]{.default_value}*, *[[dEdx_alpha]{.pre}]{.n}[[=]{.pre}]{.o}[[\[426600.0,]{.pre} [442500.0,]{.pre} [456000.0,]{.pre} [468700.0,]{.pre} [480800.0,]{.pre} [504300.0,]{.pre} [526800.0,]{.pre} [548500.0,]{.pre} [569400.0,]{.pre} [589600.0,]{.pre} [609100.0,]{.pre} [628100.0,]{.pre} [673000.0,]{.pre} [715000.0,]{.pre} [754500.0,]{.pre} [791900.0,]{.pre} [827400.0,]{.pre} [861200.0,]{.pre} [893600.0,]{.pre} [924700.0,]{.pre} [983500.0,]{.pre} [1038000.0,]{.pre} [1090000.0,]{.pre} [1139000.0,]{.pre} [1185000.0,]{.pre} [1229000.0,]{.pre} [1271000.0,]{.pre} [1311000.0,]{.pre} [1350000.0,]{.pre} [1387000.0,]{.pre} [1423000.0,]{.pre} [1457000.0,]{.pre} [1490000.0,]{.pre} [1522000.0,]{.pre} [1668000.0,]{.pre} [1792000.0,]{.pre} [1900000.0,]{.pre} [1993000.0,]{.pre} [2074000.0,]{.pre} [2145000.0,]{.pre} [2206000.0,]{.pre} [2258000.0,]{.pre} [2340000.0,]{.pre} [2398000.0,]{.pre} [2436000.0,]{.pre} [2457000.0,]{.pre} [2465000.0,]{.pre} [2463000.0,]{.pre} [2453000.0,]{.pre} [2436000.0,]{.pre} [2415000.0,]{.pre} [2389000.0,]{.pre} [2361000.0,]{.pre} [2331000.0,]{.pre} [2299000.0,]{.pre} [2267000.0,]{.pre} [2111000.0,]{.pre} [1962000.0,]{.pre} [1819000.0,]{.pre} [1684000.0,]{.pre} [1561000.0,]{.pre} [1457000.0,]{.pre} [1369000.0,]{.pre} [1292000.0,]{.pre} [1164000.0,]{.pre} [1061000.0,]{.pre} [977400.0,]{.pre} [907300.0,]{.pre} [847700.0,]{.pre} [796300.0,]{.pre} [751400.0,]{.pre} [711800.0,]{.pre} [676600.0,]{.pre} [645000.0\]]{.pre}]{.default_value}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.stoppingpowerA "Permalink to this definition"){.headerlink}

:   Estimation of the stopping power of alpha particles using tabulated
    values form the ASTAR code ref: <https://dx.doi.org/10.18434/T4NC7P>

    ::: {#id21 .section}
    ### Parameters[](#id21 "Permalink to this heading"){.headerlink}

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
    :::

    ::: {#id22 .section}
    ### Returns[](#id22 "Permalink to this heading"){.headerlink}

    float

    :   Interpolated ASTAR estimation of the stopping power.
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[tic]{.pre}]{.sig-name .descname}[(]{.sig-paren}[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.tic "Permalink to this definition"){.headerlink}

:   Records a time in TicToc, marks the beginning of a time interval

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[toc]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[tempBool]{.pre}]{.n}[[=]{.pre}]{.o}[[True]{.pre}]{.default_value}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.toc "Permalink to this definition"){.headerlink}

:   Prints the time difference yielded by generator instance TicToc

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[transf_name]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[rad]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.transf_name "Permalink to this definition"){.headerlink}

:   ::: {#id23 .section}
    ### PARAMETRE[](#id23 "Permalink to this heading"){.headerlink}

    rad -- type: str par exemple '108AG'
    :::

    ::: {#id24 .section}
    ### RETURN[](#id24 "Permalink to this heading"){.headerlink}

    RAD -- type: str par exemple 'AG108' qui correspond à la structure
    de fille de PenNuc
    :::

```{=html}
<!-- -->
```

[[tdcrpy.TDCR_model_lib.]{.pre}]{.sig-prename .descclassname}[[writeEffcurves]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[x]{.pre}]{.n}*, *[[y]{.pre}]{.n}*, *[[uy]{.pre}]{.n}*, *[[rad]{.pre}]{.n}*, *[[p]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[SDT]{.pre}]{.n}*[)]{.sig-paren}[](#tdcrpy.TDCR_model_lib.writeEffcurves "Permalink to this definition"){.headerlink}

:   
:::

::: {#module-tdcrpy.TDCRoptimize .section}
[]{#tdcrpy-tdcroptimize-module}

## tdcrpy.TDCRoptimize module[](#module-tdcrpy.TDCRoptimize "Permalink to this heading"){.headerlink}

Created on Wed Jul 5 10:04:53 2023

\@author: romain.coulon, jialin.hu

[[tdcrpy.TDCRoptimize.]{.pre}]{.sig-prename .descclassname}[[eff]{.pre}]{.sig-name .descname}[(]{.sig-paren}*[[TD]{.pre}]{.n}*, *[[TAB]{.pre}]{.n}*, *[[TBC]{.pre}]{.n}*, *[[TAC]{.pre}]{.n}*, *[[Rad]{.pre}]{.n}*, *[[pmf_1]{.pre}]{.n}*, *[[kB]{.pre}]{.n}*, *[[mode2]{.pre}]{.n}*, *[[N]{.pre}]{.n}[[=]{.pre}]{.o}[[1000]{.pre}]{.default_value}*, *[[L]{.pre}]{.n}[[=]{.pre}]{.o}[[1]{.pre}]{.default_value}*[)]{.sig-paren}[](#tdcrpy.TDCRoptimize.eff "Permalink to this definition"){.headerlink}

:   Caclulation of the efficiency of a TDCR system based on the model
    TDCRPy

    ::: {#id25 .section}
    ### Parameters[](#id25 "Permalink to this heading"){.headerlink}

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

    mode2[string]{.classifier}

    :   "sym" for symetrical model, "asym" for symetrical model.

    N[interger, optional]{.classifier}

    :   number of Monte-Carlo trials. The default is 1000.

    L[float, optional]{.classifier}

    :   free parameter(s) as initial guess. The default is 1.
    :::

    ::: {#id26 .section}
    ### Returns[](#id26 "Permalink to this heading"){.headerlink}

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
:::

::: {#tdcrpy-decay-module .section}
## tdcrpy.decay module[](#tdcrpy-decay-module "Permalink to this heading"){.headerlink}
:::

::: {#module-tdcrpy .section}
[]{#module-contents}

## Module contents[](#module-tdcrpy "Permalink to this heading"){.headerlink}
:::
:::
:::
:::

------------------------------------------------------------------------

::: {role="contentinfo"}
© Copyright 2023, Romain Coulon, Jialin Hu.
:::

Built with [Sphinx](https://www.sphinx-doc.org/) using a
[theme](https://github.com/readthedocs/sphinx_rtd_theme) provided by
[Read the Docs](https://readthedocs.org).
:::
:::
:::
:::
