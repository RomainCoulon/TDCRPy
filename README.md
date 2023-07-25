# TDCRPy

`TDCRPy` is a Python code to calculate detection efficiency of a liquide scintillation counter using 3-photomultiplier tubes.
The calculation is based on the photo-physical model called of the Triple-to-Double-Coincidence-Ratio method (TDCR) [[1]](#1) and a Monte-Carlo sampling allowing to adress complexe decay schemes and radionuclide mixtures.

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

## Documentation

The full documentation for this project can be found in [docs/_build/html/index.html](docs/_build/html/index.html).
