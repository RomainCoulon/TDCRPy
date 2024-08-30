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
