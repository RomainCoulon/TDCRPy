# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
import os

VERSION = "2.20.12"
DESCRIPTION = "TDCR model — Monte Carlo efficiency estimation for liquid scintillation counting"

# Read long description from README
with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="TDCRPy",
    version=VERSION,
    author="Romain Coulon",
    author_email="romain.coulon@bipm.org",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RomainCoulon/TDCRPy",   # ← repo, not PyPI
    project_urls={
        "Documentation": "https://github.com/RomainCoulon/TDCRPy/",
        "Bug Tracker":   "https://github.com/RomainCoulon/TDCRPy/issues",
    },
    packages=find_packages(),
    # py_modules removed — redundant with find_packages()
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "tqdm",
        "numba",
        "setuptools",   # only needed if you use pkg_resources at runtime
        # configparser and importlib.resources are stdlib in Python ≥3.11
    ],
    extras_require={
        "dev": ["pytest", "sphinx", "sphinx-rtd-theme"],
    },
    python_requires=">=3.11",
    keywords=["TDCR", "Monte-Carlo", "radionuclide", "liquid scintillation", "counting"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Natural Language :: French",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",   # ← simpler than listing all three
        "Topic :: Scientific/Engineering :: Physics",
    ],
    include_package_data=True,
    package_data={
        "": [
            "config.toml",
            "configDefault.toml",
            "decayData/All-nuclides_PenNuc.zip",
            "decayData/All-nuclides_BetaShape.zip",
            "decayData/All-nuclides_Ensdf.zip",
            "decayData/photo-ENDF.zip",
            "decayData/atom-ENDF-VII0.zip",
            "Quenching/alpha_toulene.txt",
            "Quenching/TandataUG.txt",
            "Quenching/QuenchEnergy*.txt",
            "Quenching/inputVecteur*.txt",
            "MCNP-MATRIX/matrice/fichier/*.txt",
            "MCNP-MATRIX/Spectra_for_analytical_model/*.txt",
            "EfficiencyCurves/*/*.txt",
            "Micelle/*.csv",
        ]
    },
)