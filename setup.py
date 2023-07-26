from setuptools import setup, find_packages
import codecs
import os

VERSION = "0.0.50"
DESCRIPTION = "TDCR model"

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name = "TDCRPy",
    version = VERSION,
    author = "RomainCoulon (Romain Coulon)",
    author_email = "<romain.coulon@bipm.org>",
    description = DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RomainCoulon/TDCRPy",
    packages = find_packages(exclude=["EfficiencyProfils"]),
    install_requires = ["numpy","tqdm","setuptools","scipy","configparser","importlib.resources"],
    keywords = ["python","TDCR","Monte-Carlo","radionuclide","scintillation","counting"],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Natural Language :: French",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    include_package_data = True,
    package_data = {'': [
	    'decayData/All-nuclides_PenNuc.zip',
        'decayData/All-nuclides_BetaShape.zip',
		'decayData/All-nuclides_Ensdf.zip',
        'Quenching/alpha_toulene.txt',
		'Quenching/TandataUG.txt',
		'MCNP-MATRIX/matrice/fichier/*.txt',
		]},
)