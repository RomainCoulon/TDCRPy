from setuptools import setup, find_packages
import codecs
import os

VERSION = "2.15.2"

DESCRIPTION = "TDCR model"

md_files = ["README.md"]

long_description = ""
for md_file in md_files:
    with open(md_file, "r") as f:
        long_description += f.read() + "\n\n"

setup(
    name = "TDCRPy",
    version = VERSION,
    author = "RomainCoulon (Romain Coulon)",
    author_email = "<romain.coulon@bipm.org>",
    description = DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://pypi.org/project/TDCRPy/",
    py_modules=['TDCRPy'],
    project_urls={'Documentation': 'https://github.com/RomainCoulon/TDCRPy/',},
    # packages = find_packages(exclude=["tdcrpy.EfficiencyProfils","tdcrpy.decay","tdcrpy.Activity_TDCR"], include=["tdcrpy.TDCR_model_lib","tdcrpy.TDCRPy", "tdcrpy.test.test_tdcrpy"]),
    packages = find_packages(),
    install_requires = ["numpy","tqdm","setuptools","scipy","configparser","importlib.resources","matplotlib"],
    python_requires='>=3.11',
    keywords = ["Python","TDCR","Monte-Carlo","radionuclide","scintillation","counting"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
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
		'config.toml',
	    	'decayData/All-nuclides_PenNuc.zip',
        	'decayData/All-nuclides_BetaShape.zip',
		'decayData/All-nuclides_Ensdf.zip',
		'decayData/photo-ENDF.zip',
		'decayData/atom-ENDF-VII0.zip',
        	'Quenching/alpha_toulene.txt',
		'Quenching/TandataUG.txt',
		'Quenching/QuenchEnergy*.txt',
		'Quenching/inputVecteur*.txt',
		'MCNP-MATRIX/matrice/fichier/*.txt',
		'MCNP-MATRIX/Spectra_for_analytical_model/*.txt',
		'EfficiencyCurves/*/*.txt',
                'Micelle/*.csv',
		]},
)
