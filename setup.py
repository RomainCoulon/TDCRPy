from setuptools import setup, find_packages
import codecs
import os

VERSION = "0.0.6"
DESCRIPTION = "TDCR model"

setup(
    name = "TDCRPy",
    version = VERSION,
    author = "RomainCoulon (Romain Coulon)",
    author_email = "<romain.coulon@bipm.org>",
    description = DESCRIPTION,
    packages = find_packages(),
    install_requires = ["numpy","scipy","sys","time","urllib","zipfile","re"],
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
    ]
)