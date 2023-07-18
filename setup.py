from setuptools import setup, find_packages
import codecs
import os

VERSION = "0.0.8"
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