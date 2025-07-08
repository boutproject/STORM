import os
from setuptools import setup, find_packages


# Utility function to read the README file.
# Used for the long_description.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="xstorm",
    version="0.1",
    url="https://gitlab.com/CCFE_SOL_Transport/xstorm/tree/master",
    author="Thomas Nicholas",
    author_email="thomas.nicholas@york.ac.uk",
    description="Analyse data from STORM simulations using xarray",
    license="Apache",
    long_description=read("README.md"),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache License",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    install_requires=["xarray", "xbout", "xrft"],
    packages=find_packages(),
    include_package_data=True,
)
