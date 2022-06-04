#!/usr/bin/env python
from setuptools import setup, find_packages
import versioneer

# dependencies
install_reqs = [
    "docopt>=0.6.2",
    "numpy>=1.20.1",
    "pandas>=1.2.2",
    "scipy>=0.17",
    "configobj>=5.0.6",
    "biopython>=1.68",
    "pyfastx>=0.8.1",
    "schema>=0.7.0",
    "ruamel.yaml>=0.16",
    "tqdm>=4.48",
    "pyfaidx>=0.5.9",
    "biolib>=0.1.6",
    "sorted-nearest<=0.0.33",
    "pyranges>=0.0.112",
    "ray>=0.8.7",
    "MGSIM>=0.2.2",
    "pandarallel>=1.5.2",
    "taxopy>=0.9.2",
    "ujson>=5.1.0",
]

## install main application
desc = "Ancient metagenome simulation of multiple synthetic communities"
setup(
    setup_requires=[
        # Setuptools 18.0 properly handles Cython extensions.
        "setuptools>=18.0",
        "Cython>=0.29.21",
    ],
    name="aMGSIM",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="aMGSIM: simulate ancient metagenomes for multiple synthetic communities",
    long_description=desc + "\n See README for more information.",
    author="Antonio Fernandez-Guerra",
    author_email="antonio@metagenomics.eu",
    entry_points={"console_scripts": ["aMGSIM = aMGSIM.__main__:main"]},
    install_requires=install_reqs,
    license="MIT license",
    packages=find_packages(),
    package_dir={"aMGSIM": "aMGSIM"},
    url="https://github.com/genomewalker/aMGSIM",
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
)
