#!/usr/bin/env python
from setuptools import setup, find_packages
import os
import glob
import numpy

    
# dependencies
install_reqs = [
    'docopt>=0.6.2',
    'numpy>=1.10',
    'pandas>=0.18',
    'scipy>=0.17',
    'configobj>=5.0.6',
    'biopython>=1.68',
    'pyfastx>=0.6.10',
    'schema>=0.7.0',
    'ruamel.yaml>=0.16',
    'jsonpickle>=1.4.1',
    'tqdm>=4.48',
    'pyfaidx>=0.5.9',
    'biolib>=0.1.6',
    'Cython>=0.29.21',
    'pyranges>=0.0.84',
    'ray>=0.8.7'
]

## install main application
desc = 'Anciengt metagenome simulation of multiple synthetic communities'
setup(
    name = 'aMGSIM',
    version = '0.1',
    description = desc,
    long_description = desc + '\n See README for more information.',
    author = 'Antonio Fernandez-Guerra',
    author_email = 'antonio@metagenomics.eu',
    entry_points={
        'console_scripts': [
            'aMGSIM = aMGSIM.__main__:main'
        ]
    },
    install_requires = install_reqs,
    include_dirs = [numpy.get_include()],
    license = "MIT license",
    packages = find_packages(),
    package_dir={'aMGSIM':
                 'aMGSIM'},
    url = 'https://github.com/genomewalker/aMGSIM'
)




