![MGSIM](https://github.com/genomewalker/aMGSIM/workflows/aMGSIM/badge.svg)

aMGSIM
=====

This is an extension of [MGSIM](https://github.com/nick-youngblut/MGSIM/) to create ancient metagenome read simulation of multiple synthetic communities. It integrates the methods in [Gargammel](https://github.com/grenaud/gargammel) and provides flexibility to create different experimental scenarios.

## Sections

- [REFERENCE](#reference)
- [INSTALLATION](#installation)
- [TUTORIALS](#tutorials)
- [SIMULATION WORKFLOW](#simulation_workflow)
- [CHANGE LOG](#changelog)
- [LICENSE](#license)


# REFERENCE
Please cite the original [MGSIM](https://github.com/nick-youngblut/MGSIM/)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3696891.svg)](https://doi.org/10.5281/zenodo.3696891)


# INSTALLATION

## Dependencies

aMGSIM depends on the binaries included in Gargammel. You can install it as:

`conda install -c bioconda gargammel`

## Install

#### via pip

`pip install MGSIM`

#### via `setup.py`

`python setpy.py install`

# HOW-TO

aMGSIM inetgrates all the programs from MGSIM, check the documentation of MGSIM for more details

See all subcommands:

`MGSIM --list`

## Download genomes

`MGSIM genome_download -h`

## Simulate communities

`MGSIM communities -h`

## Simulate reads for each genome in each community

### Simulating Illumina reads

`MGSIM reads -h`

### Simulating haplotagging reads

`MGSIM ht_reads -h`


# LICENSE

See [LICENSE](./LICENSE)


