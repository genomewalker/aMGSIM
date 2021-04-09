# aMGSIM

[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/aMGSIM?include_prereleases&label=version)](https://github.com/genomewalker/aMGSIM/releases) [![aMGSIM](https://github.com/genomewalker/aMGSIM/workflows/aMGSIM_ci/badge.svg)](https://github.com/genomewalker/aMGSIM/actions) [![PyPI](https://img.shields.io/pypi/v/aMGSIM)](https://pypi.org/project/aMGSIM/) [![Conda](https://img.shields.io/conda/v/genomewalker/aMGSIM)](https://anaconda.org/genomewalker/aMGSIM)


This is an extension of [MGSIM](https://github.com/nick-youngblut/MGSIM/) to create simulated ancient metagenome reads in multiple microbial synthetic communities. It integrates the methods in [Gargammel](https://github.com/grenaud/gargammel) and provides flexibility to create different experimental scenarios. Furthermore, aMGSIM also has the ability to track damage over the codons of the predicted proteins from the microbial genomes. 



## INSTALLATION

We recommend to have [**conda**](https://docs.conda.io/en/latest/) installed to manage the virtual environments

### Using pip

First we create a conda virtual environment with:

```bash
wget https://raw.githubusercontent.com/genomewalker/aMGSIM/master/environment.yml
conda env create -f environment.yml
```

Then we proceed to install using pip:

```bash
pip install aMGSIM
```

### Using conda

```bash
conda install -c conda-forge -c bioconda -c genomewalker aMGSIM
```

### Install from source to use the development version

By cloning in a dedicated conda environment

```bash
git clone git@github.com:genomewalker/aMGSIM.git
cd aMGSIM
conda env create -f environment.yml
conda activate aMGSIM
pip install -e .
```

## Usage

aMGSIM uses as starting point the `communities` subcommand from [MGSIM](https://github.com/nick-youngblut/MGSIM/) and the integrates three new subcommands:

- **ancient-genomes**: Estimate coverage, depth and other properties for each genome in each synthetic community
- **ancient-reads**: Simulate ancient reads for each taxon in each synthetic community
- **protein-analysis**: Tracking damage to the codon positions of each simulated read. 

You can access to the list of commands by:

```bash
$ aMGSIM --list

Available Commands:
communities | ancient-genomes | ancient-reads | protein-analysis
```



First we need to generate a synthetic community using the subcommand `communities`
```
aMGSIM communities --n-comm 5 examples/data/genome_list.txt test
```

where [genome_table.tsv](examples/data/genome_list.txt) requires the following columns:
```
    * "Taxon" = taxon name
    * "Fasta" = genome fasta file path
```

Once we have the composition of the synthetic communities we can 

```
aMGSIM ancient-genomes  ag-config.yaml
aMGSIM ancient-reads examples/data/genome_list.txt ar-global-config.yaml
aMGSIM protein-analysis --cpus 3 --procs 16 example_abund_read-files.json
```

We need a file that specifies 

```

```
## LICENSE
See [LICENSE](./LICENSE)

## REFERENCE
Please cite the original [MGSIM](https://github.com/nick-youngblut/MGSIM/)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3696891.svg)](https://doi.org/10.5281/zenodo.3696891)
