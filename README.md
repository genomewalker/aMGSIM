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
aMGSIM communities --n-comm 3 examples/data/genome_list.txt example
```

where [genome_table.tsv](examples/data/genome_list.txt) requires the following columns:

- **Taxon**: taxon name
- **Fasta**: genome fasta file path


The subcommand `communities` will generate three synthetic communities. From the [MGSIM]() documentation, the description of the files generated:

- **_example_\_abund.txt**: taxon relative abundances for each community this is the relative number of genome copies for each taxon
- **_example_\__wAbund.txt**: taxon relative abundances weighted by genome size this is the fraction of the DNA pool for each taxon
- **_example_\__beta-div.txt**: beta diversity among communities see the beta-div parameter for selecting beta-diversity measures


Once we have the composition of the synthetic communities we can create the ancient and modern partitions in each sample. We will use the subcommand `ancient-genomes`:

```
$ aMGSIM ancient-genomes -h
ancient-genomes: Estimate coverage, depth and other properties for each genome in each synthetic community

Usage:
  ancient-genomes [options] <config>
  ancient-genomes -h | --help
  ancient-genomes --version

Options:
  <config>       Config parameters
  -d --debug     Debug mode (no subprocesses; verbose output)
  -h --help      Show this screen.
  --version      Show version.

Description:
  Simulating ancient reads for each taxon in each synthetic community

  config
  ------
  * tab-delimited
  * must contain 3 columns
    * "Community" = community ID (ie., sample ID)
    * "Taxon" = taxon name
    * "Perc_rel_abund" = percent relative abundance of the taxon

  Output
  ------
  * A JSON file with the read properties for the selected ancient genomes
  * A TSV file with the read abundances for the ancient/modern genomes in each
    synthetic community
```

The subcommand takes as input a YAML [config](examples/config/ag-config.yml) file with the different parameters that will be use to create the set of genomes in each synthetic community. Some of the important parameters are `genome_table` (the same used as in the subcommand `communities`); `abund_table` the abundance table created with the subcommand `communities` or we can define our own values creating a table like:
```
Community   Taxon                                Perc_rel_abund   Rank
1           Escherichia_coli_K-12_MG1655         59.052261390     1
1           Methanosarcina_barkeri_MS            32.277831367     2
1           Clostridium_perfringens_ATCC_13124   8.669907243      3
2           Escherichia_coli_K-12_MG1655         84.638859282     1
2           Methanosarcina_barkeri_MS            13.605601555     2
2           Clostridium_perfringens_ATCC_13124   1.755539163      3
3           Escherichia_coli_K-12_MG1655         67.034789353     1
3           Methanosarcina_barkeri_MS            29.298568225     2
3           Clostridium_perfringens_ATCC_13124   3.666642422      3
```

And the `genome_comp`, where we can define for each sample which genome will be ancient and the coverage. [Here](examples/data/genome-comp.tsv) you can find an example of the file. It is a TSV file with four columns: 

```
Taxon                                Community   Coverage   onlyAncient
Escherichia_coli_K-12_MG1655         1           1.0        True
Methanosarcina_barkeri_MS            2           0.5        False
Clostridium_perfringens_ATCC_13124   3           20.0       True
```

The column `onlyAncient` can take one of the following values:

- **True**: The genome in this sample will only contain ancient reads with the coverage specified, 1.0 in the example
- **False**: The genome in this sample will contain a mixture of modern and ancient reads. The ancient coverage will be the one specified in the file, 10.0 in the example. In the case the coverage asked will produce

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
