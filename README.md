# aMGSIM

[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/aMGSIM?include_prereleases&label=version)](https://github.com/genomewalker/aMGSIM/releases) [![aMGSIM](https://github.com/genomewalker/aMGSIM/workflows/aMGSIM_ci/badge.svg)](https://github.com/genomewalker/aMGSIM/actions) [![PyPI](https://img.shields.io/pypi/v/aMGSIM)](https://pypi.org/project/aMGSIM/) [![Conda](https://img.shields.io/conda/v/genomewalker/aMGSIM)](https://anaconda.org/genomewalker/aMGSIM)


aMGSIM is an extension of [MGSIM](https://github.com/nick-youngblut/MGSIM/) to create simulated ancient metagenome reads in multiple microbial synthetic communities. It integrates the methods in [Gargammel](https://github.com/grenaud/gargammel) and provides flexibility to create different experimental scenarios. Furthermore, aMGSIM also can track damage over the codons of the predicted proteins from the microbial genomes. 



## INSTALLATION

We recommend having [**conda**](https://docs.conda.io/en/latest/) installed to manage the virtual environments

### Using pip

First, we create a conda virtual environment with:

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

aMGSIM uses as a starting point the `communities` subcommand from [MGSIM](https://github.com/nick-youngblut/MGSIM/) and then integrates three new subcommands:

- **ancient-genomes**: Estimate coverage, depth and other properties for each genome in each synthetic community
- **ancient-reads**: Simulate ancient reads for each taxon in each synthetic community
- **protein-analysis**: Tracking damage to the codon positions of each simulated read. 

You can access the list of commands by:

```bash
$ aMGSIM --list

Available Commands:
communities | ancient-genomes | ancient-reads | protein-analysis
```

First, we need to generate some synthetic communities. aMGSIM wraps the subcommand `communities` from [MGSIM](https://github.com/nick-youngblut/MGSIM/) to 

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


Once we have the synthetic communities' composition, we can create the ancient and modern partitions for each sample. We will use the subcommand `ancient-genomes`:

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

The subcommand takes as input a YAML [config](examples/config/ag-config.yml) file with the different parameters used to create the set of genomes in each synthetic community. A description of the different parameters can be found [here](examples/config/ag-config.yml). Some of the essential parameters are `genome_table` (the same used as in the subcommand `communities`); `abund_table` the abundance table created with the subcommand `communities` or we can define our values creating a table like:
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

And the `genome_comp` where we can define for each sample which genome will be ancient and the coverage. This file can be used to spike-in specific genomes at known coverages in the a sample. [Here](examples/data/genome-comp.tsv) you can find an example of the file. It is a TSV file with four columns: 

```
Taxon                                Community   Coverage   onlyAncient
Escherichia_coli_K-12_MG1655         1           1.0        True
Methanosarcina_barkeri_MS            2           0.5        False
Clostridium_perfringens_ATCC_13124   3             0        None
```

The column `onlyAncient` can take one of the following values:

- **True**: The genome in this sample will only contain ancient reads with the coverage specified, 1.0 in the example
- **False**: The genome in this sample will contain a mixture of modern and ancient reads. The ancient coverage will be the one specified in the file, 0.5 in the example. 
- **None**: The genome in this sample will only contain modern reads.

The values of coverage will be always downsized to fit the _maximum coverage allowed_ by the proportion of the taxon in the sample, the number of reads and the size of the genome. If `onlyAncient` is `True` and the `Coverage` value exceeds the _maximum coverage allowed_, the `Coverage` will be set to the _maximum coverage allowed_. In the case where `onlyAncient` is `False` and the `Coverage` value exceeds the _maximum coverage allowed_, the `Coverage` will be set to a random value defined by the limits defined in the config file by the `coverage` parameter. In case the random value is higher than the _maximum coverage allowed_ it will take the _maximum coverage allowed_ value.

The subcommand `ancient-genomes` also creates the fragment size distribution for each genome in each sample. At the moment aMGSIM only accepts `lognormal` as the distribution where to sample the fragment lengths. We can specify the mode of the distribution for modern (`mode-len-modern`) and ancient (`mode-len-ancient`) and aMGSIM will estimate the parameters to achieve an approximate distribution of the fragment lenghts. The estimates values will be used by `fragSim` downstream.

The output of this subcommand will produce a TSV file (_PREFIX_`_read-abundances.tsv`) with the number of reads that will be generated for each taxon in each sample; and a JSON file with the details to generate the synthetic data of each taxon in each sample.

```
{
    "data": [
        {
            "comm": 2,
            "taxon": "Clostridium_perfringens_ATCC_13124",
            "rel_abund": 1.7555391630000001,
            "genome_size": 3256683,
            "onlyAncient": true,
            "fragments_ancient": {
                "fragments": {
                    "length": [
                        30,
                        31,
                        32,
                        33,
                        34,
                        ...
                        115,
                        116,
                        117,
                        120,
                        122
                    ],
                    "freq": [
                        0.022599003395321643,
                        0.026239989701394066,
                        0.030409852333009715,
                        0.0337998101194531,
                        0.037248770335723824,
                        ...
                        2.1455429028122703e-06,
                        1.0727714514061352e-06,
                        2.1455429028122703e-06,
                        1.0727714514061352e-06,
                        2.1455429028122703e-06
                    ]
                },
                "dist_params": {
                    "mu": 3.740959125904414,
                    "sigma": 0.22820971011435498,
                    "rnd_seed": 27473
                },
                "avg_len": 43.92603991782571,
                "seq_depth": 87776,
                "seq_depth_original": 80713,
                "fold": 2.344874217109863,
                "fold_original": 2.1562061766527476
            },
            "fragments_modern": null,
            "coverage_enforced": false,
            "seq_depth": 87776
        }
    ],
    "experiment": {
        "library": "pe",
        "read_length": 150,
        "seqSys": "HS25",
        "n_reads": 5000000.0,
        "date": "2021-04-12 07:49:09"
    }
}
```



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
