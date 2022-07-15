# aMGSIM

[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/aMGSIM?include_prereleases&label=version)](https://github.com/genomewalker/aMGSIM/releases) [![aMGSIM](https://github.com/genomewalker/aMGSIM/workflows/aMGSIM_ci/badge.svg)](https://github.com/genomewalker/aMGSIM/actions) [![PyPI](https://img.shields.io/pypi/v/aMGSIM)](https://pypi.org/project/aMGSIM/) [![Conda](https://img.shields.io/conda/v/genomewalker/aMGSIM)](https://anaconda.org/genomewalker/aMGSIM)


aMGSIM is an extension of [MGSIM](https://github.com/nick-youngblut/MGSIM/) to create simulated ancient metagenome reads in multiple microbial synthetic communities. It integrates the methods in [Gargammel](https://github.com/grenaud/gargammel) and provides flexibility to create different experimental scenarios. Furthermore, aMGSIM also can track damage over the codons of the predicted proteins from the microbial genomes. 

aMGSIM can produce random metagenomes (A) or it can estimate the parameters from a real sample (B) using [bam-filter](https://github.com/genomewalker/bam-filter) and [metaDMG](https://metadmg-dev.github.io/metaDMG-core/)

<img src='extra/aMGSIM-wf.png' width='100%'>

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

aMGSIM incorporates three new subcommands to generate synthetic metagenomic reads:

- **estimate**: Estimate coverage, depth and other properties for each genome in sample processed with filterBAM
- **ancient-genomes**: Estimate coverage, depth and other properties for each genome in each synthetic community
- **ancient-reads**: Simulate ancient reads for each taxon in each synthetic community
- **protein-analysis**: Tracking damage to the codon positions of each simulated read. 

You can access the list of commands by:

```bash
$ aMGSIM --list

Available Commands:
communities | estimate | ancient-genomes | ancient-reads | protein-analysis
```
We have two different approaches to generate the synthetic communities:

1. We can use empirical results from [filterBAM]() and [metaDMG](#) as shown [here](#) with the subcommand `estimate`.
2. We can use the subcommand `communities` from [MGSIM](https://github.com/nick-youngblut/MGSIM/)


> We recommend to use the estimate subcommand


To use the empirical data from filterBAM and metaDMG, you will need the following files:

- A file containing the paths to the genomes used for mapping the reads. A TSV file with the `name` and `location of the genomes` fasta files like the file [genome-paths-list.tsv](examples/data/genome-paths-list.tsv) 
- A file containing the mapping stats after filtering the BAM files with [bam-filter](https://github.com/genomewalker/bam-filter).
- The metaDMG CSV ouput for the sample and the misincorporation file.
- The taxonomy dump files (names and nodes) used for filterBAM and metaDMG.


Then you can use the command as follows:

```bash
aMGSIM estimate f2b-config.yaml
```

And example of the config file can be found here [**aw-config.yaml**](examples/conifg/aw-config.yaml) associated to the files described above are:

```yaml
# TSV file with the name of the genome and the fasta location of the references used in filterBAM, the file has no header and has the following format: name\tfasta_location
genome_paths: genome-paths-list.tsv
# Nodes and names dumps used for the filterBAM analysis
names_dmp: names-mdmg.dmp
nodes_dmp: nodes-mdmg.dmp
# Sample name
sample_name: "0e59dd2155"
# Statistics of the mapping results from the filterBAM analysis
filterBAM_stats: 0e59dd2155.filterBAM.dedup_stats-filtered.tsv.gz
# Filter options for the filetrBAM table
filterBAM_filter_conditions: { "breadth_exp_ratio": 0.5, "coverage_mean": 0.1 }
# filterBAM abundance table
filterBAM_profiling_results: 0e59dd2155.filterBAM-none.tsv.gz
# Results of the metaDMG analysis
mdmg_results: 0e59dd2155.filterBAM-mdmg.weight-0.csv.gz
# Filter conditions to identify damaged references
mdmg_filter_conditions: { "d_max": 0.1, "phi": 500, "q": 0.25 }
# Taxonomic scope of the metaDMG analysis
taxonomic_rank: "species"
# Number of non-damaged genomes
max_genomes_nondamaged: 10
# How to select the non-damaged genomes: random, most_abundant, least_abundant
max_genomes_nondamaged_selection: "random"
# Number of damaged genomes
max_genomes_damaged: 10
# How to select the non-damaged genomes: random, most_abundant, least_abundant
max_genomes_damaged_selection: "random"
# Should we reestimate the abundances after the genome selection?
use_restimated_proportions: True
# Number of cpus. [Default: 1]
cpus: 16

```
This will create all the files need for the next step:

- **0e59dd2155.communities.tsv**: This file contains the estimated abundance table for the synthetic community. Columns are `Community`, `Taxon`, `Rank` and `Perc_rel_abund`
- **0e59dd2155.filepaths.tsv**: This file contains the file paths for the reference genomes used in the filterBAM analysis. Columns are `Taxon`, `TaxId`, `Accession` and `Fasta`
- **0e59dd2155.genome-compositions.tsv**: This file is the one controlling how the genomes are going to be processed and will contain the empirical values. We can change them to spike certain taxa to very high or very low coverage values. Anything added to this file (all genomes by default) will avoid the random generation of coverage, sequencing depth and other vlaues. Columns are `Taxon`, `Community`, `Coverage` and `onlyAncient`
  

If you want to use the `communities` subcommand:

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

When these files are in place we can call the second command that will generate the synthetic community. This command will use the config file `ag-config.yaml` will set the basic parameters to generate the fragment needed for the simulations. Here is where you need to decide if we use paired/single end, which machine, read length, modal length of the damaged/non-damaged fraction and the libprep method. We will use the subcommand `ancient-genomes`:

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

The output of this subcommand will produce a TSV file (_PREFIX_`_read-abundances.tsv`) with the number of reads that will be generated for each taxon in each sample; and a JSON file with the details to generate the synthetic data of each taxon in each sample. The JSON produced will be used by the subcommand `ancient-reads` and has the following structure:

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

The subcommand `ancient-reads` is the one creating all the fasta and fastq files with the reads for the genomes in each community:

```
ancient-reads: Simulate ancient reads for each taxon in each synthetic
               community

Usage:
  ancient-reads [options] <genome_table> <config>
  ancient-reads -h | --help
  ancient-reads --version

Options:
  <genome_table>      Taxon genome info.
  <config>            Config parameters
  -h --help           Show this screen.
  -d --debug          Debug mode (no subprocesses; verbose output)
  --version           Show version.

Description:
  Simulating ancient reads for each taxon in each synthetic community

  genome_table
  ------------
  * tab-delimited
  * must contain 2 columns
    * "Taxon" = taxon name
    * "Fasta" = genome fasta file path
  * other columns are allowed

  config
  ------
  * YAML with config parameters 
    (Check https://github.com/genomewalker/aMGSIM for details)

  Output
  ------
  * A set of read files for each sample (fragSim, deamSim, ART)
    - read sequences are named by the taxon they originate from
    - directory structure: OUTPUT_DIR/COMMUNITY/ancient-read_files
  * A JSON file with the location of each file type (fragSim, deamSim, ART)
  ```

It takes as input the file with the [genome information](examples/data/genome_list.txt) used by the subcommand `communities` and a config file that specifies the details related to the tools in gargammel. A description of the different parameters can be found [here](examples/config/ar-config.yml). The config file has five main categories:
- **global**: Here one can configure general parameters like the location of the executables, the number of cpus to use or the library preparation
- **fragSim**: Here one can specify the parameters for fragSim and where one can use the JSON produced by the `ancient-genomes` subcommand 
- **deamSim**: One can define the propertied for the damage, it can be a `misincorporation` file from mapDamage or the parameters for the [Briggs](https://www.pnas.org/content/104/37/14616) model
- **adptSim**: Where to set the adapters to be used
- **art**: One can specify custom error models in case they are not included in [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)


The ouput of `ancient-reads` are the read files for each community after each step (fragSim -> deamSim -> adptSim -> ART) and a JSON file with the location of each file. This JSON file can be used for the last subcommand in aMGSIM, `protein-analysis`:

```
protein-analysis: Tracking damage to the codon positions of each simulated read

Usage:
  protein-analysis [options] <files>
  protein-analysis -h | --help
  protein-analysis --version

Options:
  <files>                Read files
  -p --cpus=<p>          Cpus
                         [default: 1]
  -n --procs=<p>         processes
                         [default: 1]
  -m --min-len=<l>       Minimum ORF size
                         [default: 30]
  -t --tmp=<d>           Tmp dir
                         [default: .tmp]
  -o --out-dir=<d>       Output directory
                         [default: protein-analysis]
  -d --debug             Debug mode (no subprocesses; verbose output)
  -h --help              Show this screen
  --version              Show version

Description:
  Tracking damage to the codon positions of each simulated read. It predicts
  genes with Prodigal. Only for Archaea/Bacteria/Virus

  files
  -----
  * A JSON file generated by the subcommand ancient-reads
  
  Output
  ------
  * A TSV file with the damage information at the codon level
```

This subcommand will predict the genes in the genome with [Prodigal](https://github.com/hyattpd/Prodigal) and will match the reads and damages to the codon positions in the predicted genes. Depending the number of reads simulated it can take some time to produce all the results. One can run the subcommand like this:

```
aMGSIM protein-analysis --cpus 3 --procs 16 example_abund_read-files.json
```

Here, `protein-analysis` will start 3 jobs and use 16 cores per each job, and will produce a rich TSV with the following columns:

- **index**: 14
  - Read numnber
- **chromosome**: Methanosarcina_barkeri_MS__seq-0
  - Genome/contig where the read and gene are originated
- **End**: 2114466
  - End position of the read/gene intersection
- **End_gene**: 2114583
  - End position of the gene in the genome/contig
- **End_read**: 2114466
  - End position of the read in the genome/contig
- **read_name**: 1___Methanosarcina_barkeri_MS__seq-0---1027401:ancient:+:2114436:2114466:30:4
  - Name of the read
- **Overlap**: 30
  - Overlap between the gene and the read
- **Start**: 2114436
  - Start position of the read/gene intersection
- **Start_gene**: 2114416
  - Start position of the gene in the genome/contig
- **Start_read**: 2114436
  - Start position of the read in the genome/contig
- **Strand_gene**: +
  - Strand of the predicted gene
- **Strand_read**: +
  - Strand of the simulated read
- **damage**: 4
  - Damage position in the read
- **damage_aaseq_diffs**: 1:L>F
  - Amino acid differences owing to the damage (Original>Damaged)
- **damage_codon_diffs**: 3:CTT>TTT:L>F
  - Codon and amino acid differences owing to the damage (Original>Damaged)
- **damage_inframe_ag**: 
  - Damage position in the open reading frame for A-to-G (0-index)
- **damage_inframe_ct**: 3
  - Damage position in the open reading frame for C-to-T (0-index)
- **damage_intersection**: 4
  - Damage position in the intersection between read/gene in frame (1-index)
- **damaged_seq**: AACTTTTACAAATGTGAGATAAAAAATAGT
  - Damage sequence in present in read
- **damaged_seq_inframe_aa**: NFYKCEIKN
  - Damaged amino acid sequence (in frame)
- **damaged_seq_inframe_nt**: AACTTTTACAAATGTGAGATAAAAAAT
  - Damaged nucleotide sequence (in frame)
- **damaged_seq_length**: 30
  - Length of the damaged sequence
- **gene_name**: Methanosarcina_barkeri_MS__seq-0_1774
  - Name of the predicted gene
- **intersect_length**: 30
  - Length of the read/gene intersection
- **intersect_seq**: AACCTTTACAAATGTGAGATAAAAAATAGT
  - Sequence of the intersection
- **intersect_seq_inframe_aa**: NLYKCEIKN
  - Amino acid sequence of the intersection without damage
- **intersect_seq_inframe_nt**: AACCTTTACAAATGTGAGATAAAAAAT
  - Nucleotide sequence of the intersection without damage
- **nondamaged_seq**: AACCTTTACAAATGTGAGATAAAAAATAGT
  - Nucleotide non-damaged sequence
- **read_length**: 30
  - Read length
- **sequence**: ATGCCTGAAGTAGAAGTCAAGAACCTTTACAAATGTGAGATAAAAAATAGTGAGATAAAAAATAGTGAGATAAAAAATAGTGAGATAAAAAATAGTGAGATAAAAAATAGTGAGATAAAAAATAGTGAGATAAAAAATAGTGAGATAAAAAATAGGTCTAAAAATTGA
  - Gene nucletide sequence
- **sequence_aa**: MPEVEVKNLYKCEIKNSEIKNSEIKNSEIKNSEIKNSEIKNSEIKNSEIKNRSKN*
  - Gene amino acid sequence
- **type**: ancient
  - Type of read
- **community**: 1
  - Community where it belongs


## LICENSE
See [LICENSE](./LICENSE)

## REFERENCE
Please cite the original [MGSIM](https://github.com/nick-youngblut/MGSIM/)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3696891.svg)](https://doi.org/10.5281/zenodo.3696891)
