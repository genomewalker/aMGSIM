#!/usr/bin/env python

"""
protein-analysis: simulating ancient reads

Usage:
  protein-analysis [options] <genomes> <files>
  protein-analysis -h | --help
  protein-analysis --version

Options:
  <genomes>              Config parameters
  <files>                read files
  -p --cpus=<p>          cpus
                         [default: 1]
  -d --debug             Debug mode (no subprocesses; verbose output)
  -h --help              Show this screen.
  --version              Show version.

Description:
  Simulating ancient reads for each taxon in each synthetic community

  abund_table
  -----------
  * tab-delimited
  * must contain 3 columns
    * "Community" = community ID (ie., sample ID)
    * "Taxon" = taxon name
    * "Perc_rel_abund" = percent relative abundance of the taxon

  Output
  ------
  * A JSON file with the properties of the selected ancient genomes
    * directory structure: OUTPUT_DIR/COMMUNITY/ancient-read_files
    * read sequences are named by the taxon they originate from
"""
# %%
# import
# batteries
from docopt import docopt
import logging
import pandas as pd
from functools import partial
from multiprocessing import Pool
import numpy as np
import json
from pathlib import Path
import sys
# application
from MGSIM import SimReads
from aMGSIM.library import defaults as d
from aMGSIM.library import functions as f
import jsonpickle
import datetime
import tqdm
import os
from biolib.external.prodigal import Prodigal, ProdigalGeneFeatureParser
import Bio.Data.CodonTable
from itertools import product
import json

# Codon functions

# codons = Bio.Data.CodonTable.standard_dna_table.forward_table
# codons_stop = Bio.Data.CodonTable.standard_dna_table.stop_codons
# codons_stop = {el: '*' for el in codons_stop}
# codons = {**codons, **codons_stop}

# # Get codons that can get damaged
# # C -> T

# filtered_dict = {k: v for (k, v) in codons.items() if "C" in k}
# codon_list = (list(filtered_dict.keys()))

# d = {
#     "A": "A",
#     "T": "T",
#     "G": "G",
#     "C": ["T", "C"],
#     }

# d_rev = {
#     "T": "T",
#     "C": "C",
#     "A": "A",
#     "G": ["G", "A"],
#     }

# codons_damage = {}
# for k in codon_list:
#     c_list = list(map("".join, product(*map(d.get, k))))
#     c_list.remove(str(k))
#     codons_damage[k] = {"aa": codons[str(k)],
#                         "mods": {el: codons[str(el)] for el in c_list}
#                         }


# # A -> G
# filtered_dict_rev = {k: v for (k, v) in codons.items() if "G" in k}
# codon_list_rev = (list(filtered_dict_rev.keys()))
# codons_damage_rev = {}
# for k in codon_list_rev:
#     c_list = list(map("".join, product(*map(d_rev.get, k))))
#     c_list.remove(str(k))
#     codons_damage_rev[k] = {"aa": codons[str(k)],
#                             "mods": {el: codons[str(el)] for el in c_list}
#                             }
# print(json.dumps(codons_damage_rev, indent=4))


debug = None


def exceptionHandler(exception_type, exception, traceback, debug_hook=sys.__excepthook__):
    '''Print user friendly error messages normally, full traceback if DEBUG on.
       Adapted from http://stackoverflow.com/questions/27674602/hide-traceback-unless-a-debug-flag-is-set
    '''
    if debug:
        print("\n*** Error:")
        # raise
        debug_hook(exception_type, exception, traceback)
    else:
        print("{}: {}".format(exception_type.__name__, exception))


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def obj_dict(obj):
    return obj.__dict__

# Read files and get coordinates


def main(args):
    global debug
    if args['--debug']:
        debug = True
    else:
        debug = None

    sys.excepthook = exceptionHandler
        
    nproc = int(args['--cpus'])
    genome_table = SimReads.load_genome_table(args['<genomes>'])
    with open(args['<files>'], "r") as json_file:
        filename = json_file
        files = json.load(json_file)

    tmp_dir = ".tmp"

    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)

    logging.info('Predicting genes from genomes...')
    output_dir = os.path.join(tmp_dir, "gene_prediction")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # func = partial(predict_genes,
    #                cpus=1,
    #                verbose = True,
    #                output_dir=output_dir,
    #                called_genes=False,
    #                translation_table=None,
    #                meta=False,
    #                closed_ends=False)

    # if debug is True:
    #     gene_predictions = list(map(func, genomes))
    # else:
    #     p = Pool(nproc)
    #     gene_predictions = list(
    #         tqdm.tqdm(p.imap_unordered(func, genomes), total=len(genomes)))
    #     # ancient_genomes_data = p.map(func, genomes)

    prodigal = Prodigal(cpus=nproc, verbose=True)
    genes = prodigal.run(
        genome_files=set(genome_table.Fasta.values),
        output_dir=output_dir,
        called_genes=False,
        translation_table=None,
        meta=False,
        closed_ends=False)

    gffs = {}
    for i in genes:
        print(i)
        gffs[i] = ProdigalGeneFeatureParser(genes[i].gff_file)

#%%
import pyranges as pr
import pandas as pd
from Bio import SeqIO
from mimetypes import guess_type
import itertools
import gzip
import re

gff_file = "/vol/cloud/antonio/projects/anc-com-sim/sandbox/test_data/.tmp/gene_prediction/GCA_000007185.1.gff"
gff = pr.read_gtf(gff_file, as_df = False)

deam_file = "/vol/cloud/antonio/projects/anc-com-sim/sandbox/test_data/test_out3/1/comm-1_deamSim.fa"
recs = []
encoding = guess_type(deam_file)[1]  # uses file extension

_open = partial(
                gzip.open, mode='rt') if encoding == 'gzip' else open

pattern = re.compile("(\S+)__(\S+)--(\d+):(\S+):([+\-]):(\d+):(\d+):(\d+):(.*)")
with _open(deam_file) as f:
    records = enumerate(SeqIO.parse(f, 'fasta'))

    recs = []
    for i, record in records:
        m1 = re.match(pattern, record.id)
        read = {"Chromosome"   | Start     | End       | Name                                   | Score     | Strand 
            
        }
    print(recs[0])

#%%    
    print(gffs)
    #print(json.dumps(gffs, indent=4))

def opt_parse(args=None):
    if args is None:
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)

# %%
