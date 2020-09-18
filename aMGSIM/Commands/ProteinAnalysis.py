#!/usr/bin/env python

"""
protein-analysis: simulating ancient reads

Usage:
  protein-analysis [options] <files>
  protein-analysis -h | --help
  protein-analysis --version

Options:
  <files>                read files
  -p --cpus=<p>          cpus
                         [default: 1]
  -m --minlen=<l>        Minimum ORF size
                         [default: 30]
  -t --tmp=<d>           Tmp dir
                         [default: .tmp]
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
import re
import pyranges as pr
import gzip
from mimetypes import guess_type
from biolib.external.prodigal import Prodigal
from Bio import SeqIO

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
from aMGSIM.library import defaults as d
from aMGSIM.library import functions as f
from aMGSIM.library import pa_functions as pa

import jsonpickle
import datetime
import tqdm
import os
import json

# Codon functions


debug = None


def exceptionHandler(
    exception_type, exception, traceback, debug_hook=sys.__excepthook__
):
    """Print user friendly error messages normally, full traceback if DEBUG on.
    Adapted from http://stackoverflow.com/questions/27674602/hide-traceback-unless-a-debug-flag-is-set
    """
    if debug:
        print("\n*** Error:")
        # raise
        debug_hook(exception_type, exception, traceback)
    else:
        print("{}: {}".format(exception_type.__name__, exception))


logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


# logging
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def obj_dict(obj):
    return obj.__dict__


# Read files and get coordinates


def opt_parse(args=None):
    if args is None:
        args = docopt(__doc__, version="0.1")
    else:
        args = docopt(__doc__, version="0.1", argv=args)
    main(args)


def main(args):
    global debug
    if args["--debug"]:
        debug = True
    else:
        debug = None

    sys.excepthook = exceptionHandler

    nproc = int(args["--cpus"])

    with open(args["<files>"], "r") as json_file:
        filename = json_file
        files = json.load(json_file)

    genomes = files["genomes"]

    tmp_dir = args(["--tmp"])

    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)

    logging.info("Predicting genes from genomes...")
    output_dir = os.path.join(tmp_dir, "gene_prediction")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    prodigal = Prodigal(cpus=nproc, verbose=True)
    genes = prodigal.run(
        genome_files=genomes,
        output_dir=output_dir,
        called_genes=False,
        translation_table=None,
        meta=False,
        closed_ends=False,
    )

    deam_file = "/vol/cloud/antonio/projects/anc-com-sim/sandbox/test_data/test_out3/1/comm-1_deamSim.fa"
    fragsim_file = "/vol/cloud/antonio/projects/anc-com-sim/sandbox/test_data/test_out3/1/comm-1_fragSim.fa"
    fasta_file = "/vol/cloud/antonio/projects/anc-com-sim/sandbox/test_data/test_out3/genomes/GCA_000007185.fasta"

    # Genes
    fna = "/vol/cloud/antonio/projects/anc-com-sim/sandbox/test_data/.tmp/gene_prediction/GCA_000007185_genes.fna"
    faa = "/vol/cloud/antonio/projects/anc-com-sim/sandbox/test_data/.tmp/gene_prediction/GCA_000007185_genes.faa"

    recs = []
    encoding = guess_type(deam_file)[1]  # uses file extension
    _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open

    pattern = re.compile("(\S+)___(\S+)---(\d+):(\S+):([+\-]):(\d+):(\d+):(\d+):(.*)")
    with _open(deam_file) as f:
        records = enumerate(SeqIO.parse(f, "fasta"))
        recs = pa.get_headers(records=records, pattern=pattern)
        df = pd.DataFrame(recs)

    # Get reads
    df_reads = pr.PyRanges(df)

    # Get genes info
    genes = pa.get_gene_coordinates(fna=fna)

    # Get intersections
    # Will rop intersects that are too short for coding AAs
    r2g_intersections = pa.get_intersections(
        reads=df_reads, genes=genes, genome=fasta_file, min_len=min_len
    )

    names = {"Start_b": "Start_gene", "End_b": "End_gene", "Strand": "Strand_gene"}
    r2g_intersections = r2g_intersections.join(genes.drop("type"), report_overlap=True)
    r2g_intersections = pr.PyRanges(
        r2g_intersections.df[
            r2g_intersections.df["Overlap"] == r2g_intersections.df["intersect_length"]
        ].rename(columns=names)
    )
    r2g_intersections_df = r2g_intersections.df
    gn = r2g_intersections_df["Name"]
    read_multi_span = r2g_intersections_df[gn.isin(gn[gn.duplicated()])].sort_values(
        "Name"
    )
    # len(read_multi_span.index)

    # Get sequence from interval
    r2g_intersections.nondamaged_seq = pa.get_nondamaged_seqs(
        intersections=r2g_intersections, genome=fasta_file
    )

    r2g_intersections = pa.get_damaged_seqs(
        intersections=r2g_intersections, deam_file=deam_file
    )

    r2g_intersections.damage_intersection = r2g_intersections.df.apply(
        pa.find_damage, axis=1
    )

    aa_damage = r2g_intersections.df.merge(
        pa.fasta_to_dataframe(fna)[["name", "sequence"]].rename(
            columns={"name": "gene_name"}
        )
    ).apply(pa.get_seqs_inframe, axis=1)
    aa_damage = aa_damage.reset_index()

    # test[['Strand_read','damage_aaseq_diffs', 'damage_codon_diffs']]
    # test[~test.damage_aaseq_diffs.isnull()][['Strand_read','damage_aaseq_diffs', 'damage_codon_diffs']]

    # test_aa = test[['Name', 'gene_name', 'intersect_seq_inframe_aa', 'Strand_read', 'Strand_gene']].merge(
    #     fasta_to_dataframe(faa)[["name", "sequence"]].rename(columns={"name": "gene_name", "sequence": "sequence_aa"}))
    # test_aa_filt = test[test.apply(
    #     lambda x: x.intersect_seq_inframe_aa not in x.sequence_aa, axis=1)].reset_index()
