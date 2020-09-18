#!/usr/bin/env python

"""
ancient-reads: simulating ancient reads

Usage:
  ancient-reads [options] <genome_table> <config>
  ancient-reads -h | --help
  ancient-reads --version

Options:
  <genome_table>      Taxon genome info.
  <config>            Config parameters
  -h --help           Show this screen.
  -d --debug                  Debug mode (no subprocesses; verbose output)
  --version           Show version.

Description:
  Simulating ancient reads for each taxon in each synthetic community

  abund_table
  -----------
  * tab-delimited
  * must contain 3 columns
    * "Community" = community ID (ie., sample ID)
    * "Taxon" = taxon name
    * "Perc_rel_abund" = percent relative abundance of the taxon

  genome_table
  ------------
  * tab-delimited
  * must contain 2 columns
    * "Taxon" = taxon name
    * "Fasta" = genome fasta file path
  * other columns are allowed

  ancient_genome_table
  ------------
  * tab-delimited
  * must contain 2 columns
    * "Taxon" = taxon name
    * "Fasta" = genome fasta file path
  * other columns are allowed

  art_config
  ------------
  * YAML with ART config parameters

  Output
  ------
  * A set of read files for each sample
    * directory structure: OUTPUT_DIR/COMMUNITY/ancient-read_files
    * read sequences are named by the taxon they originate from
"""

# import
# batteries
from collections import defaultdict
from docopt import docopt
import sys
import os
import re
import logging
from functools import partial
from multiprocessing import Pool
from schema import Schema, And, Optional, Or, Use, SchemaError

# application
from aMGSIM import SimReadsAncient
from MGSIM import SimReads
from aMGSIM.library import defaults as d
from aMGSIM.library import functions as f
import json
import pandas as pd
import pyfaidx
from pathlib import Path
import subprocess

# import tqdm
from tqdm.contrib.concurrent import process_map  # or thread_map
import tqdm
from dataclasses import dataclass
from Bio import SeqIO
import itertools
import gzip
from mimetypes import guess_type
from biolib.external.prodigal import Prodigal


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


def main(args):

    global debug
    if args["--debug"]:
        debug = True
    else:
        debug = None

    sys.excepthook = exceptionHandler

    # simulating reads
    args = f.validate_schema(args, d.schema_init_ar, debug)

    config = f.get_config(
        config=args["<config>"], schema=d.ar_schema_config, debug=debug
    )


def opt_parse(args=None):
    if args is None:
        args = docopt(__doc__, version="0.1")
    else:
        args = docopt(__doc__, version="0.1", argv=args)
    main(args)


# %%
