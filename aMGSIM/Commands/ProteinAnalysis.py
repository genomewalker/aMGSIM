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
  -n --procs=<p>         processes
                         [default: 1]
  -m --min-len=<l>        Minimum ORF size
                         [default: 30]
  -t --tmp=<d>           Tmp dir
                         [default: .tmp]
  -o --out-dir=<d>           Output directory
                         [default: protein-analysis]
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
from biolib.external.prodigal import Prodigal
from itertools import product
import glob
import multiprocessing.pool
from docopt import docopt
import logging
import pandas as pd
from functools import partial
from pathlib import Path
import sys
import tqdm

# application
from aMGSIM.library import pa_functions as pa

import tqdm
import os
import json
from aMGSIM import __version__

# Codon functions


debug = None


# From
class NoDaemonProcess(multiprocessing.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass


class NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess


# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
    def __init__(self, *args, **kwargs):
        kwargs["context"] = NoDaemonContext()
        super(MyPool, self).__init__(*args, **kwargs)


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


def combine_files(x, tmp_damage, out_dir):
    out_suffix = ".tsv.gz"
    fname = "{}_aa-damage".format(x)
    outfile = Path(out_dir, fname).with_suffix(out_suffix)
    files = glob.glob(str(Path(tmp_damage, x + "*")))
    li = []

    for file in files:
        df = pd.read_csv(file, index_col=None, header=0, sep="\t")
        li.append(df)

    df = pd.concat(li, axis=0, ignore_index=True)
    df.to_csv(
        path_or_buf=outfile,
        sep="\t",
        header=True,
        index=False,
        compression="gzip",
    )


def obj_dict(obj):
    return obj.__dict__


# Read files and get coordinates


def opt_parse(args=None):
    version = "Version: " + __version__
    if args is None:
        args = docopt(__doc__, version=version)
    else:
        args = docopt(__doc__, version=version, argv=args)
    main(args)


def main(args):
    global debug
    if args["--debug"]:
        debug = True
    else:
        debug = None

    sys.excepthook = exceptionHandler

    nproc = int(args["--procs"])
    ncpu = int(args["--cpus"])
    min_len = int(args["--min-len"])

    with open(args["<files>"], "r") as json_file:
        filename = json_file
        files = json.load(json_file)

    genome_files = files["genomes"]
    read_files = files["reads"]

    comms = list(files["reads"].keys())

    tmp_dir = args["--tmp"]
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)

    tmp_damage = os.path.join(tmp_dir, "damage-aa")
    if not os.path.isdir(tmp_damage):
        os.makedirs(tmp_damage, exist_ok=True)

    out_dir = args["--out-dir"]
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    p_procs = nproc * ncpu

    if p_procs > len(genome_files):
        p_procs = len(genome_files)

    logging.info("Predicting genes from genomes usin {} processes...".format(p_procs))
    output_dir = os.path.join(tmp_dir, "gene_prediction")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    prodigal = Prodigal(cpus=p_procs, verbose=True)
    gene_preds = prodigal.run(
        genome_files=genome_files,
        output_dir=output_dir,
        called_genes=False,
        translation_table=None,
        meta=False,
        closed_ends=False,
    )
    gene_pred = {}

    for k in gene_preds.keys():
        gene_pred[k] = {
            "faa": gene_preds[k].aa_gene_file,
            "fna": gene_preds[k].nt_gene_file,
            "translation_table": gene_preds[k].best_translation_table,
        }

    genome_ids = list(gene_pred.keys())

    func = partial(
        pa.analyze_proteins,
        files=files,
        gene_predictions=gene_pred,
        min_len=min_len,
        outdir=tmp_damage,
        debug=debug,
        nproc=nproc,
    )

    logging.info("Finding damage in codons...")
    comm_files = list(product(comms, genome_ids))

    # if p_procs > len(comm_files):
    #     p_procs = len(comm_files)

    if debug is True:
        data = list(map(func, comm_files))
    else:
        p = MyPool(ncpu)
        data = list(
            tqdm.tqdm(
                p.imap_unordered(func, comm_files),
                total=len(comm_files),
            )
        )
    logging.info("Combining files...")

    p_procs = nproc * ncpu
    if p_procs > len(comms):
        p_procs = len(comms)

    func = partial(combine_files, tmp_damage=tmp_damage, out_dir=out_dir)
    if debug is True:
        ofiles = list(map(func, comms))
    else:
        p = MyPool(p_procs)
        ofiles = list(
            tqdm.tqdm(
                p.imap_unordered(func, comms),
                total=len(comms),
            )
        )
    # for comm in comms:
    #     out_suffix = ".tsv.gz"
    #     fname = "{}_aa-damage".format(comm)
    #     outfile = Path(out_dir, fname).with_suffix(out_suffix)
    #     files = glob.glob(str(Path(tmp_damage, comm + "*")))
    #     li = []

    #     for file in files:
    #         df = pd.read_csv(file, index_col=None, header=0, sep="\t")
    #         li.append(df)

    #     df = pd.concat(li, axis=0, ignore_index=True)
    #     df.to_csv(
    #         path_or_buf=outfile,
    #         sep="\t",
    #         header=True,
    #         index=False,
    #         compression="gzip",
    #     )

    logging.info("Protein analysis done.")
