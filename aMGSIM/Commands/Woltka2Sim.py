#!/usr/bin/env python

"""
woltka2sim: Estimate coverage, depth and other properties for each genome
                 in sample processed with WOLTKA

Usage:
  woltka2sim [options] <config>
  woltka2sim -h | --help
  woltka2sim --version

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
"""

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
import re

# application
from aMGSIM.library import defaults as d
from aMGSIM.library import functions as f
from aMGSIM.library import woltka_functions as w
import datetime
import tqdm
from collections import OrderedDict
from aMGSIM import __version__
import taxopy as txp

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


def generate_genome_compositions(df, taxdb, taxonomic_rank, sample_name):
    """Function to generate a genome composition file compatible with aMGSIM.

    Args:
        df (pandas.DataFrame): A dataframe with the genome information.
        taxdb (dict): A dict with the taxonomic information
        taxonomic_rank (str): Selected taxonomic scope
        sample_name (str): Sample name

    Returns:
        pandas.DataFrame: A table with the genome information.
    """
    df_selected = df[[taxonomic_rank, "coverage_mean", "is_damaged"]].copy()
    df_selected["Taxon"] = df_selected[taxonomic_rank].apply(
        lambda x: re.sub("[^0-9a-zA-Z]+", "_", taxdb.taxid2name[x])
    )
    df_selected.loc[:, "Community"] = sample_name
    df_selected = df_selected[["Taxon", "Community", "coverage_mean", "is_damaged"]]
    df_selected.columns = ["Taxon", "Community", "Coverage", "onlyAncient"]
    # convert onlyAncient column to string
    df_selected.loc[:, "onlyAncient"] = df_selected["onlyAncient"].astype(str)

    return df_selected


def generate_community_file(
    df, sample_name, taxonomic_rank, taxdb, use_restimated_proportions
):
    """Function to create a community file compatible with aMGSIM.

    Args:
        df (pandas.DataFrame): A dataframe with the genome information.
    """
    # create a new DataFrame with the selected taxa
    if use_restimated_proportions:
        df_selected = df[[taxonomic_rank, "proportion_new"]].copy()
        df_selected.loc[:, "Community"] = sample_name
        df_selected["Taxon"] = df_selected[taxonomic_rank].apply(
            lambda x: re.sub("[^0-9a-zA-Z]+", "_", taxdb.taxid2name[x])
        )
        df_selected["Rank"] = np.arange(len(df_selected)) + 1
        df_selected = df_selected[["Community", "Taxon", "Rank", "proportion_new"]]
        df_selected.columns = ["Community", "Taxon", "Rank", "Perc_rel_abund"]

    else:
        df_selected = df[[taxonomic_rank, "proportion"]].copy()
        df_selected.loc[:, "Community"] = sample_name
        df_selected["Taxon"] = df_selected[taxonomic_rank].apply(
            lambda x: re.sub("[^0-9a-zA-Z]+", "_", taxdb.taxid2name[x])
        )
        df_selected["Rank"] = np.arange(len(df_selected)) + 1
        df_selected = df_selected[["Community", "Taxon", "Rank", "proportion"]]
        df_selected.columns = ["Community", "Taxon", "Rank", "Perc_rel_abund"]

    return df_selected


def generate_genome_paths(df, genome_paths, taxdb, taxonomic_rank):
    """Function to generate a genome paths file compatible with aMGSIM.

    Args:
        df (pandas.DataFrame): A dataframe with the genome information.
        genome_paths (pandas.DataFrame): A dataframe with the path information.
        taxdb (dict): A dictionary with the taxonomic information
        taxonomic_rank (str): Selected taxonomic scope

    Returns:
        pandas.DataFrame: A dataframe with the genome paths files
    """
    # Taxon	Accession	Fasta
    df_selected = (
        df[["reference", taxonomic_rank]].copy().merge(genome_paths, on="reference")
    )
    df_selected["Taxon"] = df_selected[taxonomic_rank].apply(
        lambda x: re.sub("[^0-9a-zA-Z]+", "_", taxdb.taxid2name[x])
    )
    df_selected = df_selected[["Taxon", taxonomic_rank, "reference", "path"]]
    df_selected.columns = ["Taxon", "TaxId", "Accession", "Fasta"]
    return df_selected


def select_taxa(df, is_damaged=False, mode="most_abundant", n_taxa=10):
    """Function to select taxa to be simulated.

    Args:
        df (pandas.DataFrame): [description]
        is_damaged (bool, optional): [description]. Defaults to False.
        mode (str, optional): [description]. Defaults to 'most_abundant'.
        n_taxa (int, optional): [description]. Defaults to 100.

    Raises:
        ValueError: Fails if the mode is not valid.
    """
    if mode == "most_abundant":
        if is_damaged:
            df = df[df["is_damaged"] != True].sort_values(by="rank", ascending=True)
        else:
            df = df[df["is_damaged"] == True].sort_values(by="rank", ascending=True)
        # select n rows from pandas DataFrame
        df = df.head(n_taxa)
    elif mode == "least_abundant":
        if is_damaged:
            df = df[df["is_damaged"] != True].sort_values(by="rank", ascending=False)
        else:
            df = df[df["is_damaged"] == True].sort_values(by="rank", ascending=False)
        df = df.head(n_taxa)
    elif mode == "random":
        if is_damaged:
            df = df[df["is_damaged"] != True].sample(n=n_taxa)
        else:
            df = df[df["is_damaged"] == True].sample(n=n_taxa)
    else:
        raise ValueError(
            'mode should be one of "most_abundant", "random", "least_abundant"'
        )
    return df


def create_output_files(prefix):
    # create output files
    out_files = {
        "paths": f"{prefix}.filepaths.tsv",
        "compositions": f"{prefix}.genome-compositions.tsv",
        "communities": f"{prefix}.communities.tsv",
    }
    return out_files


def main(args):
    # simulating reads
    global debug
    if args["--debug"]:
        debug = True
    else:
        debug = None

    sys.excepthook = exceptionHandler

    # simulating reads
    args = f.validate_schema(args, d.schema_init_w, debug)
    config = f.get_config(
        config=args["<config>"], schema=d.w_schema_config, debug=debug
    )
    taxonomic_rank = config["taxonomic_rank"]
    # load tables
    logging.info("Loading genome paths file...")
    genome_paths = w.load_genome_paths(config["genome_paths"])

    # load taxonomy
    taxdb = txp.TaxDb(
        nodes_dmp=config["nodes_dmp"], names_dmp=config["names_dmp"], keep_files=True
    )

    logging.info("Loading woltka stats file...")
    woltka_stats = w.load_woltka_stats_table(config["woltka_stats"])
    woltka_stats_n = woltka_stats.shape[0]
    woltka_stats = w.filter_woltka_taxa(
        df=woltka_stats,
        filter_conditions=config["woltka_filter_conditions"],
    )
    logging.info(f"Kept {woltka_stats.shape[0]} taxa from {woltka_stats_n}")

    logging.info("Loading woltka profiling results...")
    woltka_profiling_results = w.load_woltka_profiling_table(
        config["woltka_profiling_results"]
    )
    woltka_profiling_results = woltka_profiling_results[
        woltka_profiling_results["reference"].isin(
            woltka_stats["reference"].unique().tolist()
        )
    ]
    logging.info("Loading metaDMG results...")
    mdmg_results = w.load_mdmg_results(config["mdmg_results"])
    # find which taxon are damaged
    damaged_taxa = w.filter_damaged_taxa(
        df=mdmg_results,
        filter_conditions=config["mdmg_filter_conditions"],
        taxonomic_rank=taxonomic_rank,
    )
    logging.info(f"Found {damaged_taxa.shape[0]} damaged taxa")

    logging.info("Getting proportion at rank ...")
    genome_abundance = w.aggregate_woltka_results_by_rank(
        df=woltka_profiling_results, rank=taxonomic_rank
    )
    # add column to genome_abundance with damaged status
    genome_abundance["is_damaged"] = (
        genome_abundance[[taxonomic_rank]]
        .isin(damaged_taxa["tax_id"].to_list())
        .astype(str)
    )
    genome_abundance["is_damaged"].replace(
        {"None": None, "True": True, "False": None, np.nan: None}, inplace=True
    )
    # Select genomes that are not damaged based on the config
    non_damaged_genomes = (
        select_taxa(
            genome_abundance,
            is_damaged=False,
            mode=config["max_genomes_nondamaged_selection"],
            n_taxa=config["max_genomes_nondamaged"],
        )
        .merge(woltka_profiling_results.drop("counts", axis=1), on=taxonomic_rank)
        .merge(woltka_stats, on="reference")
    )
    damaged_genomes = (
        select_taxa(
            genome_abundance,
            is_damaged=True,
            mode=config["max_genomes_damaged_selection"],
            n_taxa=config["max_genomes_damaged"],
        )
        .merge(woltka_profiling_results.drop("counts", axis=1), on=taxonomic_rank)
        .merge(woltka_stats, on="reference")
    )

    logging.info("Re-estimating proportions...")

    genomes = f.concat_df([damaged_genomes, non_damaged_genomes])

    genomes["proportion_new"] = 100 * (genomes["counts"] / genomes["counts"].sum())
    logging.info("Generating files for aMGSIM")
    # generate aMGSIM input files
    # community file
    community_file = generate_community_file(
        genomes,
        sample_name=config["sample_name"],
        taxonomic_rank=taxonomic_rank,
        taxdb=taxdb,
        use_restimated_proportions=config["use_restimated_proportions"],
    )
    genome_compositions = generate_genome_compositions(
        df=genomes,
        sample_name=config["sample_name"],
        taxdb=taxdb,
        taxonomic_rank=taxonomic_rank,
    )
    genome_paths_file = generate_genome_paths(
        df=genomes,
        taxdb=taxdb,
        genome_paths=genome_paths,
        taxonomic_rank=taxonomic_rank,
    )
    logging.info(f"Writing results")
    out_files = create_output_files(prefix=config["sample_name"])
    community_file.to_csv(out_files["communities"], sep="\t", index=False)
    genome_paths_file.to_csv(out_files["paths"], sep="\t", index=False)
    genome_compositions.to_csv(out_files["compositions"], sep="\t", index=False)


def opt_parse(args=None):
    version = "Version: " + __version__
    if args is None:
        args = docopt(__doc__, version=version)
    else:
        args = docopt(__doc__, version=version, argv=args)
    main(args)
