#!/usr/bin/env python
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
from aMGSIM.library import filterBAM_functions as w
import datetime
import tqdm
from collections import OrderedDict
from aMGSIM import __version__
import taxopy as txp
from aMGSIM.library import cli as c


def exceptionHandler(
    exception_type, exception, traceback, debug_hook=sys.__excepthook__
):
    """Print user friendly error messages normally, full traceback if DEBUG on.
    Adapted from http://stackoverflow.com/questions/27674602/hide-traceback-unless-a-debug-flag-is-set
    """
    debug = c.is_debug()
    if debug:
        print("\n*** Error:")
        # raise
        debug_hook(exception_type, exception, traceback)
    else:
        print("{}: {}".format(exception_type.__name__, exception))
        log.error("Please use --debug to see full traceback.")


def generate_genome_compositions(df, sample_name):

    """Function to generate a genome composition file compatible with aMGSIM.

    Args:
        df (pandas.DataFrame): A dataframe with the genome information.
        sample_name (str): Sample name

    Returns:
        pandas.DataFrame: A table with the genome information.
    """
    df_selected = df[
        [
            "Taxon",
            "coverage_mean",
            "read_length_mode",
            "read_length_std",
            "read_length_min",
            "read_length_max",
            "is_damaged",
        ]
    ].copy()
    df_selected.loc[:, "Community"] = sample_name
    df_selected = df_selected[
        [
            "Taxon",
            "Community",
            "coverage_mean",
            "read_length_mode",
            "read_length_std",
            "read_length_min",
            "read_length_max",
            "is_damaged",
        ]
    ]
    df_selected.columns = [
        "Taxon",
        "Community",
        "Coverage",
        "Read_length",
        "Read_length_std",
        "Read_length_min",
        "Read_length_max",
        "onlyAncient",
    ]
    # convert onlyAncient column to string
    df_selected.loc[:, "onlyAncient"] = df_selected["onlyAncient"].astype(str)

    return df_selected


def generate_community_file(df, sample_name, use_restimated_proportions):
    """Function to create a community file compatible with aMGSIM.

    Args:
        df (pandas.DataFrame): A dataframe with the genome information.
    """
    # create a new DataFrame with the selected taxa
    if use_restimated_proportions:
        df_selected = df[
            [
                "reference",
                "Taxon",
                "proportion_new",
                "read_length_mode",
                "read_length_min",
                "read_length_max",
                "read_length_std",
            ]
        ].copy()
        df_selected.loc[:, "Community"] = sample_name
        df_selected["Rank"] = np.arange(len(df_selected)) + 1
        df_selected = df_selected[
            [
                "Community",
                "Taxon",
                "Rank",
                "read_length_mode",
                "read_length_std",
                "read_length_min",
                "read_length_max",
                "proportion_new",
            ]
        ]
        df_selected.columns = [
            "Community",
            "Taxon",
            "Rank",
            "Read_length",
            "Read_length_std",
            "Read_length_min",
            "Read_length_max",
            "Perc_rel_abund",
        ]
    else:
        df_selected = df[
            [
                "reference",
                "Taxon",
                "proportion",
                "read_length_mode",
                "read_length_std",
                "read_length_min",
                "read_length_max",
            ]
        ].copy()
        df_selected.loc[:, "Community"] = sample_name
        df_selected["Rank"] = np.arange(len(df_selected)) + 1
        df_selected = df_selected[
            [
                "Community",
                "Taxon",
                "Rank",
                "read_length_mode",
                "read_length_std",
                "read_length_min",
                "read_length_max",
                "proportion",
            ]
        ]
        df_selected.columns = [
            "Community",
            "Taxon",
            "Rank",
            "Read_length",
            "Read_length_std",
            "Read_length_min",
            "Read_length_max",
            "Perc_rel_abund",
        ]

    return df_selected


def generate_genome_paths(df, genome_paths):
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
        df[["reference", "taxid", "Taxon"]].copy().merge(genome_paths, on="reference")
    )
    df_selected = df_selected[["Taxon", "taxid", "reference", "path"]]
    df_selected.columns = ["Taxon", "TaxId", "Accession", "Fasta"]
    return df_selected


def create_output_files(prefix):
    # create output files
    out_files = {
        "paths": f"{prefix}.filepaths.tsv",
        "compositions": f"{prefix}.genome-compositions.tsv",
        "communities": f"{prefix}.communities.tsv",
    }
    return out_files


log = logging.getLogger("my_logger")


def estimate(args):
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    # simulating reads
    if args.debug:
        debug = True
    else:
        debug = False

    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )

    sys.excepthook = exceptionHandler
    # simulating reads
    cfg = vars(args)
    cfg = f.validate_schema(cfg, d.schema_init_w, debug)
    config = f.get_config(config=cfg["config"], schema=d.w_schema_config, debug=debug)
    taxonomic_rank = config["taxonomic-rank"]
    # load tables
    log.info("Loading genome paths file...")
    genome_paths = w.load_genome_paths(config["genome-paths"])
    log.info("Loading filterBAM stats file...")
    filterBAM_stats = w.load_filterBAM_stats_table(config["filterBAM-stats"])
    filterBAM_stats_n = filterBAM_stats.shape[0]
    log.info("Generating taxonomic information...")
    taxdb = txp.TaxDb(
        nodes_dmp=config["nodes-dmp"], names_dmp=config["names-dmp"], keep_files=True
    )

    taxonomy_info, tax_ranks = w.get_taxonomy_info(
        refids=filterBAM_stats["reference"].values,
        taxdb=taxdb,
        acc2taxid=config["acc2taxid"],
        nprocs=config["cpus"],
    )
    # get list of keys for taxonomy
    taxonomy_keys = list(taxonomy_info.keys())

    # filter table to contain taxonomy keys as reference
    filterBAM_stats = filterBAM_stats[filterBAM_stats["reference"].isin(taxonomy_keys)]

    filterBAM_stats["taxid"] = filterBAM_stats["reference"].apply(
        lambda x: taxonomy_info[x]["taxid"]
    )
    filterBAM_stats["taxrank"] = filterBAM_stats["reference"].apply(
        lambda x: taxonomy_info[x][taxonomic_rank]
    )

    filterBAM_stats = w.filter_filterBAM_taxa(
        df=filterBAM_stats,
        filter_conditions=config["filterBAM-filter-conditions"],
    )
    log.info(f"Kept {filterBAM_stats.shape[0]} taxa from {filterBAM_stats_n}")

    log.info("Loading metaDMG results...")
    mdmg_results = w.load_mdmg_results(config["mdmg-results"])
    # find which taxon are damaged

    damaged_taxa = w.filter_damaged_taxa(
        df=mdmg_results,
        filter_conditions=config["mdmg-filter-conditions"],
        taxonomic_rank=taxonomic_rank,
    )
    # damaged_taxa["taxrank"] = damaged_taxa["reference"].apply(
    #     lambda x: taxonomy_info[x][taxonomic_rank]
    # )

    # add column to genome_abundance with damaged status
    filterBAM_stats["is_damaged"] = (
        filterBAM_stats["reference"]
        .isin(damaged_taxa["reference"].to_list())
        .astype(str)
    )

    filterBAM_stats["is_damaged"].replace(
        {"None": False, "True": True, np.nan: False, "False": False}, inplace=True
    )
    filterBAM_stats["proportion"] = 100 * (
        filterBAM_stats["tax_abund_read"] / filterBAM_stats["tax_abund_read"].sum()
    )

    log.info(f"Found {damaged_taxa.shape[0]} damaged taxa")

    log.info("Getting proportion at rank ...")
    rank_abundance = w.aggregate_filterBAM_results_by_rank(
        df=filterBAM_stats, rank=taxonomic_rank
    )

    # Select genomes that are not damaged based on the config
    if config["max-genomes-nondamaged"] > 0:
        log.info("Selecting genomes that are not damaged")
        non_damaged_genomes = w.select_taxa(
            rank_abundance,
            is_damaged=False,
            mode=config["max-genomes-nondamaged-selection"],
            n_taxa=config["max-genomes-nondamaged"],
            tax_rank=taxonomic_rank,
            stats=filterBAM_stats,
        )
        non_damaged_genomes["is_damaged"] = None
    else:
        log.info("Skipping selection of genomes that are not damaged")
        non_damaged_genomes = pd.DataFrame()

    log.info("Selecting genomes that are damaged")
    damaged_genomes = w.select_taxa(
        rank_abundance,
        is_damaged=True,
        mode=config["max-genomes-damaged-selection"],
        n_taxa=config["max-genomes-damaged"],
        tax_rank=taxonomic_rank,
        stats=filterBAM_stats,
    )
    # Check that both tables have entries, if not exit
    if (non_damaged_genomes.shape[0] == 0) and (damaged_genomes.shape[0] == 0):
        log.error("No genomes found")
        sys.exit(1)

    # merge tables
    log.info("Combining tables and re-estimating proportions...")
    if non_damaged_genomes.shape[0] > 0:
        genomes = f.concat_df([damaged_genomes, non_damaged_genomes])
    else:
        genomes = damaged_genomes

    genomes["proportion_new"] = 100 * (
        genomes["tax_abund_read"] / genomes["tax_abund_read"].sum()
    )

    genomes["Taxon"] = genomes["reference"].apply(
        lambda x: re.sub(
            "[^0-9a-zA-Z]+", "_", re.sub("s__", "", taxonomy_info[x]["species"])
        )
    )
    genomes["Taxon"] = genomes[["Taxon", "reference"]].apply(
        lambda row: "----".join(row.values.astype(str)), axis=1
    )
    log.info("Generating files for aMGSIM")
    # generate aMGSIM input files
    # community file
    community_file = generate_community_file(
        genomes,
        sample_name=config["sample-name"],
        use_restimated_proportions=config["use-restimated-proportions"],
    )

    genome_compositions = generate_genome_compositions(
        df=genomes,
        sample_name=config["sample-name"],
    )
    genome_compositions["reference"] = genome_compositions["Taxon"].str.split(
        "----", expand=True
    )[1]

    if "Bayesian_D_max" in mdmg_results.columns:
        genome_compositions = genome_compositions.merge(
            mdmg_results[["reference", "Bayesian_D_max"]], on="reference"
        )
        # rename column Bayesian_D_max to D_max
        genome_compositions.rename(columns={"Bayesian_D_max": "D_max"}, inplace=True)
        # drop reference column
    else:
        genome_compositions = genome_compositions.merge(
            mdmg_results[["reference", "D_max"]], on="reference"
        )

    genome_compositions.drop(columns=["reference"], inplace=True)

    genome_compositions.Read_length = genome_compositions.Read_length.round(0).astype(
        int
    )
    genome_paths_file = generate_genome_paths(
        df=genomes,
        genome_paths=genome_paths,
    )

    log.info(f"Writing results")
    out_files = create_output_files(prefix=config["sample-name"])
    community_file.to_csv(out_files["communities"], sep="\t", index=False)
    genome_paths_file.to_csv(out_files["paths"], sep="\t", index=False)
    genome_compositions.to_csv(out_files["compositions"], sep="\t", index=False)
