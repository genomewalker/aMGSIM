import logging
from tkinter import N
import pandas as pd
from itertools import takewhile
import numpy as np
import taxopy as txp
from multiprocessing import Pool
import tqdm
from functools import partial
from aMGSIM.library import functions as f
from aMGSIM.library import cli as c
import sys

log = logging.getLogger("my_logger")

# From https://stackoverflow.com/a/59737793/15704171
def largestRemainderMethod(pd_series, decimals=1):

    floor_series = ((10**decimals * pd_series).astype(np.int)).apply(np.floor)
    diff = 100 * (10**decimals) - floor_series.sum().astype(np.int)
    series_decimals = pd_series - floor_series / (10**decimals)
    series_sorted_by_decimals = series_decimals.sort_values(ascending=False)

    for i in range(0, len(series_sorted_by_decimals)):
        if i < diff:
            series_sorted_by_decimals.iloc[[i]] = 1
        else:
            series_sorted_by_decimals.iloc[[i]] = 0

    out_series = (
        (floor_series + series_sorted_by_decimals) / (10**decimals)
    ).sort_values(ascending=False)

    return out_series


def get_comments(filename):
    headiter = takewhile(lambda s: s.startswith("#"), filename)
    # you may want to process the headers differently,
    # but here we just convert it to a list
    header = list(headiter)
    return header


def load_genome_paths(file_path):
    """Function to read a genome path file to a pandas dataframe

    Args:
        file_path (str): A file path pointing to a genome path file

    Returns:
        pandas.DataFrame: A pandas dataframe containing the genome path
    """
    genome_path = pd.read_csv(
        file_path, sep="\t", header=None, index_col=None, names=["reference", "path"]
    )
    return genome_path


def load_genome_mapping_statistics(file_path):
    """Function to read a genome mapping statistics from bam-filter to a pandas dataframe

    Args:
        file_path (str): A file path pointing to a tsv file with genome mapping statistics

    Returns:
        pandas.DataFrame: A pandas dataframe containing the genome mapping statistics
    """
    genome_mapping_statistics = pd.read_csv(file_path, sep="\t", index_col=None)
    return genome_mapping_statistics


def load_mdmg_results(file_path):
    """Function to read a mdmg results file to a pandas dataframe

    Args:
        file_path (str): A file path pointing to a mdmg results file

    Returns:
        pandas.DataFrame: A pandas dataframe containing the mdmg results
    """
    mdmg_results = pd.read_csv(file_path, sep=",", index_col=None)
    mdmg_results.rename(columns={"tax_id": "reference"}, inplace=True)
    return mdmg_results


def load_misincorporation_results(file_path):
    """Function to read a misincorporation results file to a pandas dataframe

    Args:
        file_path (str): A file path pointing to a misincorporation results file
    Returns:
        pandas.DataFrame: A pandas dataframe containing the misincorporation results
    """
    misincorporation_comments = get_comments(file_path)
    misincorporation_results = pd.read_csv(file_path, index_col=None, comment="#")
    return misincorporation_comments, misincorporation_results


def load_filterBAM_stats_table(file_path):
    """Function to read a filterBAM table to a pandas dataframe

    Args:
        file_path (str): A file path pointing to a filterBAM table
    Returns:
        pandas.DataFrame: A pandas dataframe containing the filterBAM table
    """
    filterBAM_results = pd.read_csv(file_path, sep="\t", index_col=None)
    return filterBAM_results


def load_filterBAM_profiling_table(file_path):
    """Function to read a filterBAM table to a pandas dataframe

    Args:
        file_path (str): A file path pointing to a filterBAM table
    Returns:
        pandas.DataFrame: A pandas dataframe containing the filterBAM table
    """
    filterBAM_results = pd.read_csv(
        file_path,
        sep="\t",
        index_col=None,
        header=0,
        names=["reference", "counts", "name", "rank", "lineage"],
    )
    filterBAM_results[
        ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    ] = [list(map(int, i.split(";"))) for i in filterBAM_results["lineage"]]
    return filterBAM_results.filter(
        [
            "reference",
            "counts",
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]
    )


def aggregate_filterBAM_results_by_rank(df, rank):
    """Function to aggregate filterBAM results by rank

    Args:
        df (panda.DataFrame): DataFrame containing filterBAM results
        rank (str): Taxonomic scope

    Returns:
        panda.DataFrame containing: A dataframe where abundances are aggregated by rank
    """
    # select columns from dataframe by rank
    df_sub = df.filter(["reference", "tax_abund_read", "taxrank", "is_damaged"])
    df_agg = (
        df_sub.groupby(["taxrank", "is_damaged"])
        # .agg({"counts": "sum"})
        .agg({"taxrank": "count", "tax_abund_read": "sum"})
        .rename(columns={"taxrank": "n_genomes"})
        .reset_index(["taxrank", "is_damaged"])
    )

    # calculate percentage of total counts
    df_agg["proportion"] = 100 * (
        df_agg["tax_abund_read"] / df_agg["tax_abund_read"].sum()
    )
    # sort columns by proportion
    df_agg = df_agg.sort_values(by="proportion", ascending=False)
    df_agg["rank"] = np.arange(len(df_agg)) + 1
    return df_agg


def filter_damaged_taxa(df, filter_conditions, taxonomic_rank):
    """Function to filter damaged taxa

    Args:
        df (panda.DataFrame): A dataframe containing metaDMG results
        filter_conditions (dict): A dictionary with filter conditions
        taxonomic_rank (str): Select the taxonomic rank to filter

    Returns:
        pandas.DataFrame: A filtered dataframe containing metaDMG results
    """
    # filter rows by d_max, phi, qvalue
    # mdmg_results = df.loc[
    #     (df["D_max"] >= filter_conditions["d_max"])
    #     & (df["phi"] >= filter_conditions["phi"])
    #     & (df["q"] >= filter_conditions["q"])
    #     & (df["tax_rank"] == taxonomic_rank)
    # ]

    mdmg_results = df.loc[
        (df[list(filter_conditions)] >= pd.Series(filter_conditions)).all(axis=1)
    ].copy()

    return mdmg_results


def filter_filterBAM_taxa(df, filter_conditions):
    """Function to filter filterBAM taxa

    Args:
        df (pandas.DataFrame): A dataframe containing filterBAM results
        filter_conditions (dict): A dictionary with filter conditions

    Returns:
        pandas.DataFrame: A filtered filterBAM results dataframe
    """
    # filter rows by d_max, phi, qvalue
    # filterBAM_results = df.loc[
    #     (df["breadth_exp_ratio"] >= filter_conditions["breadth_exp_ratio"])
    #     & (df["coverage_mean"] >= filter_conditions["coverage_mean"])
    # ]
    filterBAM_results = df.loc[
        (df[list(filter_conditions)] >= pd.Series(filter_conditions)).all(axis=1)
    ]
    return filterBAM_results


def get_tax(ref, parms):
    taxdb = parms["taxdb"]
    acc2taxid = parms["acc2taxid"]
    if ref in acc2taxid:
        taxid = acc2taxid[ref]
        # taxid = txp.taxid_from_name(ref, taxdb)[0]
        taxonomy_info = txp.Taxon(taxid, taxdb).rank_name_dictionary
        taxonomy_info["taxid"] = taxid
        taxonomy_info["ref"] = ref
    else:
        log.debug(f"No taxid found for {ref}")
        taxonomy_info = None
    return taxonomy_info


def get_taxonomy_info(refids, taxdb, acc2taxid, nprocs=1):
    """Function to get the references taxonomic information for a given taxonomy id

    Args:
        taxids (list): A list of taxonomy ids
        taxdb (taxopy.TaxonomyDB): A taxopy DB

    Returns:
        dict: A list of taxonomy information
    """

    acc2taxid_df = pd.read_csv(acc2taxid, sep="\t", index_col=None)[
        ["accession", "taxid"]
    ].rename(columns={"accession": "reference"}, inplace=False)
    # Filter rows in refids from dataframe
    acc2taxid_df = acc2taxid_df.loc[acc2taxid_df["reference"].isin(refids)]
    acc2taxid_dict = acc2taxid_df.set_index("reference").T.to_dict("records")

    debug = c.is_debug()
    logging.getLogger("my_logger").setLevel(logging.DEBUG if debug else logging.INFO)

    parms = {"taxdb": taxdb, "acc2taxid": acc2taxid_dict[0]}
    func = partial(get_tax, parms=parms)
    if debug is True or len(refids) < 100000:
        taxonomy_info = list(map(func, refids))
    else:
        p = Pool(nprocs, initializer=f.initializer, initargs=(parms,))
        c_size = f.calc_chunksize(nprocs, len(refids))
        taxonomy_info = list(
            tqdm.tqdm(
                p.imap_unordered(func, refids, chunksize=c_size),
                total=len(refids),
                leave=False,
                ncols=100,
                desc=f"References processed",
            )
        )
        p.close()
        p.join()
    taxonomy_info = list(filter(None, taxonomy_info))
    exclude = ["taxid", "ref"]
    tax_ranks = []

    for k in taxonomy_info[0].keys():
        if k not in exclude:
            tax_ranks.append(k)

    taxonomy_info = {i["ref"]: i for i in taxonomy_info}
    return taxonomy_info, tax_ranks


def check_if_genomes_left(df):
    if df.shape[0] == 0:
        log.error(f"::: No damaged genomes found. Exiting.")
        sys.exit(1)


def check_if_nrows_is_smaller_than_ntaxa(df, n_taxa):
    if n_taxa > df.shape[0] and df.shape[0] > 0:
        n_taxa = len(df)
        log.info(f"::: The number of taxa is too high. Setting to {n_taxa}")
    return n_taxa


def get_missing_genomes(df, stats_df, n_taxa, n_taxa_first, tax_rank):
    log.info(
        f"::: We tried to get {n_taxa:,} genomes at {tax_rank} level but only found {n_taxa_first:,}."
    )
    log.info(f"::: Resampling the genome table.")
    # we will sample the table with genomes to get as closer as we can
    # while mainting a diverse set of genomes
    # get an initial set of genomes as big as we can
    stats_df = stats_df.loc[~stats_df["reference"].isin(df["reference"])]
    n_taxa_left = n_taxa - n_taxa_first
    if n_taxa_left <= stats_df.shape[0]:
        df2 = stats_df.sample(n=n_taxa_left)
    else:
        log.info(
            f"::: Unfortunately there are not enough genomes, selecting {stats_df.shape[0]:,} genomes left."
        )
        df2 = stats_df
    df = f.concat_df([df, df2])
    return df


def select_taxa(
    df, is_damaged=False, mode="most_abundant", n_taxa=10, tax_rank=None, stats=None
):
    """Function to select taxa to be simulated.
    We will select the target taxonomic rank, and choose by random or by most abundant or
    by least abundance. Then we will try to get one genome out of each taxonomic rank, if
    the selected number of genomes is not fulfilled we will look for another rank.

    Args:
        df (pandas.DataFrame): [description]
        is_damaged (bool, optional): [description]. Defaults to False.
        mode (str, optional): [description]. Defaults to 'most_abundant'.
        n_taxa (int, optional): [description]. Defaults to 100.
        stats (pandas.DataFrame): [description]. Defaults to None.

    Raises:
        ValueError: Fails if the mode is not valid.
    """
    n_taxa_first = 0
    if mode == "most_abundant":
        log.info(f"::: Selecting most abundant taxa")
        if is_damaged:
            df = df[df["is_damaged"] == True]
            check_if_genomes_left(df)
            stats_df = stats[stats["is_damaged"] == True]
        else:
            df = df[df["is_damaged"] != True]
            if df.shape[0] == 0:
                return df
            stats_df = stats[stats["is_damaged"] != True]
        # select n rows from pandas DataFrame
        n_taxa_first = check_if_nrows_is_smaller_than_ntaxa(df, n_taxa)
        if n_taxa_first == n_taxa:
            df = df.head(n_taxa)
            # Filter the genoems in stats that are found in df
            stats_df = stats_df.loc[stats_df["taxrank"].isin(df["taxrank"])]
            df = (
                stats_df.sort_values("tax_abund_read", ascending=False)
                .groupby("taxrank")
                .nth(0)
                .sort_values("tax_abund_read", ascending=False)
            )
        else:
            df1 = (
                stats_df.sort_values("tax_abund_read", ascending=False)
                .groupby("taxrank")
                .nth(0)
                .sort_values("tax_abund_read", ascending=False)
            )
            df = get_missing_genomes(
                df=df1,
                stats_df=stats_df,
                n_taxa=n_taxa,
                n_taxa_first=n_taxa_first,
                tax_rank=tax_rank,
            )

    elif mode == "least_abundant":
        log.info(f"::: Selecting least abundant taxa")
        if is_damaged:
            df = df[df["is_damaged"] == True]
            check_if_genomes_left(df)
            stats_df = stats[stats["is_damaged"] == True]
        else:
            df = df[df["is_damaged"] != True]
            if df.shape[0] == 0:
                return df
            stats_df = stats[stats["is_damaged"] != True]
        # select n rows from pandas DataFrame
        n_taxa_first = check_if_nrows_is_smaller_than_ntaxa(df, n_taxa)
        if n_taxa_first == n_taxa:
            df = df.head(n_taxa)
            # Filter the genoems in stats that are found in df
            stats_df = stats_df.loc[stats_df["taxrank"].isin(df["taxrank"])]
            df = (
                stats_df.sort_values("tax_abund_read", ascending=True)
                .groupby("taxrank")
                .nth(0)
                .sort_values("tax_abund_read", ascending=True)
            )
        else:
            df1 = (
                stats_df.sort_values("tax_abund_read", ascending=True)
                .groupby("taxrank")
                .nth(0)
                .sort_values("tax_abund_read", ascending=True)
            )
            df = get_missing_genomes(
                df=df1,
                stats_df=stats_df,
                n_taxa=n_taxa,
                n_taxa_first=n_taxa_first,
                tax_rank=tax_rank,
            )
    elif mode == "random":
        log.info("::: Selecting random taxa")
        if is_damaged:
            df = df[df["is_damaged"] == True]
            check_if_genomes_left(df)
            n_taxa_first = check_if_nrows_is_smaller_than_ntaxa(df, n_taxa)
            df = df.sample(n=n_taxa_first)
            stats_df = stats[stats["is_damaged"] == True]
        else:
            df = df[df["is_damaged"] != True]
            if df.shape[0] == 0:
                return df
            n_taxa_first = check_if_nrows_is_smaller_than_ntaxa(df, n_taxa)
            df = df.sample(n=n_taxa_first)
            stats_df = stats[stats["is_damaged"] != True]
        stats_df = stats_df.loc[stats_df["taxrank"].isin(df["taxrank"])]

        if n_taxa_first == n_taxa:
            df = stats_df.groupby("taxrank").sample(n=1)
        else:
            df1 = (
                stats_df.groupby("taxrank")
                .sample(n=1)
                .sort_values("tax_abund_read", ascending=False)
            )
            df = get_missing_genomes(
                df=df1,
                stats_df=stats_df,
                n_taxa=n_taxa,
                n_taxa_first=n_taxa_first,
                tax_rank=tax_rank,
            )
    else:
        raise ValueError(
            'mode should be one of "most_abundant", "random", "least_abundant"'
        )
    return df
