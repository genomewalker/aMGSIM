import pandas as pd
from itertools import takewhile
import numpy as np

# From https://stackoverflow.com/a/59737793/15704171
def largestRemainderMethod(pd_series, decimals=1):

    floor_series = ((10 ** decimals * pd_series).astype(np.int)).apply(np.floor)
    diff = 100 * (10 ** decimals) - floor_series.sum().astype(np.int)
    series_decimals = pd_series - floor_series / (10 ** decimals)
    series_sorted_by_decimals = series_decimals.sort_values(ascending=False)

    for i in range(0, len(series_sorted_by_decimals)):
        if i < diff:
            series_sorted_by_decimals.iloc[[i]] = 1
        else:
            series_sorted_by_decimals.iloc[[i]] = 0

    out_series = (
        (floor_series + series_sorted_by_decimals) / (10 ** decimals)
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


def load_woltka_stats_table(file_path):
    """Function to read a woltka table to a pandas dataframe

    Args:
        file_path (str): A file path pointing to a woltka table
    Returns:
        pandas.DataFrame: A pandas dataframe containing the woltka table
    """
    woltka_results = pd.read_csv(file_path, sep="\t", index_col=None)
    return woltka_results


def load_woltka_profiling_table(file_path):
    """Function to read a woltka table to a pandas dataframe

    Args:
        file_path (str): A file path pointing to a woltka table
    Returns:
        pandas.DataFrame: A pandas dataframe containing the woltka table
    """
    woltka_results = pd.read_csv(
        file_path,
        sep="\t",
        index_col=None,
        header=0,
        names=["reference", "counts", "name", "rank", "lineage"],
    )
    woltka_results[
        ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    ] = [list(map(int, i.split(";"))) for i in woltka_results["lineage"]]
    return woltka_results.filter(
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


def aggregate_woltka_results_by_rank(df, rank):
    """Function to aggregate woltka results by rank

    Args:
        df (panda.DataFrame): DataFrame containing woltka results
        rank (str): Taxonomic scope

    Returns:
        panda.DataFrame containing: A dataframe where abundances are aggregated by rank
    """
    # select columns from dataframe by rank
    df_sub = df.filter(["reference", "counts", rank])
    df_agg = (
        df_sub.groupby(rank)
        # .agg({"counts": "sum"})
        .agg({rank: "count", "counts": "sum"})
        .rename(columns={rank: "n_genomes"})
        .reset_index()
    )
    # calculate percentage of total counts
    df_agg["proportion"] = 100 * (df_agg["counts"] / df_agg["counts"].sum())
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
    mdmg_results = df.loc[
        (df["D_max"] >= filter_conditions["d_max"])
        & (df["phi"] >= filter_conditions["phi"])
        & (df["q"] >= filter_conditions["q"])
        & (df["tax_rank"] == taxonomic_rank)
    ]

    return mdmg_results


def filter_woltka_taxa(df, filter_conditions):
    """Function to filter woltka taxa

    Args:
        df (pandas.DataFrame): A dataframe containing woltka results
        filter_conditions (dict): A dictionary with filter conditions

    Returns:
        pandas.DataFrame: A filtered woltka results dataframe
    """
    # filter rows by d_max, phi, qvalue
    woltka_results = df.loc[
        (df["breadth_exp_ratio"] >= filter_conditions["breadth_exp_ratio"])
        & (df["coverage_mean"] >= filter_conditions["coverage_mean"])
    ]
    return woltka_results