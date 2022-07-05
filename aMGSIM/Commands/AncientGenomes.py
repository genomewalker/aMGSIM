#!/usr/bin/env python

from distutils.log import info
from statistics import mode
from docopt import docopt
import logging
import pandas as pd
from functools import partial
from multiprocessing import Pool
import numpy as np
import json
from pathlib import Path
import sys
import math
from scipy.stats import lognorm

# application
from MGSIM import SimReads
from aMGSIM.library import defaults as d
from aMGSIM.library import functions as f
import datetime
import tqdm
from collections import OrderedDict
from aMGSIM import __version__

debug = False


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
        print(f"{exception_type.__name__}: {exception}")
        print("\n*** Error: Please use --debug to see full traceback.")


logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s ::: %(asctime)s ::: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def obj_dict(obj):
    return obj.__dict__


def process_genome(
    x,
    dist,
    mode_len_ancient,
    mode_len_modern,
    min_len,
    max_len,
    coverage,
    library,
    read_length_freqs,
    selected_genomes,
    genome_comp,
    n_reads,
):
    genome = Genome(
        genome=x,
        dist=dist,
        mode_len_ancient=mode_len_ancient,
        mode_len_modern=mode_len_modern,
        min_len=min_len,
        max_len=max_len,
        coverage=coverage,
        library=library,
        read_length_freqs=read_length_freqs,
        # seqSys=seqSys,
        selected_genomes=selected_genomes,
        genome_comp=genome_comp,
        n_reads=n_reads,
    )
    return genome


class Genome:
    def __init__(
        self,
        genome,
        dist,
        mode_len_ancient,
        mode_len_modern,
        selected_genomes,
        genome_comp,
        min_len,
        max_len,
        coverage,
        library,
        read_length_freqs,
        n_reads,
        *args,
        **kwargs,
    ):
        # Community, Taxon, Perc_rel_abund, Rank, Fasta, Genome_size
        rnd_seed = int(np.random.RandomState().randint(0, 99999, size=1))
        np.random.seed(rnd_seed)
        # ["Community", "Taxon", "Perc_rel_abund", "Rank", "Genome_size"]
        self.comm = genome["Community"]
        self.taxon = genome["Taxon"]
        if "Perc_rel_abund" in genome:
            self.rel_abund = genome["Perc_rel_abund"]
        else:
            self.rel_abund = None
        self.genome_size = genome["Genome_size"]
        # test if key is in dictionary

        if "Accession" in genome.keys():
            self.accession = genome["Accession"]
        else:
            self.accession = None

        # if "Read_length" in genome_comp.keys():
        #     rl = genome_comp[genome_comp["Taxon"] == self.taxon].iloc[0]["Read_length"]
        #     mode_len_ancient = rl
        #     mode_len_modern = rl

        # if "Read_length_std" in genome_comp.keys():
        #     stddev = genome_comp[genome_comp["Taxon"] == self.taxon].iloc[0][
        #         "Read_length_std"
        #     ]
        #     stddev_len_ancient = stddev
        #     stddev_len_modern = stddev
        # else:
        #     stddev_len_ancient = mode_len_ancient - min_len
        #     stddev_len_modern = mode_len_modern - min_len

        if "Read_length_min" in genome_comp.keys():
            min_len = genome_comp[genome_comp["Taxon"] == self.taxon].iloc[0][
                "Read_length_min"
            ]
        if "Read_length_max" in genome_comp.keys():
            max_len = genome_comp[genome_comp["Taxon"] == self.taxon].iloc[0][
                "Read_length_max"
            ]
        # self.library = library
        # self.read_length = max_len
        # self.seqSys = seqSys
        # self.n_reads = n_reads

        # Get the genome from the selected genomes table
        selected_genomes_filt = selected_genomes[
            (selected_genomes["Community"] == self.comm)
            & (selected_genomes["Taxon"] == self.taxon)
        ].copy()

        # selected_genomes_filt["onlyAncient"].replace(
        #     {"None": None, "True": True, "False": False}, inplace=True
        # )
        if genome_comp.empty:
            taxons = pd.DataFrame()
        else:
            taxons = genome_comp[
                (genome_comp["Taxon"] == self.taxon)
                & (genome_comp["Community"] == self.comm)
            ]

        # Is onlyAncient, mix or only moders
        # if len(selected_genomes_filt.index) > 0:
        #     self.onlyAncient = selected_genomes_filt.iloc[0]["onlyAncient"]
        # else:
        #     self.onlyAncient = None
        # If it is not modern get the coverage
        self.onlyAncient = selected_genomes_filt.iloc[0]["onlyAncient"]
        # if taxons.empty:
        #     if self.onlyAncient is not None:
        #         coverage_ancient = selected_genomes_filt.iloc[0]["coverage_ancient"]
        # else:
        #     coverage_ancient = selected_genomes_filt.iloc[0]["coverage_ancient"]
        # Get fragments
        # Get average size
        coverage_ancient = selected_genomes_filt.iloc[0]["coverage_ancient"]
        # if str(self.onlyAncient) != "None":
        #     self.onlyAncient = bool(self.onlyAncient)
        # else:
        #     self.onlyAncient = None
        if self.onlyAncient == True:
            self.onlyAncient = True
            self.fragments_ancient = self.random_fragments(
                min_len=min_len,
                max_len=max_len,
                dist=dist,
                read_length_freqs=read_length_freqs,
                mode_len=mode_len_ancient,
                rnd_seed=rnd_seed,
            )
            self.fragments_modern = None

            if library == "se":
                paired = 1
            else:
                paired = 2
            # logging.debug(f"Taxon:{self.taxon}")
            read_length_ancient = int(self.fragments_ancient["avg_len"] * paired)
            # logging.debug(f"read_length_ancient:{read_length_ancient}")
            # Number max of reads genome in sample
            seq_depth = int(n_reads * ((self.rel_abund / 100)))
            # logging.debug(f"seq_depth:{seq_depth}")
            # Seq depth based based on the size
            seq_depth_ancient = int(
                (coverage_ancient * self.genome_size) / read_length_ancient
            )
            # logging.debug(f"seq_depth_ancient:{seq_depth_ancient}")
            if taxons.empty:
                seq_depth_ancient_original = seq_depth_ancient
                seq_depth_ancient = seq_depth
            else:
                seq_depth_ancient_original = seq_depth_ancient
                seq_depth = seq_depth_ancient

            # logging.debug(f"seq_depth_ancient_original:{seq_depth_ancient_original}")
            # logging.debug(f"seq_depth:{seq_depth}")

            fold_ancient = (seq_depth_ancient * read_length_ancient) / self.genome_size
            # logging.debug(f"fold_ancient:{fold_ancient}")
            # logging.debug(f"coverage_ancient:{coverage_ancient}")

            if seq_depth_ancient > 0:
                self.fragments_ancient["seq_depth"] = seq_depth_ancient
                self.fragments_ancient[
                    "seq_depth_original"
                ] = seq_depth_ancient_original
                self.fragments_ancient["fold"] = fold_ancient
                self.fragments_ancient["fold_original"] = coverage_ancient
                if taxons.empty:
                    self.coverage_enforced = False
                else:
                    self.coverage_enforced = True
            else:
                self.fragments_ancient = None

            self.seq_depth = seq_depth
        elif self.onlyAncient is None:
            self.fragments_ancient = None
            self.fragments_modern = self.random_fragments(
                min_len=min_len,
                max_len=max_len,
                dist=dist,
                mode_len=mode_len_modern,
                read_length_freqs=read_length_freqs,
                rnd_seed=rnd_seed,
            )
            if library == "se":
                paired = 1
            else:
                paired = 2
            read_length_modern = int(self.fragments_modern["avg_len"] * paired)
            # Number max of reads genome in sample
            if taxons.empty:
                seq_depth_modern = int(n_reads * ((self.rel_abund / 100)))
                fold_modern = (seq_depth_modern * read_length_modern) / self.genome_size
            else:
                """
                One can define to have a modern organism with a certain coverage, if the defined
                coverage is larger that the maximum allowed, it gets downscaled to the max
                allowed. In case it is not, it will use the defined coverage, but the defined
                relative abundances will not be correct anymore. Also the final number of reads
                might be different than the expected.
                """

                # Calculate max allowed based on the number of reads
                seq_depth = int(n_reads * ((self.rel_abund / 100)))

                # Calculate defined coverage
                coverage_modern = taxons.iloc[0]["coverage_ancient"]
                seq_depth_modern = int(
                    (coverage_modern * self.genome_size) / read_length_modern
                )

                # if defined seq depth is larger than the max allow, downscale it
                if seq_depth_modern > seq_depth:
                    seq_depth_modern = seq_depth
                fold_modern = (seq_depth_modern * read_length_modern) / self.genome_size

            if seq_depth_modern > 0:
                self.fragments_modern["seq_depth"] = seq_depth_modern
                self.fragments_modern["fold"] = fold_modern
                self.fragments_modern["fold_original"] = coverage_modern

            else:
                self.fragments_modern = None
            if taxons.empty:
                self.coverage_enforced = False
            else:
                self.coverage_enforced = True
            self.seq_depth = seq_depth_modern

        else:
            self.onlyAncient = False
            self.fragments_ancient = self.random_fragments(
                min_len=min_len,
                max_len=max_len,
                dist=dist,
                read_length_freqs=read_length_freqs,
                mode_len=mode_len_ancient,
                rnd_seed=rnd_seed,
            )
            self.fragments_modern = self.random_fragments(
                min_len=min_len,
                max_len=max_len,
                dist=dist,
                read_length_freqs=read_length_freqs,
                mode_len=mode_len_modern,
                rnd_seed=rnd_seed,
            )
            if library == "se":
                paired = 1
            else:
                paired = 2
            read_length_ancient = int(self.fragments_ancient["avg_len"] * paired)
            read_length_modern = int(self.fragments_ancient["avg_len"] * paired)

            # Number max of reads genome in sample
            seq_depth = int(n_reads * ((self.rel_abund / 100)))

            seq_depth_ancient = int(
                (coverage_ancient * self.genome_size) / read_length_ancient
            )

            if seq_depth_ancient < seq_depth:
                seq_depth_modern = int(seq_depth - seq_depth_ancient)
            else:
                seq_depth_ancient = seq_depth
                seq_depth_modern = 0

            fold_modern = (seq_depth_modern * read_length_modern) / self.genome_size
            fold_ancient = (seq_depth_ancient * read_length_ancient) / self.genome_size
            if seq_depth_ancient > 0:
                self.fragments_ancient["seq_depth"] = seq_depth_ancient
                self.fragments_ancient["fold"] = fold_ancient
                self.fragments_ancient["fold_original"] = coverage_ancient
            else:
                self.fragments_ancient = None
            if seq_depth_modern > 0:
                self.fragments_modern["seq_depth"] = seq_depth_modern
                self.fragments_modern["fold"] = fold_modern
            else:
                self.fragments_modern = None

            if taxons.empty:
                self.coverage_enforced = False
            else:
                self.coverage_enforced = True
            self.seq_depth = seq_depth

    def random_fragments(
        self, min_len, max_len, dist, read_length_freqs, mode_len, rnd_seed
    ):
        if read_length_freqs:
            lengths = read_length_freqs[self.accession]["length"]
            freqs = read_length_freqs[self.accession]["freq"]
            lens = np.multiply(lengths, freqs).sum()
            avg_len = float(lens / sum(freqs))
            frags = {
                "fragments": {
                    "length": lengths,
                    "freq": freqs,
                },
                "dist_params": {
                    "scale": None,
                    "sigma": None,
                    "rnd_seed": None,
                },
                "avg_len": avg_len,
            }
        else:
            n_frags = int(1e6)
            stddev = mode_len - min_len
            sigma, scale = f.lognorm_params(mode=mode_len, stddev=stddev)
            mu = np.log(scale)
            dist_params = {"mean": mu, "sigma": sigma}
            # drawing relative abundances from the user-defined distribution
            freq_dist = f.get_freq_dist(dist, dist_params)
            dist_freqs = freq_dist(size=n_frags)
            dist_freqs = [int(x) for x in dist_freqs]
            # freqs = np.sort(dist_freqs / sum(dist_freqs) * 100)[::-1]
            # making a series for the taxa
            lengths = pd.Series(dist_freqs)
            lengths = lengths[(lengths >= min_len) & (lengths <= max_len)]
            lengths = lengths.value_counts().sort_index()
            freqs = list(lengths / sum(lengths))
            lens = np.multiply(list(lengths.index), freqs).sum()
            avg_len = float(lens / sum(freqs))
            # print("Lengths: {}", list(lengths.index))
            # print("Freqs: {}", list(lengths/sum(lengths)))
            # print("Sizes: {}".format(sizes))
            # print("Avg size: {}".format(avg_size))
            # exit()
            frags = {
                "fragments": {"length": list(lengths.index), "freq": freqs},
                "dist_params": {
                    "scale": float(scale),
                    "sigma": float(sigma),
                    "rnd_seed": rnd_seed,
                },
                "avg_len": avg_len,
            }
        return frags

    def __str__(self):
        return (
            str(self.__class__)
            + "\n"
            + "\n".join(
                (
                    str(item) + " = " + str(self.__dict__[item])
                    for item in sorted(self.__dict__)
                )
            )
        )

    def __repr__(self):
        return str(self)

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__)


def select_genomes(self, parms):
    prop = parms["prop"]
    coverage = parms["coverage"]
    n_reads = parms["n_reads"]
    library = parms["library"]
    mode_len_ancient = parms["mode_len_ancient"]
    cov_max = coverage["max"]
    read_len = parms["read_len"]
    rank_cutoff = parms["rank_cutoff"]
    comm = self["Community"]
    genome_comp = parms["genome_comp"]

    # we identify which genomes should have a specific ancient content in each commmunity
    genome_comp = genome_comp.merge(self)[["Taxon", "coverage_ancient", "onlyAncient"]]

    # check if we have Perc_rel_abund in self
    if "Perc_rel_abund" in self:
        rel_abund = True
    else:
        rel_abund = False
        self["Perc_rel_abund"] = None
        cov_anc = self["coverage_ancient"]

    if library == "se":
        paired = 1
    else:
        paired = 2

    # Estimate the ancient read length
    read_length_ancient = mode_len_ancient * paired
    # Get max and min coverage values
    min_cov = coverage["min"]
    max_cov = coverage["max"]

    # How many genomes we want to use
    # If not defined it will randomly pick some based on how many of them we want
    # to have in each community
    # Genome comp file defines spike-in taxa to be ancient with specific coverage
    k = int(round(prop * self.shape[0]))
    k_filt = 0

    if k > 0:
        while k_filt < k:
            rnd_seed = int(np.random.RandomState().randint(0, 99999, size=1))
            np.random.seed(rnd_seed)

            coverage_ancient = np.round(
                np.random.uniform(min_cov, max_cov, size=self.shape[0]), 1
            )
            if rel_abund:
                self["seq_depth"] = n_reads * (self["Perc_rel_abund"] / 100)
                self["seq_depth"] = self["seq_depth"].astype(int)
                self["max_ancient_cov_allowed"] = (
                    self["seq_depth"] * read_length_ancient
                ) / self["Genome_size"]
            else:
                self["max_ancient_cov_allowed"] = cov_anc
                self["seq_depth"] = (
                    cov_anc * self["Genome_size"]
                ) / read_length_ancient
                self["seq_depth"] = self["seq_depth"].astype(int)

            self["coverage_ancient_rnd"] = coverage_ancient.tolist()
            self["coverage_ancient"] = 0
            self["onlyAncient"] = self.apply(rand_isAncient, axis=1)
            self["onlyAncient"] = self["onlyAncient"].astype(str)

            if genome_comp.empty:
                # self['Rank_perd'] = self['Rank']/self['Rank'].max()
                self.loc[
                    (self["Rank"] / self["Rank"].max()) <= rank_cutoff, "onlyAncient"
                ] = bool(False)

            self.set_index("Taxon", inplace=True)
            self.update(genome_comp.set_index("Taxon"))
            self.reset_index(level=0, inplace=True)

            withAncient = self[(self["onlyAncient"].notnull())]
            k_filt = withAncient.shape[0]
    else:
        rnd_seed = int(np.random.RandomState().randint(0, 99999, size=1))
        np.random.seed(rnd_seed)
        coverage_ancient = np.round(
            np.random.uniform(min_cov, max_cov, size=self.shape[0]), 1
        )
        if rel_abund:
            self["seq_depth"] = n_reads * (self["Perc_rel_abund"] / 100)
            self["seq_depth"] = self["seq_depth"].astype(int)
            self["max_ancient_cov_allowed"] = (
                self["seq_depth"] * read_length_ancient
            ) / self["Genome_size"]
        else:
            self["max_ancient_cov_allowed"] = cov_anc
            self["seq_depth"] = (cov_anc * self["Genome_size"]) / read_length_ancient
            self["seq_depth"] = self["seq_depth"].astype(int)

        self["coverage_ancient_rnd"] = coverage_ancient.tolist()
        self["coverage_ancient"] = 0
        self["onlyAncient"] = self.apply(rand_isAncient, axis=1)

        # self['Rank_perd'] = self['Rank']/self['Rank'].max()
        if genome_comp.empty:
            self.loc[
                (self["Rank"] / self["Rank"].max()) <= rank_cutoff, "onlyAncient"
            ] = bool(False)
        self.set_index("Taxon", inplace=True)
        self.update(genome_comp.set_index("Taxon"))
        self.reset_index(level=0, inplace=True)
        withAncient = self[(self["onlyAncient"].notnull())]
        k_filt = withAncient.shape[0]

    pop = withAncient["Taxon"].tolist()
    # Estimate max number of reads as function of coverage and rel abun

    # Number max of reads genome in sample
    if rel_abund:
        w = pd.DataFrame(
            withAncient["Perc_rel_abund"].apply(lambda x: np.arcsinh(100 - x))
        )
        max_w = w["Perc_rel_abund"].sum()
        w = w["Perc_rel_abund"].apply(lambda x: x / max_w)

    if genome_comp.empty:
        random_genomes = np.random.choice(pop, replace=False, p=w, size=k)
        self.loc[~self["Taxon"].isin(random_genomes), "onlyAncient"] = None
        self["onlyAncient"].replace(
            {"None": None, "True": True, "False": False, np.nan: None}, inplace=True
        )
        self.loc[
            self["onlyAncient"].isnull(),
            "coverage_ancient",
        ] = 0
        self.loc[
            (self["onlyAncient"] == True)
            & (self["max_ancient_cov_allowed"] > self["coverage_ancient_rnd"]),
            "coverage_ancient",
        ] = self["coverage_ancient_rnd"]
        self.loc[
            (self["onlyAncient"] == True)
            & (self["max_ancient_cov_allowed"] <= self["coverage_ancient_rnd"]),
            "coverage_ancient",
        ] = self["max_ancient_cov_allowed"]
        self.loc[
            (self["onlyAncient"] == False),
            "coverage_ancient",
        ] = self["coverage_ancient_rnd"]
    else:
        random_genomes = genome_comp["Taxon"].to_list()
        self["onlyAncient"].replace(
            {"None": None, "True": True, "False": False, np.nan: None}, inplace=True
        )
        mapping = dict(genome_comp[["Taxon", "coverage_ancient"]].values)

        """
        if onlyancient we use the max ancient cov instead of the random
        one if the max ancient cov is larger
        """
        self.loc[
            (self["onlyAncient"] == True)
            & (self["max_ancient_cov_allowed"] > self["coverage_ancient_rnd"]),
            "coverage_ancient",
        ] = self["max_ancient_cov_allowed"]

        """
        if onlyancient we use get rnd ancient cov instead of the max cov
        one if the max ancient cov is smaller
        """
        self.loc[
            (self["onlyAncient"] == True)
            & (self["max_ancient_cov_allowed"] <= self["coverage_ancient_rnd"]),
            "coverage_ancient",
        ] = self["coverage_ancient_rnd"]

        """
        If onlyancient is False we use the random cov
        """
        self.loc[
            (self["onlyAncient"] == False),
            "coverage_ancient",
        ] = self["coverage_ancient_rnd"]

        """
        if we have a file specifiying cov and type of cov we use it
        """
        self.loc[
            self["Taxon"].isin(random_genomes), "coverage_ancient"
        ] = self.Taxon.map(mapping)

        """
        If the assigned coverage for onlyAncient is larger than the max
        cov we assign max
        """
        self.loc[
            (self["onlyAncient"] == True)
            & (self["coverage_ancient"] > self["max_ancient_cov_allowed"]),
            "coverage_ancient",
        ] = self["max_ancient_cov_allowed"]

        """
        If the assigned coverage for onlyAncient == False is larger than the max
        cov we assign the max
        """
        self.loc[
            (self["onlyAncient"] == False)
            & (self["coverage_ancient"] > self["max_ancient_cov_allowed"]),
            "coverage_ancient",
        ] = self["max_ancient_cov_allowed"]

        self.loc[
            self["onlyAncient"].isnull(),
            "coverage_ancient",
        ] = 0
    # random_genomes = withAncient[withAncient['Taxon'].isin(random_genomes)].copy()
    # self.loc[self["Taxon"].isnull(), "coverage_ancient"] = 0

    # print(self)
    return self


def applyParallel(dfGrouped, func, nproc, parms):
    with Pool(nproc) as p:
        func = partial(func, parms=parms)
        ret_list = p.map(func, [group for name, group in dfGrouped])
    return pd.concat(ret_list)


def _get_read_abundances(self, n_reads):
    # df = pd.DataFrame()
    dfs = []
    data_dict = OrderedDict()
    for p in self:
        if p.fragments_ancient is None:
            data_dict["comm"] = str(p.comm)
            data_dict["taxon"] = p.taxon
            data_dict["frag_type"] = "ancient"
            data_dict["seq_depth"] = 0
            _df = pd.DataFrame([data_dict], columns=data_dict.keys())
            # df = df.append(_df, ignore_index=True, sort=False)
            # df = f.concat_df([df, _df])
            dfs.append(_df)
        else:
            data_dict["comm"] = str(p.comm)
            data_dict["taxon"] = p.taxon
            data_dict["frag_type"] = "ancient"
            data_dict["seq_depth"] = p.fragments_ancient["seq_depth"]
            _df = pd.DataFrame([data_dict], columns=data_dict.keys())
            # df = df.append(_df, ignore_index=True, sort=False)
            # df = f.concat_df([df, _df])
            dfs.append(_df)
        if p.fragments_modern is None:
            data_dict["comm"] = str(p.comm)
            data_dict["taxon"] = p.taxon
            data_dict["frag_type"] = "modern"
            data_dict["seq_depth"] = 0
            _df = pd.DataFrame([data_dict], columns=data_dict.keys())
            # df = df.append(_df, ignore_index=True, sort=False)
            # df = f.concat_df([df, _df])
            dfs.append(_df)
        else:
            data_dict["comm"] = str(p.comm)
            data_dict["taxon"] = p.taxon
            data_dict["frag_type"] = "modern"
            data_dict["seq_depth"] = p.fragments_modern["seq_depth"]
            _df = pd.DataFrame([data_dict], columns=data_dict.keys())
            # df = df.append(_df, ignore_index=True, sort=False)
            # df = f.concat_df([df, _df])
            dfs.append(_df)
    df = f.concat_df(dfs)
    df["seq_depth_rounded"] = df.groupby("comm")["seq_depth"].transform(
        round_to_nreads, n_reads=n_reads
    )
    return df


def refine_abundances(data, n_reads):
    df = _get_read_abundances(data, n_reads=n_reads)
    df_taxo = (
        df.groupby(["comm", "taxon"])["seq_depth_rounded"]
        .sum()
        .reset_index(name="seq_depth")
    )

    for p in data:
        if p.fragments_ancient is not None:
            df_sub = df.loc[
                (df["comm"] == str(p.comm))
                & (df["taxon"] == p.taxon)
                & (df["frag_type"] == "ancient")
            ]
            if p.coverage_enforced is not True:
                p.fragments_ancient["seq_depth"] = int(
                    df_sub["seq_depth_rounded"].values
                )
        if p.fragments_modern is not None:
            df_sub = df.loc[
                (df["comm"] == str(p.comm))
                & (df["taxon"] == p.taxon)
                & (df["frag_type"] == "modern")
            ]
            if p.coverage_enforced is not True:
                p.fragments_modern["seq_depth"] = int(
                    df_sub["seq_depth_rounded"].values
                )

        seq_depth = df_taxo.loc[
            (df_taxo["comm"] == str(p.comm)) & (df_taxo["taxon"] == p.taxon)
        ]["seq_depth"].values
        p.seq_depth = int(seq_depth)
    return data


def rand_isAncient(row):
    """
    Define if a genome in a community is onlyAncient or not
    """
    if row["max_ancient_cov_allowed"] < row["coverage_ancient"]:
        return str(None)
    else:
        return str(
            bool(np.random.choice([True, False], p=[0.4, 0.6], size=1, replace=True))
        )


# Modified from https://stackoverflow.com/q/25271388
def round_to_nreads(number_set, n_reads, digit_after_decimal=0):
    """
    This function take a list of number and return a list of percentage, which represents the portion of each number in sum of all numbers
    Moreover, those percentages are adding up to 100%!!!
    Notice: the algorithm we are using here is 'Largest Remainder'
    The down-side is that the results won't be accurate, but they are never accurate anyway:)
    """
    unround_numbers = [
        x / float(sum(number_set)) * n_reads * 10**digit_after_decimal
        for x in number_set
    ]
    decimal_part_with_index = sorted(
        [(index, unround_numbers[index] % 1) for index in range(len(unround_numbers))],
        key=lambda y: y[1],
        reverse=True,
    )
    remainder = n_reads * 10**digit_after_decimal - sum(
        [int(x) for x in unround_numbers]
    )
    index = 0
    while remainder > 0:
        unround_numbers[decimal_part_with_index[index][0]] += 1
        remainder -= 1
        index = (index + 1) % len(number_set)
    return [int(x) / float(10**digit_after_decimal) for x in unround_numbers]


log = logging.getLogger("my_logger")


def get_ancient_genomes(args):
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    # simulating reads
    global debug
    if args.debug:
        debug = True
    else:
        debug = False

    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )

    # sys.excepthook = exceptionHandler

    # simulating reads
    cfg = vars(args)
    cfg = f.validate_schema(cfg, d.schema_init_ag, debug)
    config = f.get_config(config=cfg["config"], schema=d.ag_schema_config, debug=debug)

    # Test that both genome composition of abundance table are given
    if not config["abund-table"] and not config["genome-comp"]:
        raise IOError(
            "You need to provide an abundance table and/or genome composition table"
        )

    if config["abund-table"] and not config["genome-comp"]:
        abund_table = config["abund-table"]
        genome_comp = False
        filename = abund_table
    elif config["genome-comp"] and not config["abund-table"]:
        abund_table = False
        genome_comp = config["genome-comp"]
        filename = genome_comp
    else:
        genome_comp = config["genome-comp"]
        abund_table = config["abund-table"]
        filename = abund_table

    if config["library"] == "single":
        library = "se"
    else:
        library = "pe"
    seqSys = config["seqSys"]

    read_len = config["read-len"]

    max_len = d.seqSys_props[seqSys][library]

    # log.info('Set maximum fragment length to {}nt based on the seqSys: {}...'.format(
    #     max_len, seqSys))

    # Define different variables
    nproc = int(config["cpus"])
    prop = float(config["prop-ancient-comm"])
    min_len = config["min-length"]
    dist = config["length-dist"]
    mode_len_ancient = config["mode-len-ancient"]
    mode_len_modern = config["mode-len-modern"]
    n_reads = config["num-reads"]
    coverage = {"min": config["coverage"][0], "max": config["coverage"][1]}
    read_length_freqs = config["read-length-freqs"]

    if read_len > max_len:
        raise ValueError(f"Read length longer than the one allowed by {seqSys}...")

    if read_len < min_len:
        raise ValueError(f"Read length shorter than the minimum defined ({min_len})...")

    log.info("Loading genome paths...")
    genome_table = f.load_genome_table(config["genome-table"], nproc)

    # load tables

    if genome_comp:
        log.info("Loading genome composition table...")
        genome_comp = pd.read_csv(genome_comp, sep="\t")
        genome_comp = genome_comp.rename(columns={"Coverage": "coverage_ancient"})
    else:
        genome_comp = pd.DataFrame(columns=["Taxon", "coverage_ancient", "onlyAncient"])

    if abund_table:
        log.info("Loading abundace table...")
        abund_table = SimReads.load_abund_table(abund_table)
        df = abund_table.merge(genome_table, on=["Taxon"])
        if "Accession" in df.columns:
            genomes = df[
                [
                    "Community",
                    "Accession",
                    "Taxon",
                    "Perc_rel_abund",
                    "Rank",
                    "Genome_size",
                ]
            ].to_dict("records")
        else:
            genomes = df[
                ["Community", "Taxon", "Perc_rel_abund", "Rank", "Genome_size"]
            ].to_dict("records")
    else:
        df = genome_comp.merge(genome_table, on=["Taxon"])
        if "Accession" in df.columns:
            genomes = df[
                [
                    "Community",
                    "Accession",
                    "Taxon",
                    "coverage_ancient",
                    "Genome_size",
                ]
            ].to_dict("records")
        else:
            genomes = df[
                ["Community", "Taxon", "coverage_ancient", "Genome_size"]
            ].to_dict("records")

    if read_length_freqs:
        log.info("Loading read length frequencies from JSON...")
        read_length_freqs = f.load_read_length_freqs(read_length_freqs)

    # We need to add to the taxon the proportion of ancient
    # We need to take Community,
    log.info("Selecting random genomes...")
    if "Perc_rel_abund" in list(df.columns):
        df_anc = df[
            ["Community", "Taxon", "Perc_rel_abund", "Rank", "Genome_size"]
        ].copy()
    else:
        df_anc = df[["Community", "Taxon", "coverage_ancient", "Genome_size"]].copy()

    # df_anc['Perc_rel_abund'] = df_anc.groupby('Community')['Perc_rel_abund'].transform(round_to_100_percent)
    grouped = df_anc.groupby("Community")

    # Read from config file
    # Community Taxon onlyAncient Fold

    parms = {
        "prop": prop,
        "coverage": coverage,
        "n_reads": n_reads,
        "library": library,
        "mode_len_ancient": mode_len_ancient,
        "read_len": read_len,
        "rank_cutoff": 0.1,
        "genome_comp": genome_comp,
    }

    if debug is True:
        selected_genomes = applyParallel(
            grouped, func=select_genomes, nproc=1, parms=parms
        )
    else:
        selected_genomes = applyParallel(
            grouped, func=select_genomes, nproc=nproc, parms=parms
        )

    if read_length_freqs:
        log.info(f"Generating random fragment lengths from JSON file")
        log.info(f"::: Using empirical read lengths")
    else:
        log.info(f"Generating random fragment lengths from {dist} distribution")
        log.info(
            f"::: Using read lengths with modes [{mode_len_ancient}/{mode_len_modern}]"
        )

    # If there's not abundace table lets recalcultae them using the selected genomes
    if not "Perc_rel_abund" in list(df.columns):
        log.info("Calculating abundace table...")
        selected_genomes["Perc_rel_abund"] = 100 * (
            selected_genomes["seq_depth"]
            / selected_genomes.groupby("Community")["seq_depth"].transform("sum")
        )
        for genome in genomes:
            # extract value from dataframe where taxon and community match
            rel_abun = selected_genomes[
                (selected_genomes["Community"] == genome["Community"])
                & (selected_genomes["Taxon"] == genome["Taxon"])
            ]["Perc_rel_abund"].values[0]
            genome.update({"Perc_rel_abund": rel_abun})

    func = partial(
        process_genome,
        dist=dist,
        mode_len_ancient=mode_len_ancient,
        mode_len_modern=mode_len_modern,
        min_len=min_len,
        max_len=read_len,
        coverage=coverage,
        library=library,
        read_length_freqs=read_length_freqs,
        # seqSys=seqSys,
        selected_genomes=selected_genomes,
        n_reads=n_reads,
        genome_comp=genome_comp,
    )
    if debug is True:
        ancient_genomes_data = list(map(func, genomes))
    else:
        p = Pool(nproc)
        ancient_genomes_data = list(
            tqdm.tqdm(
                p.imap_unordered(func, genomes),
                total=len(genomes),
                leave=False,
                ncols=100,
                desc=f"Genomes processed",
            )
        )
        p.close()
        p.join()
        # ancient_genomes_data = p.map(func, genomes)
    if "Perc_rel_abund" in list(df.columns):
        file_name = Path(filename.name).stem
    else:
        file_name = Path(
            filename.name.replace("genome-compositions", "communities")
        ).stem
    dir_name = Path(filename.name).parent

    log.info("Refinning read abundances...")
    ancient_genomes = {}
    ancient_genomes["data"] = refine_abundances(ancient_genomes_data, n_reads=n_reads)
    # ancient_genomes["data"] = ancient_genomes_data
    ancient_genomes["experiment"] = {
        "library": library,
        "read_length": read_len,
        "seqSys": seqSys,
        "n_reads": n_reads,
        "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    df_tsv = _get_read_abundances(ancient_genomes_data, n_reads=n_reads)

    suffix = ".json"
    out_file = Path(dir_name, f"{file_name}{suffix}")

    #
    suffix_tsv = ".tsv"
    out_file_tsv = Path(dir_name, f"{file_name}_read-abundances{suffix_tsv}")

    df_tsv.to_csv(path_or_buf=out_file_tsv, sep="\t", header=True, index=False)

    ancient_genomes_json = json.dumps(
        ancient_genomes, default=obj_dict, ensure_ascii=False, indent=4
    )
    # ancient_genomes_json = jsonpickle.encode(ancient_genomes)
    log.info(f"Saving results to: {out_file} & {out_file_tsv}...")

    with open(out_file, "w", encoding="utf-8") as outfile:
        print(ancient_genomes_json, file=outfile)
