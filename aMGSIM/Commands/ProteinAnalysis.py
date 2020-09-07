#!/usr/bin/env python

"""
ancient-genomes: simulating ancient reads

Usage:
  ancient-genomes [options] <config>
  ancient-genomes -h | --help
  ancient-genomes --version

Options:
  <config>              Config parameters
  -d --debug                 Debug mode (no subprocesses; verbose output)
  -h --help                  Show this screen.
  --version                  Show version.

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
from biolib.external.prodigal import Prodigal
# from biolib.external.execute import check_on_path
#import biolib

# Per sample generate all necessary to recreate ancient genomes

# {
#     'comm': {
#         'name': 'com1',
#         'Taxon': {
#             'name': 'name1',
#             'fragments': {
#                 'length': [10 20 30 40 60],
#                 'freq': [0.1 0.4 0.3 0.1 0.1]
#             },
#             'onlyAncient': False
#         }
#     }
# }


def exceptionHandler(exception_type, exception, traceback, debug_hook=sys.excepthook):
    '''Print user friendly error messages normally, full traceback if DEBUG on.
       Adapted from http://stackoverflow.com/questions/27674602/hide-traceback-unless-a-debug-flag-is-set
    '''
    if debug:
        print("\n*** Error:")
        # raise
        debug_hook(exception_type, exception, traceback)
    else:
        print("{}: {}".format(exception_type.__name__, exception))


sys.excepthook = exceptionHandler


# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def obj_dict(obj):
    return obj.__dict__


def process_genome(x, dist, mode_len_ancient, mode_len_modern, min_len, max_len, coverage, library, selected_genomes, n_reads):
    genome = Genome(
        genome=x,
        dist=dist,
        mode_len_ancient=mode_len_ancient,
        mode_len_modern=mode_len_modern,
        min_len=min_len,
        max_len=max_len,
        coverage=coverage,
        library=library,
        # seqSys=seqSys,
        selected_genomes=selected_genomes,
        n_reads=n_reads)
    return(genome)


class Genome:

    def __init__(self, genome, dist, mode_len_ancient, mode_len_modern, selected_genomes, min_len, max_len, coverage, library, n_reads, *args, **kwargs):
        # Community, Taxon, Perc_rel_abund, Rank, Fasta, Genome_size
        rnd_seed = int(np.random.RandomState().randint(0, 99999, size=1))
        np.random.seed(rnd_seed)

        self.comm = genome[0]
        self.taxon = genome[1]
        self.rel_abund = genome[2]
        self.genome_size = genome[5]
        # self.library = library
        # self.read_length = max_len
        # self.seqSys = seqSys
        # self.n_reads = n_reads

        selected_genomes_filt = selected_genomes[(
            selected_genomes['Community'] == self.comm) & (selected_genomes['Taxon'] == self.taxon)]

        if len(selected_genomes_filt.index) > 0:
            self.onlyAncient = bool(
                selected_genomes_filt.iloc[0]['onlyAncient'])
        else:
            self.onlyAncient = None

        if self.onlyAncient is not None:
            min_cov = coverage['min']
            max_cov = coverage['max']
            coverage_ancient = selected_genomes_filt.iloc[0]['coverage_ancient']

        # Get fragments
        # Get average size
        if self.onlyAncient:
            self.fragments_ancient = self.random_fragments(min_len=min_len,
                                                           max_len=max_len,
                                                           dist=dist,
                                                           mode_len=mode_len_ancient,
                                                           rnd_seed=rnd_seed)
            self.fragments_modern = None

            if library == "se":
                paired = 1
            else:
                paired = 2

            read_length_ancient = int(
                self.fragments_ancient['avg_len'] * paired)

            # Number max of reads genome in sample
            seq_depth = int(n_reads*((self.rel_abund/100)))

            seq_depth_ancient = int(
                (coverage_ancient * self.genome_size)/read_length_ancient)
            seq_depth_ancient_original = seq_depth_ancient
            seq_depth_ancient = seq_depth

            fold_ancient = (
                (seq_depth_ancient * read_length_ancient) / self.genome_size)

            if seq_depth_ancient > 0:
                self.fragments_ancient['seq_depth'] = seq_depth_ancient
                self.fragments_ancient['seq_depth_original'] = seq_depth_ancient_original
                self.fragments_ancient['fold'] = fold_ancient
                self.fragments_ancient['fold_original'] = coverage_ancient
            else:
                self.fragments_ancient = None

            self.seq_depth = seq_depth

        elif self.onlyAncient is None:
            self.fragments_ancient = None
            self.fragments_modern = self.random_fragments(min_len=min_len,
                                                          max_len=max_len,
                                                          dist=dist,
                                                          mode_len=mode_len_modern,
                                                          rnd_seed=rnd_seed)
            if library == "se":
                paired = 1
            else:
                paired = 2

            read_length_modern = int(
                self.fragments_modern['avg_len'] * paired)

            # Number max of reads genome in sample
            seq_depth_modern = int(n_reads*((self.rel_abund/100)))
            fold_modern = (
                (seq_depth_modern * read_length_modern) / self.genome_size)

            if seq_depth_modern > 0:
                self.fragments_modern['seq_depth'] = seq_depth_modern
                self.fragments_modern['fold'] = fold_modern
            else:
                self.fragments_modern = None
            self.seq_depth = seq_depth_modern

        else:
            self.fragments_ancient = self.random_fragments(min_len=min_len,
                                                           max_len=max_len,
                                                           dist=dist,
                                                           mode_len=mode_len_ancient,
                                                           rnd_seed=rnd_seed)
            self.fragments_modern = self.random_fragments(min_len=min_len,
                                                          max_len=max_len,
                                                          dist=dist,
                                                          mode_len=mode_len_modern,
                                                          rnd_seed=rnd_seed)
            if library == "se":
                paired = 1
            else:
                paired = 2
            read_length_ancient = int(
                self.fragments_ancient['avg_len'] * paired)
            read_length_modern = int(
                self.fragments_ancient['avg_len'] * paired)

            # Number max of reads genome in sample
            seq_depth = int(n_reads*((self.rel_abund/100)))

            seq_depth_ancient = int(
                (coverage_ancient * self.genome_size)/read_length_ancient)

            if seq_depth_ancient < seq_depth:
                seq_depth_modern = int(seq_depth - seq_depth_ancient)
            else:
                seq_depth_ancient = seq_depth
                seq_depth_modern = 0

            fold_modern = (
                (seq_depth_modern * read_length_modern) / self.genome_size)
            fold_ancient = (
                (seq_depth_ancient * read_length_ancient) / self.genome_size)
            if seq_depth_ancient > 0:
                self.fragments_ancient['seq_depth'] = seq_depth_ancient
                self.fragments_ancient['fold'] = fold_ancient
                self.fragments_ancient['fold_original'] = coverage_ancient
            else:
                self.fragments_ancient = None
            if seq_depth_modern > 0:
                self.fragments_modern['seq_depth'] = seq_depth_modern
                self.fragments_modern['fold'] = fold_modern
            else:
                self.fragments_modern = None
            self.seq_depth = seq_depth

    def random_fragments(self, min_len, max_len, dist, mode_len, rnd_seed):
        n_frags = int(1e6)
        # n_frags = max_len - min_len + 1
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
        freqs = list(lengths/sum(lengths))
        lens = np.multiply(list(lengths.index), freqs).sum()
        avg_len = float(lens/sum(freqs))
        # print("Lengths: {}", list(lengths.index))
        # print("Freqs: {}", list(lengths/sum(lengths)))
        # print("Sizes: {}".format(sizes))
        # print("Avg size: {}".format(avg_size))
        # exit()
        frags = {'fragments': {
            'length': list(lengths.index),
            'freq': freqs
        },
            'dist_params': {
            'mu': float(mu),
            'sigma': float(sigma),
            'rnd_seed': rnd_seed
        },
            'avg_len': avg_len
        }
        return(frags)

    def __str__(self):
        return str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item]) for item in sorted(self.__dict__)))

    def __repr__(self):
        return str(self)

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__)


def select_genomes(self, parms):
    prop = parms['prop']
    coverage = parms['coverage']
    n_reads = parms['n_reads']
    library = parms['library']
    mode_len_ancient = parms['mode_len_ancient']
    cov_max = coverage['max']
    read_len = parms['read_len']

    if library == "se":
        paired = 1
    else:
        paired = 2

    read_length_ancient = mode_len_ancient * paired

    rnd_seed = int(np.random.RandomState().randint(0, 99999, size=1))
    np.random.seed(rnd_seed)

    min_cov = coverage['min']
    max_cov = coverage['max']

    k = int(round(prop * self.shape[0]))
    k_filt = 0

    while (k_filt <= k):
        # print(i)

        coverage_ancient = np.round(np.random.uniform(
            min_cov, max_cov, size=self.shape[0]), 1)

        self['seq_depth'] = (n_reads * (self['Perc_rel_abund']/100))
        self['seq_depth'] = self['seq_depth'].astype(int)
        self['max_ancient_cov_allowed'] = (
            self['seq_depth'] * read_length_ancient)/self['Genome_size']
        self['coverage_ancient'] = coverage_ancient.tolist()

        self['onlyAncient'] = self.apply(rand_isAncient, axis=1)
        # seq_depth_ancient =
        # seq_depth_ancient_original = seq_depth_ancient
        # seq_depth_ancient = seq_depth

        # fold_ancient =

        withAncient = self[(self['onlyAncient'].notnull())]

        k_filt = withAncient.shape[0]

    pop = withAncient['Taxon'].tolist()
    # Estimate max number of reads as function of coverage and rel abun

    # Number max of reads genome in sample

    w = pd.DataFrame(withAncient['Perc_rel_abund'].apply(
        lambda x: np.arcsinh(100-x)))
    max_w = w['Perc_rel_abund'].sum()
    w = w['Perc_rel_abund'].apply(lambda x: x/max_w)

    random_genomes = np.random.choice(pop,
                                      replace=False,
                                      p=w,
                                      size=k)

   # random_genomes = withAncient[withAncient['Taxon'].isin(random_genomes)].copy()
    self.loc[~self['Taxon'].isin(random_genomes), 'onlyAncient'] = None
    self.loc[~self['Taxon'].isin(random_genomes), 'coverage_ancient'] = None
    #self['coverage_ancient'] = self['coverage_ancient'].where(self['onlyAncient'].notnull(), self['max_ancient_cov_allowed'])
    self['coverage_ancient'] = self['coverage_ancient'].where(
        self['onlyAncient'].isin([False]), self['max_ancient_cov_allowed'])
    #self = self.drop('Perc_rel_abund', axis=1)

    return(self)


def applyParallel(dfGrouped, func, nproc, parms):
    with Pool(nproc) as p:
        func = partial(func,
                       parms=parms)
        ret_list = p.map(func, [group for name, group in dfGrouped])
    return pd.concat(ret_list)


# def refine_abundances(x, data, selected_genomes):

#     data = data.merge(selected_genomes,  on=["Community", "Taxon"], how="left")
#     selected_genomes_filt = data[(
#         data['Community'] == x)]
#     print(selected_genomes_filt)
#     exit(0)
#     rnd_seed = int(np.random.RandomState().randint(0, 99999, size=1))
#     np.random.seed(rnd_seed)

#     if len(selected_genomes_filt.index) > 0:
#         onlyAncient = bool(
#             selected_genomes_filt.iloc[0]['onlyAncient'])
#     else:
#         self.onlyAncient = None

#     if self.onlyAncient is not None:
#         min_cov = coverage['min']
#         max_cov = coverage['max']
#         coverage_ancient = round(np.random.uniform(min_cov, max_cov), 1)

#     print(selected_genomes_filt)
#     exit(0)

# # Community, Taxon, Perc_rel_abund, Rank, Fasta, Genome_size

#     self.comm = genome[0]
#     self.taxon = genome[1]
#     self.rel_abund = genome[2]
#     self.genome_size = genome[5]
#     # self.library = library
#     # self.read_length = max_len
#     # self.seqSys = seqSys
#     # self.n_reads = n_reads

#     # data = list(set([s for s in ancient_genome_data]))
#     # for d in data:
#     #     print(d)

#     data = [d for d in ancient_genome_data if d.comm == x]
#     for d in data:
#         if d.fragments_ancient is not None:
#             if (d.fragments_ancient['fold'] > d.fragments_ancient['fold_original']):
#                 print("taxon:{} fold:{} fold_original:{} seq_depth:{} seq_depth_original:{}".format(d.taxon, d.fragments_ancient['fold'],
#                                                                                                     d.fragments_ancient[
#                                                                                                         'fold_original'], d.fragments_ancient['seq_depth'],
#                                                                                                     d.fragments_ancient['seq_depth_original']))

def rand_isAncient(row):
    if row['max_ancient_cov_allowed'] < row['coverage_ancient']:
        return None
    else:
        return bool(np.random.choice([True, False], p=[0.6, 0.4], size=1, replace=True))

# def predict_genes(x, cpus, output_dir, called_genes=False, translation_table=None, meta=False, closed_ends=False, verbose = True):
#     genome = [x[4]]
#     prodigal = Prodigal(cpus=cpus, verbose = verbose)
#     genes = prodigal.run(
#         genome_files=genome,
#         output_dir=output_dir, 
#         called_genes=False, 
#         translation_table=None, 
#         meta=False,
#         closed_ends=False)
#     return(genes)

def main(args):
    # simulating reads
    global debug
    debug = False

    if args['--debug']:
        debug = True

    # simulating reads
    args = f.validate_schema(args, d.schema_init_ag, debug)
    config = f.get_config(
        config=args['<config>'], schema=d.ag_schema_config, debug=debug)

    # Test that both genome composition of abundance table are given
    if config['abund_table'] and config['genome_comp']:
        raise IOError(
            "You cannot provide an abundance table and genome composition table at the same time.")

    if config['abund_table']:
        abund_table = config['abund_table']
        filename = abund_table
    elif config['genome_comp']:
        genome_comp = config['genome_comp']
        filename = genome_comp

    tmp_dir = config['tmp_dir']

    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)

    if config['single'] is True:
        library = 'se'
    else:
        library = 'pe'
    seqSys = config['seqSys']

    read_len = config['read-len']

    max_len = d.seqSys_props[seqSys][library]

    # logging.info('Set maximum fragment length to {}nt based on the seqSys: {}...'.format(
    #     max_len, seqSys))

    # Define different variables
    nproc = int(config['cpus'])
    prop = float(config['prop-ancient-comm'])
    min_len = config['min-length']
    dist = config['length-dist']
    mode_len_ancient = config['mode-len-ancient']
    mode_len_modern = config['mode-len-modern']
    n_reads = config['num-reads']
    coverage = {
        "min": config['coverage'][0],
        "max": config['coverage'][1]
    }

    enforce_cov = config['enforce-cov']

    if read_len > max_len:
        raise ValueError(
            "Read length longer than the one allowed by {}...".format(seqSys))

    if read_len < min_len:
        raise ValueError(
            "Read length shorter than the minimum defined ({})...".format(min_len))

    # load tables
    abund_table = SimReads.load_abund_table(abund_table)
    genome_table = SimReads.load_genome_table(config['genome_table'])
    df = abund_table.merge(genome_table, on=['Taxon'])
    genomes = df.values.tolist()
    # We need to add to the taxon the proportion of ancient
    # We need to take Community,
    logging.info('Selecting random genomes...')
    df_anc = df[['Community', "Taxon", "Perc_rel_abund", "Genome_size"]]
    grouped = df_anc.groupby('Community')

    # Read from config file
    # Community Taxon onlyAncient Fold

    parms = {
        'prop': prop,
        'coverage': coverage,
        'n_reads': n_reads,
        'library': library,
        'mode_len_ancient': mode_len_ancient,
        'read_len': read_len
    }
    selected_genomes = applyParallel(grouped,
                                     func=select_genomes,
                                     nproc=nproc,
                                     parms=parms)

    # refine abundances
    # Ancient genome coverage takes over community community composition
    # Revert to original genome coverage and split reads between all other modern taxa

    # comms = list(set(df["Community"].values))
    # print(comms)

    # func = partial(refine_abundances,
    #                data=df,
    #                selected_genomes=selected_genomes)

    # logging.info('Combining reads by sample...')

    # if debug is True:
    #     files = list(map(func, comms))
    # else:
    #     p = Pool(nproc)
    #     files = list(tqdm.tqdm(p.imap_unordered(
    #         func, comms), total=len(comms)))

    # print(df)
    # exit(0)

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

    prodigal = Prodigal(cpus = nproc, verbose = False)
    genes = prodigal.run(
        genome_files=set(df.Fasta.values),
        output_dir=output_dir, 
        called_genes=False, 
        translation_table=None, 
        meta=False,
        closed_ends=False)
    print(genes)
    exit(0)
    logging.info('Generating random fragment lengths from {} distribution with modes {} and {}...'.format(
        dist, mode_len_ancient, mode_len_modern))

    func = partial(process_genome,
                   dist=dist,
                   mode_len_ancient=mode_len_ancient,
                   mode_len_modern=mode_len_modern,
                   min_len=min_len,
                   max_len=read_len,
                   coverage=coverage,
                   library=library,
                   # seqSys=seqSys,
                   selected_genomes=selected_genomes,
                   n_reads=n_reads)

    if debug is True:
        ancient_genomes_data = list(map(func, genomes))
    else:
        p = Pool(nproc)
        ancient_genomes_data = list(
            tqdm.tqdm(p.imap_unordered(func, genomes), total=len(genomes)))
        # ancient_genomes_data = p.map(func, genomes)

    file_name = Path(filename.name).stem
    dir_name = Path(filename.name).parent

    ancient_genomes = {}
    ancient_genomes["data"] = ancient_genomes_data
    ancient_genomes["experiment"] = {
        "library": library,
        "read_length": read_len,
        "seqSys": seqSys,
        "n_reads": n_reads,
        "date": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    }

    suffix = ".json"
    out_file = Path(dir_name, file_name).with_suffix(suffix)

    ancient_genomes_json = json.dumps(
        ancient_genomes, default=obj_dict, ensure_ascii=False, indent=4)
    # ancient_genomes_json = jsonpickle.encode(ancient_genomes)
    logging.info('Saving results to: {}...'.format(out_file))

    with open(out_file, 'w', encoding='utf-8') as outfile:
        print(ancient_genomes_json, file=outfile)


def opt_parse(args=None):
    if args is None:
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
