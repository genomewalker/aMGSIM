import numpy as np
import ruamel.yaml
from collections import OrderedDict, ChainMap
from distutils.spawn import find_executable
from functools import partial, reduce
from itertools import chain
import sys
import re
from schema import Schema, SchemaError
import gzip
from mimetypes import guess_type
from Bio import SeqIO
import pandas as pd
from multiprocessing import Pool
import gzip
import shutil
import json
import ujson


def str2dict(s):
    """Parsing string (format: 'item:value,item:value')
    to create a dict object.
    """
    if hasattr(s, "split"):
        l = re.split("[:,]", s)
        try:
            return {k.lower(): float(v) for k, v in zip(l[0::2], l[1::2])}
        except TypeError:
            msg = "Coverage parameter values must be ints or floats."
            raise TypeError(msg)
    else:
        return s


def validate_schema(self, schema, debug):
    schema = Schema(schema)
    try:
        self = schema.validate(self)
    except SchemaError as e:
        if debug:
            raise SchemaError("Schema error: {}".format(e)) from None
            # traceback gets printed
        else:
            print("{}: {}".format(type(e).__name__, e))
            sys.exit(1)
    return self


def merge_args(self, conf_args, orig_args):
    """Return new dict with args, and then conf_args merged in.
    Make sure that any keys in conf_args are also present in args
    """
    args = {}
    for k in conf_args.keys():
        if k not in orig_args:
            print("ERROR: Configuration file has unknown option %s" % k)
            sys.exit(-1)

    args.update(orig_args)
    args.update(conf_args)
    return args


def ordered_load(stream, Loader=ruamel.yaml.Loader, object_pairs_hook=OrderedDict):
    """Helper function to allow yaml load routine to use an OrderedDict instead of regular dict.
    This helps keeps things sane when ordering the runs and printing out routines
    """

    class OrderedLoader(Loader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))

    OrderedLoader.add_constructor(
        ruamel.yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_mapping
    )
    return ruamel.yaml.load(stream, OrderedLoader)


def check_correct_class(key, var_type):
    if key is not None:
        if not isinstance(key, var_type):
            raise TypeError(
                "Parameter {} has wrong type {}. It should be {}".format(
                    key, type(key).__name__, var_type.__name__
                )
            )


def check_config_key_exec(key, params, default, var_type):
    check_correct_class(params[key], var_type)

    if (key in params) and (params[key] != ""):
        key = params[key]
    else:
        key = default[key]
    check_exec(key)
    return key


def check_config_key(key, params, default, var_type):
    check_correct_class(params[key], var_type)
    if (key in params) and (params[key] != "") and (params[key] is not None):
        key = params[key]
    else:
        key = default[key]
    return key


def check_config_key_file(key, params, default, var_type, schema, debug):
    check_correct_class(params[key], var_type)
    if (key in params) and (params[key] != "") and (params[key] is not None):
        file = {key: params[key]}
        file = validate_schema(file, schema, debug)
    else:
        raise IOError("{} is obligatory".format(key))
    # check schema and return dict
    return file[key]


def check_exec(self):
    if find_executable(self) == None:
        raise IOError("Cannot find executable: {}".format(self))


def get_config(config, schema, debug):
    """Helper function to get the basic config params"""
    parameters = ordered_load(config)
    config = validate_schema(dict(parameters), schema, debug)
    return config


# def get_art_params(config):
#     """ Helper function to get the parameters for ART
#         Parameters accepted:
#             paired,
#             len,
#             mflen,
#             sdev,
#             amplicon,
#             rndSeed,
#             seqSys
#     """
#     parameters = ordered_load(config)
#     new_parameters = parameters.copy()
#     art_parameters = {}

#     if (('paired' in new_parameters) and (new_parameters['paired'] is True)):
#         art_parameters['paired'] = ""
#     else:
#         new_parameters.pop('paired', None)
#     if (('seqSys' in new_parameters) and (new_parameters['seqSys'] in d.seqSys_props)):
#         art_parameters['seqSys'] = new_parameters['seqSys']
#     else:
#         print("Falling back to default seqSys HS25")
#         art_parameters['seqSys'] = "HS25"
#     if ('len' in new_parameters):
#         # Check that length is smaller than seqSys
#         if (new_parameters['paired'] is True):
#             end_type = 'pe'
#         else:
#             end_type = 'se'
#         seqSys_rl = d.seqSys_props[art_parameters['seqSys']][end_type]
#         if (new_parameters['len'] > seqSys_rl):
#             print("Read length too long, falling back to default {} read length of {}bp".format(
#                 art_parameters['seqSys'], seqSys_rl))
#             art_parameters['len'] = seqSys_rl
#         else:
#             art_parameters['len'] = new_parameters['len']
#     if (('rndSeed' in new_parameters) and (isinstance(new_parameters['rndSeed'], int))):
#         art_parameters['rndSeed'] = new_parameters['rndSeed']
#     else:
#         art_parameters['rndSeed'] = 12345

#     return(art_parameters)


# def get_fragSim_params(config, debug):
#     """ Helper function to get the parameters for fragSim
#         Parameters accepted:
#             comp: --comp,
#             dist: --dist,
#             norev: --norev,
#             case: --case,
#             frag_length_file: -l but for each sample,
#             default_length: -l all sample same length,
#             frag_distribution_file: -f for each sample,
#             frag_lognormal_distribution_file: --loc --scale for each sample,
#             default_loc: --loc all samples same location,
#             default_scale: --scale all samples same location
#     """
#     parameters = ordered_load(config)
#     new_parameters = parameters.copy()
#     fragSim_parameters = {}

#     if (("comp" in new_parameters) and (new_parameters['comp'] is not False)):
#         fragSim_parameters['comp'] = new_parameters['comp']
#     if (("dist" in new_parameters) and (new_parameters['dist'] is not False)):
#         fragSim_parameters['dist'] = new_parameters['dist']
#     if (("norev" in new_parameters) and (new_parameters['norev'] is True)):
#         fragSim_parameters['norev'] = ""
#     if (("case" in new_parameters) and (new_parameters['case'] is True)):
#         fragSim_parameters['case'] = ""
#     if (("ancient_genomes" in new_parameters) and (new_parameters['ancient_genomes'] is not False)):
#         try:
#             with open(new_parameters['ancient_genomes']) as f:
#                 fragSim_parameters['ancient_genomes'] = json.load(f)
#         except EnvironmentError:  # parent of IOError, OSError *and* WindowsError where available
#             raise EnvironmentError("Cannot open {}".format(
#                 new_parameters['ancient_genomes']))

#     # if (("default_fragment_length" in new_parameters) and (new_parameters['default_fragment_length'] is not False)):
#     #     fragSim_parameters['l'] = new_parameters['default_fragment_length']

#     return(fragSim_parameters)


# def get_deamSim_params(config, debug):
#     """ Helper function to get the parameters for deamSim
#         Parameters accepted:
#             mapdamage  Read the miscorporation file [mis.txt] produced by mapDamage,
#             damage  For the Briggs et al. 2007 model The parameters must be comma-separated e.g.: -damage 0.03,0.4,0.01,0.7
#     """
#     parameters = ordered_load(config)
#     new_parameters = parameters.copy()
#     deamSim_parameters = {}

#     if ((new_parameters['mapdamage'] is not False) and (new_parameters['damage'] is not False)):
#         raise Exception("You cannot specify -mapdamage and -damage options")

#     if (new_parameters['mapdamage'] is not False):
#         deamSim_parameters['mapdamage'] = new_parameters['mapdamage']
#     if "damage" in new_parameters is not False:
#         for k in new_parameters['damage']:
#             if k < 0:
#                 raise Exception("Briggs parameters should be positive")

#         deamSim_parameters['damage'] = ",".join(
#             map(str, new_parameters['damage']))

#     return(deamSim_parameters)


def power_neg(*args, **kwargs):
    return 1 - np.random.power(*args, **kwargs)


def get_freq_dist(dist, params):
    try:
        distFunc = getattr(np.random, dist)
    except AttributeError:
        msg = 'Distribution "{}" is not supported'.format(dist)
        raise AttributeError(msg)

    if dist == "power":
        distFunc = power_neg
    # expected_length = lognorm.expect(lambda x :max(x,min_length), args=[s],
    #                                     loc=loc, scale=scale, lb=-np.inf, ub=np.inf)
    try:
        return partial(distFunc, **params)
    except TypeError:
        param_str = [":".join([str(k), str(v)]) for k, v in params.items()]
        param_str = ",".join(param_str)
        msg = 'Params "{}" do not work with function "{}"'.format(param_str, dist)
        raise TypeError(msg)


# From https://stackoverflow.com/a/41477069


def lognorm_params(mode, stddev):
    """
    Given the mode and std. dev. of the log-normal distribution, this function
    returns the shape and scale parameters for scipy's parameterization of the
    distribution.
    """
    p = np.poly1d([1, -1, 0, 0, -((stddev / mode) ** 2)])
    r = p.roots
    sol = r[(r.imag == 0) & (r.real > 0)].real
    shape = np.sqrt(np.log(sol))
    scale = mode * sol
    return shape, scale


# functions
def tidy_taxon_names(x):
    """Remove special characters from taxon names"""
    x = re.sub(r"[()\/:;, ]+", "_", x)
    return x


def load_genome_table(in_file, nproc=1):
    """Loading genome table
    Parameters
    ----------
    in_file : str
        input file path
    """
    nproc = int(nproc)
    df = pd.read_csv(in_file, sep="\t")

    ## check headers
    diff = set(["Taxon", "Fasta"]) - set(df.columns.values)
    if len(diff) > 0:
        diff = ",".join(diff)
        raise ValueError("Cannot find table columns: {}".format(diff))

    # getting genome sizes
    if nproc > 1:
        p = Pool(nproc)
        df["Genome_size"] = p.map(_genome_size, [x for i, x in df.iterrows()])
        p.close()
        p.join()
    else:
        df["Genome_size"] = [_genome_size(x) for i, x in df.iterrows()]

    # tidy taxon names
    df["Taxon"] = df["Taxon"].astype(str).apply(tidy_taxon_names)
    return df


def _genome_size(x):
    """Get the total bp for the genome sequence of all genomes

    Parameters
    ----------
    x : pd.Series
       Series that includes 'Fasta' in the index
    """
    bp = 0
    encoding = guess_type(x["Fasta"])[1]  # uses file extension
    _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open
    with _open(x["Fasta"]) as f:
        for i, record in enumerate(SeqIO.parse(f, "fasta")):
            bp += len(record.seq)
    return bp


def load_read_length_freqs(file):
    file = ujson.load(file)
    return reduce(lambda d, src: d.update(src) or d, file, {})


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, "rb")
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = "plain"
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("Error: cannot use bzip2 format - use gzip instead")
        sys.exit("Error: cannot use zip format - use gzip instead")
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == "gz":
        return gzip.open
    else:  # plain text
        return open


def check_filter_conditions(filt_dict, default_filter_values, filters):
    """A function that checks if a filter condition is valid.

    Args:
        filt_dict (dict): A dictionary with different filter conditions.
        default_filter_values (list): A list default filter values.

    Raises:
        SchemaError: If the filter is not correct

    Returns:
        dict: A dictionary with the filter conditions.
    """
    # remove keys if they are not in the filter dictionary
    filt_dict = {k: filt_dict[k] for k in filters if k in filt_dict}
    # check if the filter values are correct
    # if filter is empty, use default values
    if not filt_dict:
        filt_dict = {**default_filter_values, **filt_dict}
    if all(value >= 0 for value in filt_dict.values()):
        return filt_dict
    else:
        raise SchemaError()


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def concat_df(frames):
    COLUMN_NAMES = frames[0].columns
    df_dict = dict.fromkeys(COLUMN_NAMES, [])
    for col in COLUMN_NAMES:
        extracted = (frame[col] for frame in frames)
        # Flatten and save to df_dict
        df_dict[col] = fast_flatten(extracted)
    df = pd.DataFrame.from_dict(df_dict)[COLUMN_NAMES]
    return df


def initializer(init_data):
    global parms
    parms = init_data


# from https://stackoverflow.com/questions/53751050/python-multiprocessing-understanding-logic-behind-chunksize/54032744#54032744
def calc_chunksize(n_workers, len_iterable, factor=4):
    """Calculate chunksize argument for Pool-methods.
    Resembles source-code within `multiprocessing.pool.Pool._map_async`.
    """
    chunksize, extra = divmod(len_iterable, n_workers * factor)
    if extra:
        chunksize += 1
    return chunksize


def unzip_files(infile, outfile):
    with gzip.open(infile, "rb") as f_in:
        with open(outfile, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
