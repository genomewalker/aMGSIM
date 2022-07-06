import argparse
import sys
import gzip
import os
import shutil
import logging
from aMGSIM import __version__
import shutil


class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            values = list(map(float, values.split(",")))
        except ValueError:
            raise argparse.ArgumentError(self, "invalid values: %r" % values)
        # check all elements in list are greater than 0
        if not all(x > 0 for x in values):
            raise argparse.ArgumentError(self, "all elements must be greater than 0")
        setattr(namespace, self.dest, values)


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


filters = ["breadth", "depth", "depth_evenness", "breadth_expected_ratio"]


# From https://stackoverflow.com/a/59617044/15704171
def convert_list_to_str(lst):
    n = len(lst)
    if not n:
        return ""
    if n == 1:
        return lst[0]
    return ", ".join(lst[:-1]) + f" or {lst[-1]}"


def check_filter_values(val, parser, var):
    value = str(val)
    if value in filters:
        return value
    else:
        parser.error(
            f"argument {var}: Invalid value {value}. Filter has to be one of {convert_list_to_str(filters)}"
        )


def check_values(val, minval, maxval, parser, var):
    value = float(val)
    if value < minval or value > maxval:
        parser.error(
            f"argument {var}: Invalid value value. Range has to be between {minval} and {maxval}!"
        )
    return value


# From: https://note.nkmk.me/en/python-check-int-float/
def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


# function to check if the input value has K, M or G suffix in it
def check_suffix(val, parser, var):
    units = ["K", "M", "G"]
    unit = val[-1]
    value = int(val[:-1])

    if is_integer(value) & (unit in units) & (value > 0):
        return val
    else:
        parser.error(
            "argument %s: Invalid value %s. Memory has to be an integer larger than 0 with the following suffix K, M or G"
            % (var, val)
        )


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


# From: https://stackoverflow.com/a/11541450
def is_valid_file(parser, arg, var):
    if not os.path.exists(arg):
        parser.error("argument %s: The file %s does not exist!" % (var, arg))
    else:
        return arg


def is_executable(parser, arg, var):
    """Check whether `name` is on PATH and marked as executable."""
    if shutil.which(arg):
        return arg
    else:
        parser.error("argument %s: %s not found in path." % (var, arg))


def create_folder_and_delete(path):
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        log.info(f"Folder {path} already exists. Deleting.")
        try:
            shutil.rmtree(path)
            os.makedirs(path)
        except OSError as e:
            print("Error: %s : %s" % (path, e.strerror))


def delete_folder(path):
    if os.path.exists(path):
        log.info(f"Deleting folder {path}.")
        try:
            shutil.rmtree(path)
        except OSError as e:
            print("Error: %s : %s" % (path, e.strerror))
    else:
        log.info(f"Path {path} does not exist.")


def create_folder(path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError as e:
            print("Error: %s : %s" % (path, e.strerror))
    else:
        log.info(f"Folder {path} already exists. Skipping creation.")


defaults = {
    "threads": 1,
    "processes": 1,
    "gene_min_length": 30,
    "tmp": ".tmp",
    "output": "protein-analysis",
    "cov_suffix": "multicov",
    "coverages": [0.1, 0.5, 1, 5, 10, 20, 50, 100],
}

help_msg = {
    "config": "Config file ",
    "genome_table": "A TSV file containing the columns Taxon and Fasta",
    "abund_read_file": "A JSON file containing the abundance of reads",
    "genome_compositions": "A TSV file containing the genome compositions",
    "processes": "Number of processes to use",
    "threads": "Number of threads to use in each process",
    "output": "Output folder",
    "cov_suffix": "Suffix added to the filename for the multi coverage file",
    "coverages": "Coverage values for the multi coverage",
    "tmp": "Temporary folder",
    "gene_min_length": "Minimum gene length predicts by Prodigal",
    "help": "Help message",
    "debug": f"Print debug messages",
    "version": f"Print program version",
}

log = logging.getLogger("my_logger")


def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="A simple tool to simulate ancient metagenomes for multiple synthetic communities",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
    )

    # Same subparsers as usual
    sub_parsers = parser.add_subparsers(
        help="positional arguments",
        dest="action",
    )

    # Create parent subparser. Note `add_help=False` and creation via `argparse.`
    parent_parser = argparse.ArgumentParser(add_help=False)
    optional = parent_parser._action_groups.pop()
    optional.add_argument(
        "--debug", dest="debug", action="store_true", help=help_msg["debug"]
    )

    # create the parser sub-commands
    parser_communities = sub_parsers.add_parser(
        "communities",
        help="This will call the MGSIM communities command to generate random taxon abundances",
        parents=[parent_parser],
    )
    parser_estimate = sub_parsers.add_parser(
        "estimate",
        help="Estimate coverage, depth and other properties for each genome in a sample processed with bam-filter",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_multicov = sub_parsers.add_parser(
        "multicov",
        help="Create different coverage versions of a genome compositions file",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_ag = sub_parsers.add_parser(
        "ancient-genomes",
        help="Estimate coverage, depth and other properties for each genome in each synthetic community",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_ar = sub_parsers.add_parser(
        "ancient-reads",
        help="Simulate ancient reads for each taxon in each synthetic community",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_pa = sub_parsers.add_parser(
        "protein-analysis",
        help="Tracking damage to the codon positions of each simulated read",
        parents=[parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    estimate_required_args = parser_estimate.add_argument_group("required arguments")
    estimate_required_args.add_argument(
        "config",
        type=lambda x: is_valid_file(parser, x, "config"),
        help=help_msg["config"],
    )

    cov_required_args = parser_multicov.add_argument_group("required arguments")
    cov_optional_args = parser_multicov.add_argument_group("optional arguments")
    cov_required_args.add_argument(
        "genome_compositions",
        type=lambda x: is_valid_file(parser, x, "genome_compositions"),
        help=help_msg["genome_compositions"],
    )
    cov_optional_args.add_argument(
        "--suffix",
        type=str,
        default=defaults["cov_suffix"],
        metavar="STR",
        dest="cov_suffix",
        help=help_msg["cov_suffix"],
    )
    cov_optional_args.add_argument(
        "--coverages",
        default=defaults["coverages"],
        dest="coverages",
        action=SplitArgs,
        help=help_msg["coverages"],
    )

    ag_required_args = parser_ag.add_argument_group("required arguments")
    ag_required_args.add_argument(
        "config",
        type=lambda x: is_valid_file(parser, x, "config"),
        help=help_msg["config"],
    )
    ar_required_args = parser_ar.add_argument_group("required arguments")
    ar_required_args.add_argument(
        "genome_table",
        type=lambda x: is_valid_file(parser, x, "genome_table"),
        help=help_msg["genome_table"],
    )
    ar_required_args.add_argument(
        "config",
        type=lambda x: is_valid_file(parser, x, "config"),
        help=help_msg["config"],
    )

    pa_required_args = parser_pa.add_argument_group("required arguments")
    pa_optional_args = parser_pa.add_argument_group("optional arguments")
    pa_required_args.add_argument(
        "abund_read_file",
        type=lambda x: is_valid_file(parser, x, "abund_read_file"),
        help=help_msg["abund_read_file"],
    )
    pa_optional_args.add_argument(
        "--threads",
        type=lambda x: int(
            check_values(x, minval=1, maxval=1000, parser=parser, var="--threads")
        ),
        dest="threads",
        default=1,
        help=help_msg["threads"],
    )
    pa_optional_args.add_argument(
        "--processes",
        type=lambda x: int(
            check_values(x, minval=1, maxval=1000, parser=parser, var="--processes")
        ),
        dest="processes",
        default=1,
        help=help_msg["processes"],
    )
    pa_optional_args.add_argument(
        "--output",
        type=str,
        default=defaults["output"],
        metavar="STR",
        dest="output",
        help=help_msg["output"],
    )
    pa_optional_args.add_argument(
        "--tmp",
        type=str,
        default=defaults["tmp"],
        metavar="STR",
        dest="tmp",
        help=help_msg["tmp"],
    )
    pa_optional_args.add_argument(
        "--gene-min-length",
        type=lambda x: int(
            check_values(
                x, minval=1, maxval=100000, parser=parser, var="--gene-min-length"
            )
        ),
        default=defaults["gene_min_length"],
        dest="gene_min_length",
        help=help_msg["gene_min_length"],
    )

    args = parser.parse_args(None if sys.argv[1:] else ["-h"])

    if args.action is not None and len(sys.argv) == 2:
        if args.action != "communities":
            args = parser.parse_args([args.action, "-h"])
    return args
