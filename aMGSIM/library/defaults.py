# Constants
from schema import And, Optional, Or, Use
from multiprocessing import cpu_count
import os
import shutil
import datetime
from aMGSIM.library import functions as f
import random
import string
import re

# schema_init_ag = {
#     '<abund_table>': [Use(open, error='abund_table should be readable')],
#     '-p': Or(None, And(Use(float), lambda n: 0 < n < 1),
#              error='Proportion (-p) should be a float 0 < p < 1'),
#     object: object
# }


def generate_random_sample_name():
    """
    Generate a random sample name
    """
    return "".join(random.choices(string.ascii_lowercase + string.digits, k=10))


debug = False

FWD_ADPT = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG"
REV_ADPT = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT"

IUPACAmbiguousDNA = {
    "G",
    "A",
    "T",
    "C",
    "R",
    "Y",
    "W",
    "S",
    "M",
    "K",
    "H",
    "B",
    "V",
    "D",
    "N",
}

default_filter_values_mdmg = {"D_max": 0.1, "phi": 500, "q": 0.25}
default_filter_values_filterBAM = {"breadth_exp_ratio": 0.5, "coverage_mean": 0.1}
tax_ranks = [
    "superkingdom",
    "lineage",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
]
genome_selection_options = ["random", "most_abundant", "least_abundant"]
# schema_init_ag = {
#     '<abund_table>': Use(open, error='<abund_table> should be readable'),
#     '<genome_table>': Use(open, error='genome_table should be readable'),
#     '--cpus': Optional(Or(None, And(Use(int), lambda n: 0 < n < cpu_count() - 1)),
#                        error='-n=N should be integer between 0 < N < {}'.format(cpu_count() - 1)),
#     '--prop': Optional(Or(None, And(Use(float), lambda n: 0 < n < 1),
#                           error='Proportion (-p) should be a float 0 < p < 1')),
#     '--min-length': Optional(Or(None, And(Use(int), lambda min_len: min_len > 0),
#                                 error='Minimum fragment length should be greater than 0')),
#     '--read-len': Optional(Or(None, And(Use(int), lambda read_len: read_len > 0),
#                               error='Read lenhgth should be greater than 0')),
#     '--num-reads': Optional(Or(None, And(Use(float), lambda min_len: min_len > 0),
#                                error='The number of reads should be longer than 0')),
#     '--length-dist': Optional(Or(None, And(Use(str), lambda dist: dist == "lognormal"),
#                                  error='We only accept lognormal at the moment')),
#     '--seqSys': Optional(Or(None, And(Use(str), lambda seqSys: seqSys in seqSys_props),
#                             error='The sequencing technology is not supported by ART')),
#     '--mode-len-ancient': Optional(Or(None, And(Use(int), lambda m: m > 0),
#                                       error="Mode fragment length should be greater than 0")),
#     '--mode-len-modern': Optional(Or(None, And(Use(int), lambda m: m > 0),
#                                      error="Mode fragment length should be greater than 0")),
#     '--coverage': Optional(Or(None, And(Use(str), lambda m: f.str2dict(m)),
#                               error='Coverage values must be integers or float')),
#     object: object
# }

schema_init_ag = {
    "<config>": Use(open, error="config should be readable"),
    object: object,
}

schema_init_w = {
    "<config>": Use(open, error="config should be readable"),
    object: object,
}

schema_init_ar = {
    "<config>": Use(open, error="config should be readable"),
    "<genome_table>": Use(open, error="genome_table should be readable"),
    object: object,
}


schema_art_config = {"art_config": Use(open, error="art_config should be readable")}

# schema_fragSim_config = {
#     "fragSim_config": Use(open, error='fragSim_config should be readable')
# }

# schema_deamSim_config = {
#     "deamSim_config": Use(open, error='deamSim_config should be readable')
# }


ag_schema_config = {
    "genome-table": And(Use(open, error="genome-table should be readable")),
    Optional("abund-table", default=False): Or(
        bool, And(Use(open)), error="Problem finding the abundance table file"
    ),
    Optional("genome-comp", default=False): Or(
        bool,
        And(Use(open)),
        error="Problem finding the custom genome composition table file",
    ),
    Optional("read-length-freqs", default=False): Or(
        bool,
        And(Use(open)),
        error="Problem finding the read length frequencies JSON file",
    ),
    Optional("cpus", default=1): Or(
        None,
        And(Use(int), lambda n: 0 < n < cpu_count() - 1),
        error="Number of CPUs should be integer between 0 < N < {}".format(
            cpu_count() - 1
        ),
    ),
    Optional("prop-ancient-comm", default=0.1): Or(
        None,
        And(Use(float), lambda n: 0 <= n < 1),
        error="Proportion should be a float 0 < p < 1",
    ),
    Optional("min-length", default=30): Or(
        None,
        And(Use(int), lambda min_len: min_len > 0),
        error="Minimum fragment length should be greater than 0",
    ),
    Optional("read-len", default=150): Or(
        None,
        And(Use(int), lambda read_len: read_len > 0),
        error="Read lenhgth should be greater than 0",
    ),
    Optional("num-reads", default=1e6): Or(
        None,
        And(Use(float), lambda min_len: min_len > 0),
        error="The number of reads should be larger than 0",
    ),
    Optional("length-dist", default="lognormal"): Or(
        None,
        And(Use(str), lambda dist: dist == "lognormal"),
        error="We only accept lognormal at the moment",
    ),
    Optional("seqSys", default="HS25"): Or(
        None,
        And(Use(str), lambda seqSys: seqSys in seqSys_props),
        error="The sequencing technology is not supported by ART",
    ),
    Optional("mode-len-ancient", default=40): Or(
        None,
        And(Use(int), lambda m: m > 0),
        error="Mode fragment length should be greater than 0",
    ),
    Optional("mode-len-modern", default=100): Or(
        None,
        And(Use(int), lambda m: m > 0),
        error="Mode fragment length should be greater than 0",
    ),
    Optional("coverage", default=[0.1, 5]): Optional(
        Or(
            None,
            And(list, lambda x: all(i > 0 for i in x) & (x[0] < x[1])),
            error="Coverage values must be integers or float, and min should be larger than max",
        )
    ),
    "library": And(
        str,
        lambda x: x in ["single", "paired"],
        error="Library has to be single or paired",
    ),
    Optional("enforce-cov", default=True): And(
        bool, error="The enforce coverage option should be True/False"
    ),
    Optional("single-cov", default=False): And(
        bool, error="The single option should be True/False"
    ),
    object: object,
}

w_schema_config = {
    "genome-paths": And(
        Use(lambda x: f.get_open_func(x)(x, "rt")),
        error="genome_paths should be readable",
    ),
    "filterBAM-stats": And(
        Use(lambda x: f.get_open_func(x)(x, "rt")),
        error="Problem finding the abundance table file",
    ),
    "mdmg-results": And(
        Use(lambda x: f.get_open_func(x)(x, "rt")),
        error="mdmg_results should be readable",
    ),
    Optional("filterBAM-stats-filter", default=default_filter_values_filterBAM,): And(
        Use(
            lambda x: f.check_filter_conditions(
                x,
                default_filter_values_filterBAM,
                filters=["breadth_exp_ratio", "coverage_mean"],
            )
        ),
        error="There's a problem with the filter. Please check the keys in the dict are valid and that all the values are >= 0",
    ),
    # "mdmg_misincorporation": And(
    #     Use(lambda x: f.get_open_func(x)(x, "rt")),
    #     error="mdmg_misincorporation should be readable",
    # ),
    Optional("mdmg-filter-conditions", default=default_filter_values_mdmg,): And(
        Use(
            lambda x: f.check_filter_conditions(
                x,
                default_filter_values_mdmg,
                filters=[
                    "D_max",
                    "lambda_LR",
                    "q",
                    "phi",
                    "Bayesian_D_max",
                    "Bayesian_z",
                ],
            )
        ),
        error="There's a problem with the filter. Please check the keys in the dict are valid and that all the values are >= 0",
    ),
    Optional("taxonomic-rank", default="genus"): And(
        str,
        lambda x: x in tax_ranks,
        error=f"Invalid taxonomic ranks. Pick one valid from {', '.join(str(x) for x in tax_ranks)}",
    ),
    Optional("max-genomes-nondamaged", default=10): Or(
        None,
        And(Use(int), lambda max_genomes_nondamaged: max_genomes_nondamaged > 0),
        error="The number of non-damaged genomes selected should be greater than 0",
    ),
    Optional("max-genomes-damaged", default=10): Or(
        None,
        And(Use(int), lambda max_genomes_damaged: max_genomes_damaged > 0),
        error="The number of damaged genomes selected should be greater than 0",
    ),
    Optional("max-genomes-nondamaged-selection", default="most_abundant"): And(
        str,
        lambda x: x in genome_selection_options,
        error=f"Invalid option. Pick one valid from {', '.join(str(x) for x in tax_ranks)}",
    ),
    Optional("cpus", default=1): Or(
        None,
        And(Use(int), lambda n: 0 < n < cpu_count() - 1),
        error="Number of CPUs should be integer between 0 < N < {}".format(
            cpu_count() - 1
        ),
    ),
    Optional("sample-name", default=generate_random_sample_name()): Or(
        None,
        And(
            Use(str),
            lambda x: re.match("^[A-Za-z0-9_-]*$", x),
            error="The sample name contains invalid characters. Only  alphanumeric characters and dash and underscore characters are allowed",
        ),
    ),
    object: object,
}


ar_schema_config = {
    And("global"): {
        "fragSim": And(
            str,
            lambda x: shutil.which(x) is not None,
            error="fragSim executable not found",
        ),
        "deamSim": And(
            str,
            lambda x: shutil.which(x) is not None,
            error="fragSim executable not found",
        ),
        "adptSim": And(
            str,
            lambda x: shutil.which(x) is not None,
            error="fragSim executable not found",
        ),
        "art": And(
            str,
            lambda x: shutil.which(x) is not None,
            error="art_illumina executable not found",
        ),
        "AdapterRemoval": And(
            str,
            lambda x: shutil.which(x) is not None,
            error="AdapterRemoval executable not found",
        ),
        "libprep": And(
            str,
            lambda x: x in ["single", "double"],
            error="Lib prep has to be single or double",
        ),
        # Optional('seq_depth', default=1e-5): Or(None, And(Use(float), lambda x: x > 0),
        #                                         error='Sequencing depth has to be larger than 0'),
        Optional(
            "tmp-dir",
            default=os.path.join(
                os.getcwd(), datetime.datetime.now().strftime(".tmp-%Y%m%d%H%M%S")
            ),
        ): And(Use(str)),
        Optional(
            "output-dir",
            default=os.path.join(
                os.getcwd(), datetime.datetime.now().strftime("sim_reads-%Y%m%d%H%M%S")
            ),
        ): And(Use(str)),
        Optional("cpus", default=1): Or(
            None,
            And(Use(int), lambda n: 0 <= n < cpu_count()),
            error="Number of CPUs should be integer between 0 < N < {}".format(
                cpu_count() - 1
            ),
        ),
        Optional("remove-adapters", default=True): And(
            bool, error="The key dist should be True/False"
        ),
    },
    And("fragSim"): {
        And("ancient-genomes"): And(
            str,
            lambda x: os.path.exists(x),
            error="Cannot open the file with the ancient genomes information.",
        ),
        And("modern"): {
            Optional("--case", default=False): And(
                bool, error="The key case should be True/False"
            ),
            Optional("--comp", default=False): And(
                bool, error="The key comp should be True/False"
            ),
            Optional("--dist", default=False): And(
                bool, error="The key dist should be True/False"
            ),
            Optional("--norev", default=False): And(
                bool, error="The key norev should be True/False"
            ),
        },
        And("ancient"): {
            Optional("--case", default=False): And(
                bool, error="The key case should be True/False"
            ),
            Optional("--comp", default=False): And(
                bool, error="The key comp should be True/False"
            ),
            Optional("--dist", default=False): And(
                bool, error="The key dist should be True/False"
            ),
            Optional("--norev", default=False): And(
                bool, error="The key norev should be True/False"
            ),
        },
    },
    And("deamSim"): {
        Optional("-mapdamage", default=False): Or(
            bool,
            And(str, lambda x: os.path.exists(x)),
            error="Problem finding the mapDamage misincorporation file",
        ),
        Optional("-damage", default=[0.03, 0.4, 0.01, 0.7]): Or(
            bool,
            And(list, lambda x: all(i > 0 for i in x)),
            error="The key damage must be a list with positive values",
        ),
    },
    And("adptSim"): {
        Optional("-f", default=FWD_ADPT): And(
            And(Use(str), lambda x: IUPACAmbiguousDNA.issuperset(x)),
            error="Forward adapter is not a valid sequence",
        ),
        Optional("-s", default=REV_ADPT): And(
            And(Use(str), lambda x: IUPACAmbiguousDNA.issuperset(x)),
            error="Reverse adapter is not a valid sequence",
        ),
    },
    And("art"): {
        Optional("--qprof1", default=False): Or(
            bool,
            And(str, lambda x: os.path.exists(x)),
            error="Problem finding the first-read quality profile",
        ),
        Optional("--qprof2", default=False): Or(
            bool,
            And(str, lambda x: os.path.exists(x)),
            error="Problem finding the second-read quality profile",
        ),
    },
    # And("AdapterRemoval"): {
    #     Optional("--adapter1", default=FWD_ADPT): And(
    #         And(Use(str), lambda x: IUPACAmbiguousDNA.issuperset(x)),
    #         error="Forward adapter is not a valid sequence",
    #     ),
    #     Optional("--adapter2", default=REV_ADPT): And(
    #         And(Use(str), lambda x: IUPACAmbiguousDNA.issuperset(x)),
    #         error="Reverse adapter is not a valid sequence",
    #     ),
    #     Optional("--minlength", default=30): Or(
    #         None,
    #         And(Use(int), lambda n: 0 <= n),
    #         error="Minimum read length should be larger than 0",
    #     ),
    #     Optional("--trimns", default=True): And(
    #         bool, error="The key trimns should be True/False"
    #     ),
    #     Optional("--collapse", default=True): And(
    #         bool, error="The key collapse should be True/False"
    #     ),
    #     Optional("--trimqualities", default=True): And(
    #         bool, error="The key trimqualities should be True/False"
    #     ),
    #     Optional("--minquality", default=20): Or(
    #         None,
    #         And(Use(int), lambda n: 0 <= n),
    #         error="Minimum read quality value should be larger than 0",
    #     ),
    #     Optional("--minadapteroverlap", default=30): Or(
    #         None,
    #         And(Use(int), lambda n: 0 < n),
    #         error="Minimum overlap length should be larger or equal than 0",
    #     ),
    #     Optional("--preserve5p", default=True): And(
    #         bool, error="The key preserve5p should be True/False"
    #     ),
    # },
}

seqSys_props = {
    "GA2": {"se": 50, "pe": 75},
    "HS20": {"se": 100, "pe": 100},
    "HS25": {"se": 125, "pe": 150},
    "HSXt": {"se": 150, "pe": 150},
    "MSv1": {"se": 250, "pe": 250},
    "MSv3": {"se": 250, "pe": 250},
    "Other": {"se": None, "pe": None},
}
