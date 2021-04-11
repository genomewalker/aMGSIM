# Constants
from schema import And, Optional, Or, Use
from multiprocessing import cpu_count
import os
import shutil
import datetime


# schema_init_ag = {
#     '<abund_table>': [Use(open, error='abund_table should be readable')],
#     '-p': Or(None, And(Use(float), lambda n: 0 < n < 1),
#              error='Proportion (-p) should be a float 0 < p < 1'),
#     object: object
# }

debug = None

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
    "genome_table": And(Use(open, error="genome_table should be readable")),
    Optional("abund_table", default=False): Or(
        bool, And(Use(open)), error="Problem finding the abundance table file"
    ),
    Optional("genome_comp", default=False): Or(
        bool,
        And(Use(open)),
        error="Problem finding the custom genome composition table file",
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
            "tmp_dir",
            default=os.path.join(
                os.getcwd(), datetime.datetime.now().strftime(".tmp-%Y%m%d%H%M%S")
            ),
        ): And(Use(str)),
        Optional(
            "output_dir",
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
    },
    And("fragSim"): {
        And("ancient_genomes"): And(
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
    And("AdapterRemoval"): {
        Optional("--adapter1", default=FWD_ADPT): And(
            And(Use(str), lambda x: IUPACAmbiguousDNA.issuperset(x)),
            error="Forward adapter is not a valid sequence",
        ),
        Optional("--adapter2", default=REV_ADPT): And(
            And(Use(str), lambda x: IUPACAmbiguousDNA.issuperset(x)),
            error="Reverse adapter is not a valid sequence",
        ),
        Optional("--minlength", default=30): Or(
            None,
            And(Use(int), lambda n: 0 <= n),
            error="Minimum read length should be larger than 0",
        ),
        Optional("--trimns", default=True): And(
            bool, error="The key trimns should be True/False"
        ),
        Optional("--collapse", default=True): And(
            bool, error="The key collapse should be True/False"
        ),
        Optional("--trimqualities", default=True): And(
            bool, error="The key trimqualities should be True/False"
        ),
        Optional("--minquality", default=20): Or(
            None,
            And(Use(int), lambda n: 0 <= n),
            error="Minimum read quality value should be larger than 0",
        ),
        Optional("--minadapteroverlap", default=30): Or(
            None,
            And(Use(int), lambda n: 0 < n),
            error="Minimum overlap length should be larger or equal than 0",
        ),
        Optional("--preserve5p", default=True): And(
            bool, error="The key preserve5p should be True/False"
        ),
    },
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
