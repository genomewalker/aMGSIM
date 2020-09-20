# Constants
from schema import And, Optional, Or, Use
from multiprocessing import cpu_count
import aMGSIM.library.functions as f
import os
import shutil
import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
            And(Use(int), lambda n: 0 < n < cpu_count() - 1),
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
}

seq = Seq("ATGCGTGCAT")
record = SeqRecord(seq, id="test", annotations={"molecule_type": "DNA"})

default_config = {
    "art_exe": "art_illumina",
    # Path to the fragSim binary
    "fragSim": "fragSim",
    # Path to the deamSim binary
    "deamSim": "deamSim",
    # Path to the adptSim binary
    "adptSim": "adptSim",
    # Sequencing depth
    "seq_depth": 1e5,
    # Temporary folder
    "tmp_dir": ".sim_reads",
    # Number of cpus
    "cpus": 1,
    # Output dir
    "output_dir": "test_out",
}

default_art_config = {
    "paired": True,
    "read_length": 150,
    "amplicon": True,
    "rndSeed": 12345,
    "seqSys": "HS25",
}


default_fragSim_config = {
    "comp": False,  # [file] Base composition for the fragments (default none)
    # [file] Distance from ends to consider  (default: 1) if this is not specified, the base composition will only reflect the chromosome file used
    "dist": False,
    # Do not reverse complement (default: rev. comp half of seqs.)
    "norev": False,
    # Do not set the sequence to upper-case (default: uppercase the seqs.)
    "case": False,
    # We define for each sample max/min
    # [file] Open file with min/max read lengths in the following format: sample[TAB]min_length[TAB]max_length
    "fragment_length_file": False,
    # [length] Generate fragments of fixed length for all samples (default: 20)
    "default_fragment_length": 60,
    # [file] Open file with lognormal location/scale following format: sample[TAB]length[TAB]freq
    "frag_distribution_file": False,
    # [file] Open file with size frequency in the following format: sample[TAB]lognorm_location[TAB]lognorm_scale
    "frag_lognormal_distribution_file": False,
    # [float] Location for lognormal distribution (default none)
    "default_loc": False,
    # [float] Scale for lognormal distribution      (default none)
    "default_scale": False,
    "gc": False,  # [float] Use GC bias factor  (default: 0)
    # [file] Open file with size distribution (one fragment length per line)
    "s_ancient": False,
}

default_deamSim_config = {"mapdamage": False, "damage": False}

seqSys_props = {
    "GA2": {"se": 50, "pe": 75},
    "HS20": {"se": 100, "pe": 100},
    "HS25": {"se": 125, "pe": 150},
    "HSXt": {"se": 150, "pe": 150},
    "MSv1": {"se": 250, "pe": 250},
    "MSv3": {"se": 250, "pe": 250},
    "Other": {"se": None, "pe": None},
}
