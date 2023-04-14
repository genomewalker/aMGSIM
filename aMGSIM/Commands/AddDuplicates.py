from tqdm import tqdm
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
from itertools import zip_longest
from mimetypes import guess_type
from functools import partial
import io
from decoratio import (
    load_cz_d_from_tsv,
    calc_duprate,
    calc_lander_waterman_depth_complexity_ratio,
    AmplModel,
    optimize_ratio,
    LOG_CLASSES_BASE,
    optimize_ratio_noise,
    cz_d_to_reads,
    generate_cz_distrib,
    InheffAmplModel,
)
import pandas as pd
import warnings
import math
import numpy as np
import _md5
import logging
import sys


warnings.filterwarnings("ignore", category=RuntimeWarning)


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
        print("{}: {}".format(exception_type.__name__, exception))
        print("\n*** Error: Please use --debug to see full traceback.")


def round_base_10(x):
    if x < 0:
        return 0
    elif x == 0:
        return 10
    return 10 ** math.ceil(math.log10(x))


def closest_multiple(x, multiple):
    if x == 0:
        return 0
    if multiple > x:
        return multiple
    return math.trunc((x / multiple)) * multiple


def get_nreads(n, prop, nreads, dict_):
    n1 = dict_[dict_["counts"] == n].explode(["record_id", "qual", "strand"]).shape[0]
    if prop == 0:
        return 0
    else:
        if n == 1:
            n1 = 0
        val = int((n * prop * nreads)) + (n1 - 1)
        return closest_multiple(val, n)


# from https://stackoverflow.com/questions/53751050/python-multiprocessing-understanding-logic-behind-chunksize/54032744#54032744
def calc_chunksize(n_workers, len_iterable, factor=4):
    """Calculate chunksize argument for Pool-methods.
    Resembles source-code within `multiprocessing.pool.Pool._map_async`.
    """
    chunksize, extra = divmod(len_iterable, n_workers * factor)
    if extra:
        chunksize += 1
    return chunksize


# From https://stackoverflow.com/a/52913128
def get_random_key(keys, size):
    L = len(keys)
    i = np.random.randint(0, L, size=size, dtype=np.int32)
    return [keys[i] for i in i]


tab = str.maketrans("ACTG", "TGAC")


def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]


# from https://groverj3.github.io/articles/2019-08-22_just-write-your-own-python-parsers-for-fastq-files.html
def parse_zip_longest(input_fastq):
    encoding = guess_type(input_fastq)[1]  # uses file extension
    _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open
    with _open(input_fastq) as input_handle:
        fastq_iterator = (line.rstrip() for line in input_handle)
        for record in zip_longest(*[fastq_iterator] * 4):
            yield record


def md5hasher(what_text):
    return _md5.md5(what_text.encode("utf")).hexdigest()


# Define function to read FASTQ file and update shared dictionary with sequences
def update_dict(dict_, counts, filename, progress_bar):
    data = parse_zip_longest(filename)
    # p0 = re.compile(r"@\S+___(\S+)__.*")

    for record in data:
        # create a md5 hash of the sequence
        # seq_hash_fwd = md5hasher(record[1])
        seq = record[1]
        # check if the sequence is already in the dictionary in forward or reverse complement orientation
        if seq in dict_:
            strand = "fwd"
        else:
            seq_rev = reverse_complement_table(record[1])
            if seq_rev in dict_:
                seq = seq_rev
                strand = "rev"
            else:
                strand = "fwd"

        if seq in dict_:
            dict_[seq]["record_id"].append(record[0])
            dict_[seq]["strand"].append(strand)
            dict_[seq]["qual"].append(record[3])
            dict_[seq]["counts"] += 1
        else:
            dict_[seq]["record_id"] = [record[0]]
            dict_[seq]["strand"] = [strand]
            dict_[seq]["qual"] = [record[3]]
            dict_[seq]["counts"] = 1

        progress_bar.update(1)


def get_pred_df(obs_cz_d, pred_cz_d, n_clones):
    obs_cz_d_reads = cz_d_to_reads(obs_cz_d)
    pred_cz_d_reads = cz_d_to_reads(pred_cz_d)
    n_reads = n_clones / (1.0 - calc_duprate(obs_cz_d))

    res = []
    for z in range(len(obs_cz_d)):
        if n_clones >= 2:
            res.append(
                {
                    "clone_size": z,
                    "n_reads": round(obs_cz_d_reads[z] * n_reads),
                    "n_clones": round(obs_cz_d[z] * n_clones),
                    "pred_n_reads": pred_cz_d_reads[z]
                    if len(pred_cz_d_reads) > z
                    else 0.0,
                    "pred_n_clones": pred_cz_d[z] if len(pred_cz_d) > z else 0.0,
                }
            )

        else:
            res.append(
                {
                    "clone_size": z,
                    "n_reads": obs_cz_d_reads[z],
                    "n_clones": obs_cz_d[z],
                    "pred_n_reads": pred_cz_d_reads[z]
                    if len(pred_cz_d_reads) > z
                    else 0.0,
                    "pred_n_clones": pred_cz_d[z] if len(pred_cz_d) > z else 0.0,
                }
            )
    df = pd.DataFrame(res)
    return df


def get_duplicates(
    seq_dict, clone_size, n_reads_dup_aprox, seq_dict_done, progress_bar
):
    # get keys from seq_dict not in seq_dict_done
    # we get the sequences with the same clone size
    df_cs = seq_dict[seq_dict["counts"] == clone_size]
    n = df_cs.shape[0]
    if df_cs.shape[0] > 0:
        df_cs = df_cs.explode(["record_id", "qual", "strand"])

    # we get the sequences with only one read
    df_uniq = seq_dict[seq_dict["counts"] == 1]

    to_extract = int(((n_reads_dup_aprox) / clone_size) - n)

    if to_extract < 0:
        to_extract = 0

    seq_dict_filt = df_uniq.sample(n=to_extract)
    seq_dict_filt = seq_dict_filt.explode(["record_id", "qual", "strand"])

    # convert
    seq_dict_filt_dups = seq_dict_filt.loc[seq_dict_filt.index.repeat(clone_size)]
    seq_dict_filt_dups = pd.concat([seq_dict_filt_dups, df_cs])
    seq_dict_filt_dups.reset_index(inplace=True, drop=True)

    # Check if string has __seq- or __dseq- and replace it

    seq_dict_filt_dups["record_id"].replace(
        "__seq-", "__dseq-", regex=True, inplace=True
    )
    seq_dict_filt_dups.loc[
        seq_dict_filt_dups.groupby("record_id")["record_id"].head(1).index,
        "record_id",
    ] = seq_dict_filt_dups.loc[
        seq_dict_filt_dups.groupby("record_id")["record_id"].head(1).index,
        "record_id",
    ].str.replace(
        "__dseq-", "__seq-"
    )

    seq_dict_filt_dups = pd.concat([seq_dict_filt_dups, df_cs])
    seq_dict = seq_dict.drop(seq_dict_filt.index)
    seq_dict_done.append(seq_dict_filt_dups)
    progress_bar.update(1)


log = logging.getLogger(__name__)
# logging.getLogger(name="matplotlib").setLevel(logging.WARNING)
# logging.getLogger(name="scipy").setLevel(logging.WARNING)
# logging.getLogger(name="scipy.stats").setLevel(logging.WARNING)
# logging.getLogger(name="scipy.optimize").setLevel(logging.WARNING)


# Define main function to read both FASTQ files in parallel and compute LCS
def add_duplicates(args):
    # Parse command-line arguments
    global debug
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    # simulating reads
    if args.debug:
        debug = True
    else:
        debug = False

    try:
        args.model = AmplModel.factory(args.model)
    except ValueError:
        sys.exit("ERROR: failed to parse model '{}'".format(args.model))
    if isinstance(args.model, InheffAmplModel):
        assert hasattr(args.model, "n_efficiency_bins")
        assert hasattr(args.model, "n_sims_per_bin")
        assert hasattr(args.model, "legacy_sims")
        if args.inheff_n_efficiency_bins:
            args.model.n_efficiency_bins = args.inheff_n_efficiency_bins
        if args.inheff_n_sims_per_bin:
            args.model.n_sims_per_bin = args.inheff_n_sims_per_bin
        if args.inheff_legacy_sims:
            args.model.legacy_sims = True
        if args.inheff_uncouple_beta_parameters:
            args.model.uncouple_beta_parameters = True

    sys.excepthook = exceptionHandler
    # Initialize shared dictionary to store sequences
    seq_dict = defaultdict(dict)
    counts = defaultdict(int)
    log.info("Reading FASTQ files...")
    # Initialize progress bars
    progress_desc1 = f"Reading {args.fastq}"
    progress_bar1 = tqdm(total=0, desc=progress_desc1, unit=" reads", leave=False)
    # Read both FASTQ files in parallel with progress bars
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for filename, progress_bar in zip([args.fastq], [progress_bar1]):
            futures.append(
                executor.submit(update_dict, seq_dict, counts, filename, progress_bar)
            )
        for future in as_completed(futures):
            pass
    # Close progress bars
    progress_bar1.close()

    log.info("Converting sequence dictionary to dataframe...")
    seq_dict = pd.DataFrame.from_dict(seq_dict, orient="index")
    seq_dict.reset_index(inplace=True)
    seq_dict.rename(columns={"index": "seq"}, inplace=True)

    len_all = seq_dict["counts"].sum()
    # calculate the duplication rate of the seq_dict

    duplicated_reads = (
        seq_dict[seq_dict["counts"] > 1]
        .groupby(["counts"])["counts"]
        .size()
        .reset_index(name="n")
    )
    duplicated_reads["duplicated_reads"] = (
        duplicated_reads["n"] * duplicated_reads["counts"]
    )
    # print(f"::: {reads_duplicated:,} reads are duplicated")
    dup_reads = duplicated_reads["n"].values.sum()
    duprate = dup_reads / len_all
    # dups = seq_dict.groupby(["counts"])["counts"].size().values
    # dups = np.insert(dups, 0, 0)
    # print(dups / len_all)
    # dups = dups / len_all
    # print(abs(sum(dups) - 1.0))
    # duprate = calc_duprate(dups)
    log.info(f"Total number of reads: {len_all:,} ({duprate:.1%} duplication rate)")

    log.info("Calculating duplication estimates...")
    with open(args.clone_size_freqs, "r") as f:
        obs_cz_d, n_clones = load_cz_d_from_tsv(f)

    obs_duprate = calc_duprate(obs_cz_d)

    lw_cr = calc_lander_waterman_depth_complexity_ratio(obs_duprate)
    log.info(f"Original data shows {obs_duprate:.1%} of duplicates")
    log.info(
        f"Lander-Waterman/noiseless amplification depth-complexity ratio estimate: {lw_cr:.3g}"
    )

    # generate the number of reads needed

    if args.ratio_max is None:
        args.ratio_max = 1 / (1 - obs_duprate)
    if not args.model.needs_optimization():
        model = args.model
        log.info("Generating the amplification factor distribution...")
        ampl_cz_d = model.get_ampl_cz_d(LOG_CLASSES_BASE)
        log.info(
            "::: The provided model had a noisiness of {:.3f}; a skewness of {:.2f}".format(
                AmplModel.calc_noisiness(ampl_cz_d, LOG_CLASSES_BASE),
                AmplModel.calc_skewness(ampl_cz_d, LOG_CLASSES_BASE),
            )
        )
        log.info("Optimizing the depth-complexity ratio...")
        ratio = optimize_ratio(obs_cz_d, ampl_cz_d, args)
    else:
        log.info(
            "Optimizing the depth-complexity ratio & amplification noise...",
        )
        ratio, model = optimize_ratio_noise(obs_cz_d, args)
        ampl_cz_d = model.get_ampl_cz_d(LOG_CLASSES_BASE)

    pred_cz_d = generate_cz_distrib(ratio, ampl_cz_d, LOG_CLASSES_BASE)
    df = get_pred_df(obs_cz_d, pred_cz_d, n_clones)

    df["n_clones_dup_aprox"] = df.apply(
        lambda x: get_nreads(
            n=x["clone_size"], prop=x["pred_n_clones"], nreads=len_all, dict_=seq_dict
        ),
        axis=1,
    )
    df_filt = df[(df["n_clones_dup_aprox"] > 0) & (df["clone_size"] > 0)]
    # we keep the one with n_clones_dup > 0
    df.to_csv("test.tsv", sep="\t", index=False)
    df_filt.to_csv("test1.tsv", sep="\t", index=False)

    seq_dict_done = []

    log.info("Generating duplicates...")
    for clone_size, n_reads_dup_aprox in tqdm(
        zip(df_filt["clone_size"], df_filt["n_clones_dup_aprox"]),
        total=df_filt.shape[0],
        desc="Processing duplicates",
        unit=" sizes",
        disable=args.quiet,
        leave=False,
    ):
        # get keys from seq_dict not in seq_dict_done
        # we get the sequences with the same clone size
        df_cs = seq_dict[seq_dict["counts"] == clone_size]
        n = df_cs.shape[0]
        if df_cs.shape[0] > 0:
            df_cs = df_cs.explode(["record_id", "qual", "strand"])

        # we get the sequences with only one read
        df_uniq = seq_dict[seq_dict["counts"] == 1]

        to_extract = int(((n_reads_dup_aprox) / clone_size) - n)

        if to_extract < 0:
            to_extract = 0
        seq_dict_filt = df_uniq.sample(n=to_extract)
        seq_dict_filt = seq_dict_filt.explode(["record_id", "qual", "strand"])
        # convert
        seq_dict_filt_dups = seq_dict_filt.loc[seq_dict_filt.index.repeat(clone_size)]
        seq_dict_filt_dups = pd.concat([seq_dict_filt_dups, df_cs])
        seq_dict_filt_dups.reset_index(inplace=True, drop=True)

        # Check if string has __seq- or __dseq- and replace it

        seq_dict_filt_dups["record_id"].replace(
            "__seq-", "__dseq-", regex=True, inplace=True
        )
        seq_dict_filt_dups.loc[
            seq_dict_filt_dups.groupby("record_id")["record_id"].head(1).index,
            "record_id",
        ] = seq_dict_filt_dups.loc[
            seq_dict_filt_dups.groupby("record_id")["record_id"].head(1).index,
            "record_id",
        ].str.replace(
            "__dseq-", "__seq-"
        )

        seq_dict_filt_dups = pd.concat([seq_dict_filt_dups, df_cs])
        seq_dict = seq_dict.drop(seq_dict_filt.index)
        seq_dict_done.append(seq_dict_filt_dups)

    seq_dict_done.append(seq_dict.explode(["record_id", "qual", "strand"]))
    seq_dict_done = pd.concat(seq_dict_done)
    # seq_dict_done = pd.concat(seq_dict_done, seq_dict_sngl)
    log.info(
        f"Total number of reads after adding duplicates: {seq_dict_done.shape[0]:,}"
    )

    # write to fastq file
    with gzip.GzipFile("test.fq.gz", mode="wb", compresslevel=6) as gz_file:
        with io.BufferedWriter(gz_file, buffer_size=100 * (2**20)) as buffer:
            for record_id, seq, qual, strand in tqdm(
                zip(
                    seq_dict_done["record_id"],
                    seq_dict_done["seq"],
                    seq_dict_done["qual"],
                    seq_dict_done["strand"],
                ),
                total=seq_dict_done.shape[0],
                desc="Writing to file",
                unit=" reads",
                disable=args.quiet,
                leave=False,
            ):
                if strand == "rev":
                    seq = reverse_complement_table(seq)
                buffer.write(f"{record_id}\n{seq}\n+\n{qual}\n".encode("utf-8"))
