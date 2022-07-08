from operator import ge
from aMGSIM.library import cli as c
import logging
import pandas as pd
import sys
from pathlib import Path
import natsort as ns


def exceptionHandler(
    exception_type, exception, traceback, debug_hook=sys.__excepthook__
):
    """Print user friendly error messages normally, full traceback if DEBUG on.
    Adapted from http://stackoverflow.com/questions/27674602/hide-traceback-unless-a-debug-flag-is-set
    """
    debug = c.is_debug()
    if debug:
        print("\n*** Error:")
        # raise
        debug_hook(exception_type, exception, traceback)
    else:
        print("{}: {}".format(exception_type.__name__, exception))
        log.error("Please use --debug to see full traceback.")


log = logging.getLogger("my_logger")


def generate_multi_cov(args):
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    log.info("Loading genome composition table...")
    genome_comp = pd.read_csv(args.genome_compositions, sep="\t")
    genome_comp["reference"] = genome_comp["Taxon"].str.split("----", expand=True)[1]
    file_name = Path(args.genome_compositions).stem
    file_name = Path(f"{file_name}.{args.cov_suffix}.tsv")
    # combine columns Community and reference
    genome_comp["Community"] = (
        genome_comp["Community"] + "___" + genome_comp["reference"]
    )
    genome_comp.drop(["reference"], axis=1, inplace=True)

    genome_comp_tmp = genome_comp.copy(deep=True)

    # split column Taxon by delimiter "----"

    # drop columns Taxon and reference

    dfs = []
    # loop over each row and add teh coverage values from coverages
    for index, row in genome_comp_tmp.iterrows():
        row_tmp = row.copy(deep=True)
        for coverage in args.coverages:
            row_tmp["Community"] = row_tmp["Community"] + "___" + str(coverage)
            row_tmp["Coverage"] = coverage
            dfs.append(row_tmp.to_dict())
            row_tmp = row.copy(deep=True)

    # concatenate all the dataframes
    genome_comp["Community"] = genome_comp["Community"] + "___raw"
    genome_comp = pd.concat([genome_comp, pd.DataFrame(dfs)])
    log.info(f"Saving new genome composition table to {file_name}")
    genome_comp.sort_values(["Community", "Taxon"]).to_csv(
        file_name, sep="\t", index=False
    )
