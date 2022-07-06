#!/usr/bin/env python

# import
## batteries
import sys
import logging
import time

## 3rd party
from docopt import docopt

# import Utils
## application
from MGSIM.Commands import Communities

from aMGSIM.Commands.Estimate import estimate
from aMGSIM.Commands.AncientGenomes import get_ancient_genomes
from aMGSIM.Commands.AncientReads import get_ancient_reads
from aMGSIM.Commands.ProteinAnalysis import do_proteins_analysis
from aMGSIM.Commands.MultiCov import generate_multi_cov

from aMGSIM import __version__
from aMGSIM.library.cli import get_arguments


def main():

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    if len(sys.argv) > 1 and sys.argv[1] == "communities":
        # args = sys.argv[2:]
        # print(args)
        docs = Communities.__doc__
        # print(docs)
        args = docopt(docs, argv=sys.argv[2:], options_first=True)
        Communities.main(args)
    else:
        args = get_arguments()
        logging.getLogger("my_logger").setLevel(
            logging.DEBUG if args.debug else logging.INFO
        )
        if args.action == "estimate":
            estimate(args)
        elif args.action == "multicov":
            generate_multi_cov(args)
        elif args.action == "ancient-genomes":
            get_ancient_genomes(args)
        elif args.action == "ancient-reads":
            get_ancient_reads(args)
        elif args.action == "protein-analysis":
            do_proteins_analysis(args)


if __name__ == "__main__":
    main()
