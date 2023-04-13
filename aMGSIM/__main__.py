#!/usr/bin/env python

# import
import sys
import logging

from docopt import docopt

# import Utils
from MGSIM.Commands import Communities

from aMGSIM.Commands.Estimate import estimate
from aMGSIM.Commands.AncientGenomes import get_ancient_genomes
from aMGSIM.Commands.AncientReads import get_ancient_reads
from aMGSIM.Commands.ProteinAnalysis import do_proteins_analysis
from aMGSIM.Commands.MultiCov import generate_multi_cov
from aMGSIM.Commands.AddDuplicates import add_duplicates
from aMGSIM.library.cli import get_arguments
import matplotlib

log = logging.getLogger(__name__)
# logging.getLogger(name="matplotlib").setLevel(logging.WARNING)
# logging.getLogger(name="scipy").setLevel(logging.WARNING)
# logging.getLogger(name="scipy.stats").setLevel(logging.WARNING)
# logging.getLogger(name="scipy.optimize").setLevel(logging.WARNING)


def main():

    if len(sys.argv) > 1 and sys.argv[1] == "communities":
        # args = sys.argv[2:]
        # print(args)
        docs = Communities.__doc__
        # print(docs)
        args = docopt(docs, argv=sys.argv[2:], options_first=True)
        Communities.main(args)
    else:
        args = get_arguments()
        logging.basicConfig(
            level=logging.DEBUG if args.debug else logging.INFO,
            format="%(levelname)s ::: %(asctime)s ::: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True,
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
        elif args.action == "add-duplicates":
            add_duplicates(args)


if __name__ == "__main__":
    main()
