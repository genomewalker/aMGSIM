#!/usr/bin/env python

# import
## batteries
import os
import sys

## 3rd party
from docopt import docopt

# import Utils
## application
from MGSIM.Commands import Communities
from aMGSIM.Commands import AncientGenomes
from aMGSIM.Commands import ReadsAncient
from aMGSIM.Commands import ProteinAnalysis


def main(args=None):
    """Main entry point for application"""
    if args is None:
        args = sys.argv[1:]

    docs = """
aMGSIM: simulate ancient metagenomes for multiple synthetic communities

Usage:
  aMGSIM <command> [<args>...]
  aMGSIM -l | --list
  aMGSIM -h | --help
  aMGSIM --version

Options:
  -l --list     List subcommands.
  -h --help     Show this screen.
  --version     Show version.

Commands:
  Use the `list` option.
Description:
  Simulate metagenomes for multiple synthetic communities.
  See the sub-command documentation for more information on features.
    """
    # arg parse
    args = docopt(docs, version="0.1", options_first=True)

    # dict of all subcommands
    cmds = {
        "communities": Communities,
        "ancient-genomes": AncientGenomes,
        "ancient-reads": ReadsAncient,
        "protein-analysis": ProteinAnalysis,
    }

    # list subcommands
    if args["--list"]:
        cmd_list = " | ".join(cmds.keys())
        print("Available Commands:")
        print(cmd_list)
        exit()

    # running subcommand
    try:
        func = cmds[args["<command>"]]
    except KeyError:
        msg = 'ERROR: command "{}" does not exist'
        print(msg.format(args["<command>"]))
        exit()
    func.opt_parse(args["<args>"])


if __name__ == "__main__":
    main()