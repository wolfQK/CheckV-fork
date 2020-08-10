import argparse
import bisect
import csv
import logging
import operator
import os
import shutil
import subprocess as sp
import sys
import time
import numpy as np
import checkv
from checkv import utility


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="end_to_end")
    parser.add_argument(
        "input", type=str, help="Input nucleotide sequences in FASTA format"
    )
    parser.add_argument("output", type=str, help="Output directory")
    parser.add_argument(
        "-d",
        dest="db",
        type=str,
        required=False,
        metavar="PATH",
        help="Reference database path. By default the CHECKVDB environment variable is used",
    )
    parser.add_argument(
        "-t",
        dest="threads",
        type=int,
        default=1,
        metavar="INT",
        help="Number of threads to use for Prodigal and DIAMOND",
    )
    parser.add_argument(
        "--restart",
        action="store_true",
        default=False,
        help="Overwrite existing intermediate files. By default CheckV continues where program left off",
    )
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Suppress logging messages",
    )


def main(args):

    # if restart flag, delete tmp only 1x
    args["tmp"] = os.path.join(args["output"], "tmp")
    if args["restart"] and os.path.exists(args["tmp"]):
        shutil.rmtree(args["tmp"])
    args["restart"] = False

    from checkv import contamination

    contamination.main(args)

    from checkv import completeness

    completeness.main(args)

    from checkv import complete_genomes

    complete_genomes.main(args)

    from checkv import quality_summary

    quality_summary.main(args)
