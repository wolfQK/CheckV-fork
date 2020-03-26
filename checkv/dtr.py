#!/usr/bin/env python

import argparse
import os
import subprocess as sp
import time

import Bio.SeqIO
from checkv import utility


class Genome:
    def __init__(self):
        self.seq = None
        self.id = None
        self.num = None
        self.length = None
        self.dtr = None


class DTR:
    def __init__(self):
        pass


class Gene:
    def __init__(self):
        self.length = None
        self.reads = 0
        self.bases = 0
        self.depth = 0.0
        self.intra_freq = None
        self.inter_freq = None
        self.alns = []


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter, usage=argparse.SUPPRESS
    )
    parser.add_argument("program", help=argparse.SUPPRESS)
    parser.add_argument(
        "-i",
        dest="in",
        type=str,
        required=True,
        metavar="PATH",
        help="""Input viral sequences in FASTA format""",
    )
    parser.add_argument(
        "-o", dest="out", type=str, required=True, metavar="PATH", help="""Output path"""
    )
    parser.add_argument(
        "--min_len",
        type=int,
        default=2000,
        metavar="INT",
        help="""Min contig length (2000)""",
    )
    parser.add_argument(
        "--min_dtr",
        type=int,
        default=20,
        metavar="INT",
        help="""Min length of DTRs (20)""",
    )
    parser.add_argument(
        "--max_lc",
        type=float,
        default=20.0,
        metavar="FLOAT",
        help="""Max %% of DTR classified as low complexity (20.0)""",
    )
    parser.add_argument(
        "--keep_rep",
        action="store_true",
        default=False,
        help="""Keep DTRs that occur more than 2x per sequence""",
    )
    return vars(parser.parse_args())


def fetch_dtr(seq, min_length=20):

    # see if minimal substring occurs in 2nd half of seq
    substring = seq[0:min_length]
    pos = seq.rfind(substring)
    if pos < len(seq) / 2:
        return ""

    # see if substring terminal repeat
    substring = seq[pos:]
    if seq[0 : len(substring)] == substring:
        return substring

    return ""


def main():

    program_start = time.time()
    args = parse_arguments()
    args["tmp"] = os.path.join(args["out"], "tmp")
    if not os.path.exists(args["out"]):
        os.makedirs(args["out"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    print("read input sequences...")
    genomes = {}
    for index, r in enumerate(Bio.SeqIO.parse(args["in"], "fasta")):
        genome = Genome()
        genome.num = str(index)
        genome.id = r.id
        genome.seq = str(r.seq).upper()
        genome.length = len(genome.seq)
        genomes[index] = genome

    print("find direct terminal repeats...")
    for genome in genomes.values():
        dtr = DTR()
        dtr.seq = fetch_dtr(genome.seq, args["min_dtr"])
        dtr.length = len(dtr.seq)
        dtr.count = genome.seq.count(dtr.seq)
        dtr.dust = 0
        genome.dtr = dtr

    print("check sequence complexity...")
    args["dtr_path"] = f"{args['tmp']}/dtr.fna"
    args["dustmaskerout"] = os.path.join(args["tmp"], "dustmasker.txt")
    with open(args["dtr_path"], "w") as f:
        for num in sorted(genomes.keys()):
            genome = genomes[num]
            if genome.dtr.length >= 10:
                f.write(">" + genome.num + "\n" + genome.dtr.seq + "\n")
    utility.run_dustmasker(args["dtr_path"], args["dustmaskerout"])
    for line in open(f"{args['tmp']}/dustmasker.txt"):
        num, start, end = line.split()
        num = int(num[1:])
        genome = genomes[num]
        genome.dtr.dust = int(end) - int(start) + 1
    print("write summary...")
    header = [
        "genome_num",
        "genome_id",
        "genome_length",
        "dtr_length",
        "dtr_count",
        "dtr_dust",
        "is_complete",
    ]
    out = open(args["out"] + "/circularity.tsv", "w")
    out.write("\t".join(header) + "\n")
    for genome_num in sorted(genomes.keys()):
        genome = genomes[genome_num]
        dtr = genome.dtr
        is_complete = "No"
        if (
            genome.length >= args["min_len"]
            and dtr.length >= args["min_dtr"]
            and (dtr.count == 2 or args["keep_rep"])
            and (100.0 * dtr.dust / dtr.length <= args["max_lc"])
        ):
            is_complete = "Yes"
        row = [
            genome_num,
            genome.id,
            genome.length,
            dtr.length,
            dtr.count,
            dtr.dust,
            is_complete,
        ]
        out.write("\t".join([str(_) for _ in row]) + "\n")

    print("done!")
