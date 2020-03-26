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


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="circularity")
    parser.add_argument(
        "input",
        type=str,
        help="Input viral sequences in FASTA format",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Output directory"
    )
    parser.add_argument(
        "--min_len",
        type=int,
        default=2000,
        metavar="INT",
        help="Min contig length",
    )
    parser.add_argument(
        "--min_dtr",
        type=int,
        default=20,
        metavar="INT",
        help="Min length of direct terminal repeats (DTRs)",
    )
    parser.add_argument(
        "--max_lc",
        type=float,
        default=20.0,
        metavar="FLOAT",
        help="Max %% of direct terminal repeats (DTRs) classified as low complexity",
    )
    parser.add_argument(
        "--keep_rep",
        action="store_true",
        default=False,
        help="Keep direct terminal repeats (DTRs) that occur more than 2x per sequence",
    )


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


def main(args):
    utility.check_executables(["dustmasker"])
    args["tmp"] = os.path.join(args["output"], "tmp")
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    print("read input sequences...")
    genomes = {}
    for index, r in enumerate(Bio.SeqIO.parse(args["input"], "fasta")):
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
    out = open(args["output"] + "/circularity.tsv", "w")
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
