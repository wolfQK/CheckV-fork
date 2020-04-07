#!/usr/bin/env python

import argparse
import logging
import os
import time
import csv
import sys
import string
import Bio.SeqIO
from checkv import utility


class Genome:
    def __init__(self):
        self.seq = None
        self.id = None
        self.length = None
        self.dtr = None


class TR:
    def __init__(self):
        pass


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="terminal_repeats")
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
        "--min_contig_len",
        type=int,
        default=2000,
        metavar="INT",
        help="Min contig length",
    )
    parser.add_argument(
        "--min_tr_len",
        type=int,
        default=20,
        metavar="INT",
        help="Min length of TR",
    )
    parser.add_argument(
        "--max_tr_count",
        type=int,
        default=5,
        metavar="INT",
        help="Max occurrences of TR per contig",
    )
    parser.add_argument(
        "--max_tr_dust",
        type=float,
        default=20.0,
        metavar="FLOAT",
        help="Longest low complexity region per TR, as %% of TR length",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Suppress logging messages",
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


def reverse_complement(seq):
    if sys.version_info > (3,0):
        trans = str.maketrans("ACTG", "TGAC")
    else:
        trans = string.maketrans("ACTG", "TGAC")
    return seq[::-1].translate(trans)


def fetch_itr(seq, min_len=20, max_len=1000):
    rev = reverse_complement(seq)
    # see if minimal substring occurs at end
    if seq[:min_len] == rev[:min_len]:
        # extend to maximum substring, up to <max_len>
        i = min_len+1
        while seq[:i] == rev[:i] and i <= max_len:
            i += 1
        return seq[:i - 1]
    # no match
    else:
        return ""

def main(args):

    program_start = time.time()
    logger = utility.get_logger(args["quiet"])
    utility.check_executables(["dustmasker"])
    args["tmp"] = os.path.join(args["output"], "tmp")
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    logger.info("[1/5] Reading input sequences...")
    genomes = {}
    for r in Bio.SeqIO.parse(args["input"], "fasta"):
        genome = Genome()
        genome.id = r.id
        genome.seq = str(r.seq).upper()
        genome.length = len(genome.seq)
        genomes[genome.id] = genome

    logger.info("[2/5] Finding direct/inverted terminal repeats...")
    for genome in genomes.values():
        dtr = TR()
        dtr.seq = fetch_dtr(genome.seq, args["min_tr_len"])
        dtr.length = len(dtr.seq)
        dtr.count = genome.seq.count(dtr.seq)
        dtr.dust = 0
        genome.dtr = dtr
        itr = TR()
        itr.seq = fetch_itr(genome.seq, args["min_tr_len"])
        itr.length = len(itr.seq)
        itr.count = genome.seq.count(itr.seq)
        itr.dust = 0
        genome.itr = itr

    logger.info("[3/5] Writing terminal repeat sequences...")
    args["tr_path"] = f"{args['tmp']}/tr.fna"
    with open(args["tr_path"], "w") as f:
        for id in sorted(genomes.keys()):
            genome = genomes[id]
            if genome.dtr.length >= args["min_tr_len"]:
                f.write(">DTR_" + genome.id + "\n" + genome.dtr.seq + "\n")
            if genome.itr.length >= args["min_tr_len"]:
                f.write(">ITR_" + genome.id + "\n" + genome.itr.seq + "\n")

    logger.info("[4/5] Running dustmasker to check repeat sequence complexity...")
    args["dustmaskerout"] = os.path.join(args["tmp"], "dustmasker.txt")
    utility.run_dustmasker(args["tr_path"], args["dustmaskerout"])
    for line in open(args["dustmaskerout"]):
        row, start, end = line[1:].split()
        type, genome_id = row.split("_", 1)
        if type == "DTR":
            genomes[genome_id].dtr.dust = int(end) - int(start) + 1
        elif type == "ITR":
            genomes[genome_id].itr.dust = int(end) - int(start) + 1

    logger.info("[5/5] Writing results...")
    out = open(args["output"] + "/terminal_repeats.tsv", "w")
    header = ["contig_id", "contig_length", "repeat_type", "repeat_length", "repeat_count", "repeat_dust_length"]
    out.write("\t".join(header) + "\n")
    for genome in genomes.values():
        if (
            genome.length >= args["min_contig_len"]
            and genome.dtr.length >= args["min_tr_len"]
            and (genome.dtr.count <= args["max_tr_count"] )
            and (100.0 * genome.dtr.dust / genome.dtr.length <= args["max_tr_dust"])
        ):
            row = [genome.id, genome.length, 'DTR', genome.dtr.length, genome.dtr.count, genome.dtr.dust]
            out.write('\t'.join([str(_) for _ in row])+'\n')
        elif (
            genome.length >= args["min_contig_len"]
            and genome.itr.length >= args["min_tr_len"]
            and (genome.itr.count <= args["max_tr_count"] )
            and (100.0 * genome.itr.dust / genome.itr.length <= args["max_tr_dust"])
        ):
            row = [genome.id, genome.length, 'ITR', genome.itr.length, genome.itr.count, genome.itr.dust]
            out.write('\t'.join([str(_) for _ in row])+'\n')

    logger.info("\nDone!")
    logger.info("Run time: %s seconds" % round(time.time()-program_start,2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(),2))

