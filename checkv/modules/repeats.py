import argparse
import csv
import logging
import os
import string
import sys
import time
import Bio.SeqIO
import numpy as np
import checkv
import collections
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
    parser.set_defaults(program="repeats")
    parser.add_argument(
        "input", type=str, help="Input viral sequences in FASTA format",
    )
    parser.add_argument("output", type=str, help="Output directory")
    parser.add_argument(
        "--min_tr_len", type=int, default=20, metavar="INT", help="Min length of TR",
    )
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Suppress logging messages",
    )


def set_defaults(args):
    key_values = [("min_tr_len", 20), ("max_tr_count", 5), ("max_tr_dust", 20.0)]
    for key, value in key_values:
        if key not in args:
            args[key] = value


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
    if sys.version_info > (3, 0):
        trans = str.maketrans("ACTG", "TGAC")
    else:
        trans = string.maketrans("ACTG", "TGAC")
    return seq[::-1].translate(trans)


def fetch_itr(seq, min_len=20, max_len=1000):
    rev = reverse_complement(seq)
    # see if minimal substring occurs at end
    if seq[:min_len] == rev[:min_len]:
        # extend to maximum substring, up to <max_len>
        i = min_len + 1
        while seq[:i] == rev[:i] and i <= max_len:
            i += 1
        return seq[: i - 1]
    # no match
    else:
        return ""


def main(args):

    program_start = time.time()
    logger = utility.get_logger(args["quiet"])
    set_defaults(args)
    utility.check_executables(["dustmasker"])
    args["tmp"] = os.path.join(args["output"], "tmp")
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    logger.info(f"\nCheckV v{checkv.__version__}: repeats")

    logger.info("[1/6] Reading input sequences...")
    genomes = {}
    for r in Bio.SeqIO.parse(args["input"], "fasta"):
        genome = Genome()
        genome.id = r.id
        genome.seq = str(r.seq).upper()
        genome.length = len(genome.seq)
        genomes[genome.id] = genome

    logger.info("[2/6] Determining genome copy number...")
    # at least 20 windows with max win size of 2000-bp
    for genome in genomes.values():
        counts = []
        win_size = min([int(len(genome.seq) / 20), 2000])
        start, end = 0, win_size
        while end <= len(genome.seq):
            counts.append(genome.seq.count(genome.seq[start:end]))
            start += win_size
            end += win_size
        genome.avg_copies = round(np.mean(counts), 2)

    logger.info("[3/6] Finding direct/inverted terminal repeats...")
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

    logger.info("[4/6] Determining repeat base composition...")
    for id in sorted(genomes.keys()):
        genome = genomes[id]
        if len(genome.dtr.seq) >= args["min_tr_len"]:
            mode_base, mode_count = collections.Counter(genome.dtr.seq).most_common(1)[
                0
            ]
            genome.dtr.mode_freq = 1.0 * mode_count / len(genome.dtr.seq)
            genome.dtr.n_freq = 1.0 * genome.dtr.seq.count("N") / len(genome.dtr.seq)
        if len(genome.itr.seq) >= args["min_tr_len"]:
            mode_base, mode_count = collections.Counter(genome.itr.seq).most_common(1)[
                0
            ]
            genome.itr.mode_freq = 1.0 * mode_count / len(genome.itr.seq)
            genome.itr.n_freq = 1.0 * genome.itr.seq.count("N") / len(genome.itr.seq)

    logger.info("[5/6] Running dustmasker to check repeat sequence complexity...")
    args["tr_path"] = f"{args['tmp']}/tr.fna"
    with open(args["tr_path"], "w") as f:
        for id in sorted(genomes.keys()):
            genome = genomes[id]
            if genome.dtr.length >= args["min_tr_len"]:
                f.write(">DTR_" + genome.id + "\n" + genome.dtr.seq + "\n")
            if genome.itr.length >= args["min_tr_len"]:
                f.write(">ITR_" + genome.id + "\n" + genome.itr.seq + "\n")
    args["dustmaskerout"] = os.path.join(args["tmp"], "dustmasker.txt")
    utility.run_dustmasker(args["tr_path"], args["dustmaskerout"])
    for line in open(args["dustmaskerout"]):
        row, start, end = line[1:].split()
        type, genome_id = row.split("_", 1)
        if type == "DTR":
            genomes[genome_id].dtr.dust = int(end) - int(start) + 1
        elif type == "ITR":
            genomes[genome_id].itr.dust = int(end) - int(start) + 1

    logger.info("[6/6] Writing results...")
    out = open(args["output"] + "/repeats.tsv", "w")
    header = [
        "contig_id",
        "contig_length",
        "genome_copies",
        "repeat_type",
        "repeat_length",
        "repeat_count",
        "repeat_dust_length",
        "repeat_mode_base_freq",
        "repeat_ambig_base_freq",
        "repeat_flagged",
        "flagged_reason",
        "repeat_sequence",
    ]
    out.write("\t".join(header) + "\n")
    for genome in genomes.values():

        if genome.dtr.length >= args["min_tr_len"]:
            tr_type = "DTR"
            tr_seq = genome.dtr.seq
            tr_length = genome.dtr.length
            tr_count = genome.dtr.count
            tr_dust = genome.dtr.dust
            tr_ambig = genome.dtr.n_freq
            tr_mode_base_freq = genome.dtr.mode_freq
            if tr_mode_base_freq >= 0.75:
                tr_flag = "Yes"
                reason = "high_single_base_freq"
            elif tr_ambig >= 0.10:
                tr_flag = "Yes"
                reason = "high_ambig_base_freq"
            elif tr_count > args["max_tr_count"]:
                tr_flag = "Yes"
                reason = "repetetive"
            elif 100.0 * tr_dust / tr_length > args["max_tr_dust"]:
                tr_flag = "Yes"
                reason = "low_complexity"
            else:
                tr_flag = "No"
                reason = "NA"

        elif genome.itr.length >= args["min_tr_len"]:
            tr_type = "ITR"
            tr_seq = genome.itr.seq
            tr_length = genome.itr.length
            tr_count = genome.itr.count
            tr_dust = genome.itr.dust
            tr_ambig = genome.itr.n_freq
            tr_mode_base_freq = genome.itr.mode_freq
            if tr_mode_base_freq >= 0.75:
                tr_flag = "Yes"
                reason = "high_single_base_freq"
            elif tr_ambig >= 0.10:
                tr_flag = "Yes"
                reason = "high_ambig_base_freq"
            elif tr_count > args["max_tr_count"]:
                tr_flag = "Yes"
                reason = "repetetive"
            elif 100.0 * tr_dust / tr_length > args["max_tr_dust"]:
                tr_flag = "Yes"
                reason = "low_complexity"
            else:
                tr_flag = "No"
                reason = "NA"

        else:
            tr_type = "NA"
            tr_seq = "NA"
            tr_length = "NA"
            tr_count = "NA"
            tr_dust = "NA"
            tr_ambig = "NA"
            tr_mode_base_freq = "NA"
            tr_flag = "NA"
            reason = "NA"

        row = [
            genome.id,
            genome.length,
            genome.avg_copies,
            tr_type,
            tr_length,
            tr_count,
            tr_dust,
            tr_mode_base_freq,
            tr_ambig,
            tr_flag,
            reason,
            tr_seq,
        ]
        out.write("\t".join([str(_) for _ in row]) + "\n")

    logger.info("Run time: %s seconds" % round(time.time() - program_start, 2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
