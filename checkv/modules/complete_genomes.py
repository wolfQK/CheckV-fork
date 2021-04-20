import csv
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
        "--tr_min_len",
        type=int,
        default=20,
        metavar="INT",
        help="Min length of TR (20)",
    )
    parser.add_argument(
        "--tr_max_count",
        type=int,
        default=8,
        metavar="INT",
        help="Max occurences of TR per contig (8)",
    )
    parser.add_argument(
        "--tr_max_ambig",
        type=float,
        default=0.20,
        metavar="FLOAT",
        help="Max fraction of TR composed of Ns (0.20)",
    )
    parser.add_argument(
        "--tr_max_basefreq",
        type=float,
        default=0.75,
        metavar="FLOAT",
        help="Max fraction of TR composed of single nucleotide (0.75)",
    )
    parser.add_argument(
        "--kmer_max_freq",
        type=float,
        default=1.5,
        metavar="FLOAT",
        help="Max kmer frequency (1.5). Computed by splitting genome into kmers, counting occurence of each kmer, and taking the average count. Expected value of 1.0 for no duplicated regions; 2.0 for the same genome repeated back-to-back",
    )
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Suppress logging messages",
    )


def set_defaults(args):
    key_values = [
        ("tr_min_len", 20),
        ("tr_max_count", 8),
        ("tr_max_ambig", 0.2),
        ("tr_max_basefreq", 0.70),
        ("kmer_max_freq", 1.5),
    ]
    for key, value in key_values:
        if key not in args:
            args[key] = value


def fetch_dtr(seq, min_length):
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


def fetch_itr(seq, min_len, max_len=1000):
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


def determine_confidence(genome):

    if genome.tr.type is not None and genome.flagged:
        return "low", genome.reason
    if genome.aai_completeness and genome.aai_completeness > 90:
        return "high", "AAI-based completeness > 90%"
    if genome.aai_completeness and genome.aai_completeness > 80:
        return "medium", "AAI-based completeness > 80%"
    elif genome.aai_completeness and genome.aai_completeness < 80:
        return "low", "AAI-based completeness < 80%"
    if genome.hmm_completeness and genome.hmm_completeness[0] > 90:
        return "high", "HMM-based completeness (lower-bound) > 90%"
    elif genome.hmm_completeness and genome.hmm_completeness[0] > 80:
        return "medium", "HMM-based completeness (lower-bound) > 80%"
    elif genome.hmm_completeness and genome.hmm_completeness[1] < 80:
        return "low", "HMM-based completeness (upper-bound) < 80%"
    return "not-determined", "NA"


def main(args):

    program_start = time.time()
    logger = utility.get_logger(args["quiet"])
    set_defaults(args)
    args["tmp"] = os.path.join(args["output"], "tmp")
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    for file in ["completeness.tsv", "contamination.tsv"]:
        path = os.path.join(args["output"], file)
        if not os.path.exists(path):
            sys.exit(f"Error: input file does not exist: {path}\n")

    logger.info(f"\nCheckV v{checkv.__version__}: complete_genomes")

    logger.info("[1/7] Reading input sequences...")
    genomes = {}
    for r in Bio.SeqIO.parse(args["input"], "fasta"):
        if len(r.seq) == 0:
            continue
        genome = Genome()
        genome.id = r.id
        genome.seq = str(r.seq).upper()
        genome.length = len(genome.seq)
        genome.aai_completeness = None
        genome.hmm_completeness = None
        genomes[genome.id] = genome

    logger.info("[2/7] Finding complete proviruses...")
    path = os.path.join(args["output"], "contamination.tsv")
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["region_types"] == "host,viral,host":
            genomes[r["contig_id"]].complete_provirus = "Yes"
        else:
            genomes[r["contig_id"]].complete_provirus = "No"

    logger.info("[3/7] Finding direct/inverted terminal repeats...")
    for genome in genomes.values():
        genome.tr = TR()
        dtr = fetch_dtr(genome.seq, args["tr_min_len"])
        itr = fetch_itr(genome.seq, args["tr_min_len"])
        if len(dtr) < args["tr_min_len"] and len(itr) < args["tr_min_len"]:
            genome.tr.type = None
        elif len(dtr) >= len(itr):
            genome.tr.type = "DTR"
            genome.tr.seq = dtr
            genome.tr.length = len(dtr)
        else:
            genome.tr.type = "ITR"
            genome.tr.seq = itr
            genome.tr.length = len(itr)
    logger.info("[4/7] Filtering terminal repeats...")
    # mode_freq, n_freq, repeat counts
    for genome in genomes.values():
        if genome.tr.type is not None:
            mode_base, mode_count = collections.Counter(genome.tr.seq).most_common(1)[0]
            genome.tr.mode_freq = 1.0 * mode_count / len(genome.tr.seq)
            genome.tr.n_freq = 1.0 * genome.tr.seq.count("N") / len(genome.tr.seq)
            genome.tr.count = genome.seq.count(genome.tr.seq)
    # write trs to fasta
    args["tr_path"] = f"{args['tmp']}/tr.fna"
    with open(args["tr_path"], "w") as f:
        for genome in genomes.values():
            if genome.tr.type is not None:
                f.write(
                    ">" + genome.tr.type + "_" + genome.id + "\n" + genome.tr.seq + "\n"
                )
    # flag repeats based on user parameters
    for genome in genomes.values():
        if genome.tr.type is not None:
            if genome.tr.n_freq > args["tr_max_ambig"]:
                genome.flagged = True
                genome.reason = "Too many ambiguous bases in TR"
            elif genome.tr.count > args["tr_max_count"]:
                genome.flagged = True
                genome.reason = "Repetetive TR sequence"
            elif genome.tr.mode_freq > args["tr_max_basefreq"]:
                genome.flagged = True
                genome.reason = "Low complexity TR"
            else:
                genome.flagged = False

    logger.info("[5/7] Checking genome for completeness...")
    path = os.path.join(args["output"], "completeness.tsv")
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["aai_confidence"] in ("medium", "high"):
            genomes[r["contig_id"]].aai_completeness = float(r["aai_completeness"])
        if int(r["hmm_num_hits"]) > 0:
            genomes[r["contig_id"]].hmm_completeness = (
                float(r["hmm_completeness_lower"]),
                float(r["hmm_completeness_upper"]),
            )

    logger.info("[6/7] Checking genome for large duplications...")
    path = os.path.join(args["output"], "completeness.tsv")
    for r in csv.DictReader(open(path), delimiter="\t"):
        genomes[r["contig_id"]].kmer_freq = float(r["kmer_freq"])
        if float(r["kmer_freq"]) > args["kmer_max_freq"]:
            genomes[r["contig_id"]].flagged = True
            genomes[r["contig_id"]].reason = "Multiple genome copies detected"

    logger.info("[7/7] Writing results...")
    out = open(args["output"] + "/complete_genomes.tsv", "w")
    header = [
        "contig_id",
        "contig_length",
        "kmer_freq",
        "prediction_type",
        "confidence_level",
        "confidence_reason",
        "repeat_length",
        "repeat_count",
        "repeat_n_freq",
        "repeat_mode_base_freq",
        "repeat_seq",
    ]
    out.write("\t".join(header) + "\n")
    for genome in genomes.values():

        if genome.tr.type is not None:
            type = genome.tr.type
            repeat_length = genome.tr.length
            repeat_count = genome.tr.count
            repeat_nfreq = genome.tr.n_freq
            repeat_modefreq = genome.tr.mode_freq
            repeat_seq = genome.tr.seq
            kmer_freq = genome.kmer_freq
        elif genome.complete_provirus == "Yes":
            type = "Provirus"
            repeat_length = "NA"
            repeat_count = "NA"
            repeat_nfreq = "NA"
            repeat_modefreq = "NA"
            repeat_seq = "NA"
            kmer_freq = genome.kmer_freq
        else:
            continue

        confidence, evidence = determine_confidence(genome)
        row = [
            genome.id,
            genome.length,
            kmer_freq,
            type,
            confidence,
            evidence,
            repeat_length,
            repeat_count,
            repeat_nfreq,
            repeat_modefreq,
            repeat_seq,
        ]
        out.write("\t".join([str(_) for _ in row]) + "\n")

    logger.info("Run time: %s seconds" % round(time.time() - program_start, 2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
