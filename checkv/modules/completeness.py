#!/usr/bin/env python

import argparse
import bisect
import csv
import logging
import operator
import os
import subprocess as sp
import sys
import time

import numpy as np
from checkv import utility


class Genome:
    def __init__(self):
        pass


class Gene:
    def __init__(self):
        pass


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="completeness")
    parser.add_argument(
        "input",
        type=str,
        help="Input nucleotide sequences in FASTA format",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Output directory"
    )
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
        "--percent_of_top_hit",
        type=float,
        default=50,
        metavar="FLOAT",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--exclude_identical",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--exclude_list",
        type=str,
        default=None,
        metavar="PATH",
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Display logging messages",
    )


def yield_query_alns(path):
    handle = utility.parse_blastp(path)
    alns = [next(handle)]
    query = alns[0]["qname"].rsplit("_", 1)[0]
    for r in handle:
        if r["qname"].rsplit("_", 1)[0] != query:
            yield query, alns
            alns = []
            query = r["qname"].rsplit("_", 1)[0]
        alns.append(r)
    yield query, alns


def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.
    If two numbers are equally close, return the smallest number.
    """
    myNumber = float(myNumber)
    pos = bisect.bisect_left(myList, float(myNumber))
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if (after - myNumber) < (myNumber - before):
        return after
    else:
        return before


def compute_aai(blastp_path, out_path, genomes, refs):
    """ Compute AAI to reference genomes
    1. loop over alignment blocks from <blastp_path> (1 per query)
    2. compute average aa identity to each reference genome
    3. write all hits to disk
    """
    with open(out_path, "w") as out:
        header = [
            "query",
            "target",
            "aligned_genes",
            "percent_genes",
            "aligned_length",
            "percent_length",
            "identity",
            "score",
        ]
        out.write("\t".join(header) + "\n")

        for query, alns in yield_query_alns(blastp_path):

            # get best hits for each gene to each target genome
            hits = {}
            for r in alns:
                target = r["tname"].rsplit("_", 1)[0]
                # store gene aln
                if target not in hits:
                    hits[target] = {}
                if r["qname"] not in hits[target]:
                    hits[target][r["qname"]] = r
                elif r["score"] > hits[target][r["qname"]]["score"]:
                    hits[target][r["qname"]] = r

            # compute aai
            aai = []
            for target, alns in hits.items():
                pids = [r["pid"] for r in alns.values()]
                lens = [r["aln"] - r["gap"] for r in alns.values()]
                aligned_length = sum(lens)
                aligned_genes = len(alns)
                identity = round(
                    sum(x * y for x, y in zip(pids, lens)) / aligned_length, 2
                )
                percent_length = round(100.0 * aligned_length / genomes[query].protlen, 2)
                percent_genes = round(100.0 * aligned_genes / genomes[query].genes, 2)
                score = round(identity * aligned_length / 100, 2)
                row = [
                    query,
                    target,
                    aligned_genes,
                    percent_genes,
                    aligned_length,
                    percent_length,
                    identity,
                    score,
                ]
                aai.append(row)
            if not aai:
                continue

            # write all hits to disk
            aai = sorted(aai, key=operator.itemgetter(-1), reverse=True)
            top_score = max(_[-1] for _ in aai)
            for row in aai:
                score = row[-1]
                out.write("\t".join([str(_) for _ in row]) + "\n")


def main(args):
    logger = utility.get_logger(args["verbose"])
    utility.check_executables(["prodigal", "diamond"])
    args["db"] = utility.check_database(args["db"])
    args["tmp"] = os.path.join(args["output"], "tmp")
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    logger.info("[1/5] Calling genes with prodigal…")
    args["faa"] = os.path.join(args["tmp"], "proteins.faa")
    if args["restart"] or not os.path.exists(args["faa"]):
        utility.call_genes(args["input"], args["output"], args["threads"])

    logger.info("[2/5] Running DIAMOND blastp search…")
    args["blastp"] = os.path.join(args["tmp"], "diamond.tsv")
    if args["restart"] or not os.path.exists(args["blastp"]):
        utility.run_diamond(args["blastp"], args["db"], args["faa"], args["threads"])

    logger.info("[3/5] Initializing queries and databas…")
    genomes = {}
    for header, seq in utility.read_fasta(args["input"]):
        genome = Genome()
        genome.id = header.split()[0]
        genome.length = len(seq)
        genome.genes = 0
        genome.protlen = 0
        genome.aai = []
        genomes[genome.id] = genome
    # read input genes
    genes = {}
    for header, seq in utility.read_fasta(args["faa"]):
        gene = Gene()
        gene.id = header.split()[0]
        gene.length = len(seq)
        gene.genome_id = gene.id.rsplit("_", 1)[0]
        genes[gene.id] = gene
        genomes[gene.genome_id].genes += 1
        genomes[gene.genome_id].protlen += gene.length
    # read reference genomes
    refs = {}
    p = os.path.join(args["db"], "checkv_refs.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        genome = Genome()
        genome.id = r["checkv_id"]
        genome.length = int(r["length"])
        refs[genome.id] = genome

    # exclude references
    exclude = set([])
    for l in open(args["db"] + "/exclude_genomes.list"):
        exclude.add(l.rstrip())
    if args["exclude_list"]:
        for l in open(args["exclude_list"]):
            exclude.add(l.rstrip())

    # estimated error rates for alignment cutoffs
    error_rates = {}
    p = os.path.join(args["db"], "error_rates.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        key = int(r["length"]), int(r["aai"]), int(r["cov"])
        error_rates[key] = float(r["error"]) if int(r["count"]) > 100 else "NA"
    error_keys = {}
    error_keys["length"] = sorted(list(set([_[0] for _ in error_rates.keys()])))
    error_keys["aai"] = sorted(list(set([_[1] for _ in error_rates.keys()])))
    error_keys["cov"] = sorted(list(set([_[2] for _ in error_rates.keys()])))

    logger.info("[4/5] Computing AAI…")
    args["aai"] = os.path.join(args["tmp"], "aai.tsv")
    if args["restart"] or not os.path.exists(args["aai"]):
        compute_aai(args["blastp"], args["aai"], genomes, refs)
    for r in csv.DictReader(open(args["aai"]), delimiter="\t"):

        # better formatting here to floats
        r["identity"] = float(r["identity"])
        r["aligned_genes"] = float(r["aligned_genes"])
        r["aligned_length"] = int(r["aligned_length"])
        r["score"] = float(r["score"])
        r["percent_length"] = float(r["percent_length"])

        if r["target"] in exclude:
            continue
        elif (
            args["exclude_identical"]
            and r["identity"] == 100
            and r["percent_length"] == 100
        ):
            continue
        elif len(genomes[r["query"]].aai) == 0:
            genomes[r["query"]].aai.append(r)
        else:
            top_score = genomes[r["query"]].aai[0]["score"]
            if 100.0 * (top_score - r["score"]) / top_score <= args["percent_of_top_hit"]:
                genomes[r["query"]].aai.append(r)

    logger.info("[5/5] Estimating completeness…")
    module_start = time.time()
    p = os.path.join(args["output"], "completeness.tsv")
    with open(p, "w") as out:
        header = [
            "genome_id",
            "genome_length",
            "num_hits",
            "ref_length_min",
            "ref_length_q1",
            "ref_length_mean",
            "ref_length_q2",
            "ref_length_max",
            "top_hit",
            "top_hit_length",
            "top_hit_aai",
            "top_hit_cov",
            "error_rate",
        ]
        out.write("\t".join(header) + "\n")

        for genome in genomes.values():

            # at least 1 hit to a reference
            if len(genome.aai) > 0:

                # fetch top hit
                top = genome.aai[0]
                key1 = take_closest(error_keys["length"], genome.length)
                key2 = take_closest(error_keys["aai"], top["identity"])
                key3 = take_closest(error_keys["cov"], top["percent_length"])
                error_rate = error_rates[key1, key2, key3]

                # estimate genome size
                scores = [_["score"] for _ in genome.aai]
                lengths = [refs[_["target"]].length for _ in genome.aai]
                avg_len = sum([l * s for l, s in zip(lengths, scores)]) / sum(scores)
                len_75, len_25 = np.percentile(lengths, [75, 25])
                min_len, max_len = min(lengths), max(lengths)
                num_hits = len(lengths)

                # write
                row = [genome.id, genome.length]
                row += [num_hits, min_len, len_25, avg_len, len_75, max_len]
                row += [
                    top["target"],
                    refs[top["target"]].length,
                    top["identity"],
                    top["percent_length"],
                    error_rate,
                ]
                out.write("\t".join([str(_) for _ in row]) + "\n")

            # no hits to any reference
            else:
                row = [genome.id, genome.length] + ["NA"] * 11
                out.write("\t".join([str(_) for _ in row]) + "\n")

    # done!
    logger.info("Done!")