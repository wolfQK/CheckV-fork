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
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--max_aai",
        type=float,
        default=None,
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
        "--exclude_circular",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--exclude_genbank",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Suppress logging messages",
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

    program_start = time.time()
    logger = utility.get_logger(args["quiet"])
    utility.check_executables(["prodigal", "diamond"])
    args["db"] = utility.check_database(args["db"])
    args["tmp"] = os.path.join(args["output"], "tmp")

    if args["restart"] and os.path.exists(args["tmp"]):
        shutil.rmtree(args["tmp"])
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    args["faa"] = os.path.join(args["tmp"], "proteins.faa")
    if not os.path.exists(args["faa"]):
        logger.info("[1/7] Calling genes with prodigal...")
        utility.call_genes(args["input"], args["output"], args["threads"])
    else:
        logger.info("[1/7] Skipping gene calling...")

    args["blastp"] = os.path.join(args["tmp"], "diamond.tsv")
    if not os.path.exists(args["blastp"]):
        logger.info("[2/7] Running DIAMOND blastp search...")
        utility.run_diamond(args["blastp"], args["db"], args["faa"], args["threads"])
    else:
        logger.info("[2/7] Skipping DIAMOND blastp search...")

    logger.info("[3/7] Initializing queries and database...")

    # read fna input
    genomes = {}
    for header, seq in utility.read_fasta(args["input"]):
        genome = Genome()
        genome.id = header.split()[0]
        genome.length = len(seq)
        genome.genes = 0
        genome.protlen = 0
        genome.viral_length = None
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

    # read prophage predictions, if exists
    p = os.path.join(args["output"], "contamination.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            if "genome_id" in r: ## temporary fix due to changing field format
                r["contig_id"] = r["genome_id"]
            if "viral" in r["region_types"]:
                genomes[r["contig_id"]].viral_length = int(r["viral_length"])
            else:
                genomes[r["contig_id"]].viral_length = None

    # read reference genomes
    refs = {}
    p = os.path.join(args["db"], "genome_db/checkv_reps.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        genome = Genome()
        genome.id = r["checkv_id"]
        genome.length = int(r["length"])
        genome.type = r["type"]
        genome.weight = 1 if r["type"] == "circular" else 2 if r["type"] == "genbank" else None
        refs[genome.id] = genome

    # exclude references
    exclude = set([])
    if args["exclude_list"]:
        for l in open(args["exclude_list"]):
            exclude.add(l.rstrip())
    for genome in refs.values():
        if args["exclude_circular"] and genome.type == "circular":
            exclude.add(genome.id)
        elif args["exclude_genbank"] and genome.type == "genbank":
            exclude.add(genome.id)

    # estimated error rates for alignment cutoffs
    error_rates = {}
    p = os.path.join(args["db"], "genome_db/checkv_error.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        key = int(r["length"]), int(r["aai"]), int(r["cov"])
        error_rates[key] = float(r["error"]) if int(r["count"]) > 100 else "NA"
    error_keys = {}
    error_keys["length"] = sorted(list(set([_[0] for _ in error_rates.keys()])))
    error_keys["aai"] = sorted(list(set([_[1] for _ in error_rates.keys()])))
    error_keys["cov"] = sorted(list(set([_[2] for _ in error_rates.keys()])))

    args["aai"] = os.path.join(args["tmp"], "aai.tsv")
    if not os.path.exists(args["aai"]):
        logger.info("[4/7] Computing AAI...")
        compute_aai(args["blastp"], args["aai"], genomes, refs)
    else:
        logger.info("[4/7] Skipping AAI computation...")

    logger.info("[5/7] Storing AAI...")
    for r in csv.DictReader(open(args["aai"]), delimiter="\t"):

        # format to floats
        r["identity"] = float(r["identity"])
        r["aligned_genes"] = float(r["aligned_genes"])
        r["aligned_length"] = int(r["aligned_length"])
        r["score"] = float(r["score"])
        r["percent_length"] = float(r["percent_length"])

        if args["max_aai"] is not None and r["identity"] > args["max_aai"]:
            continue
        elif r["target"] in exclude:
            continue
        elif r["target"] not in refs:
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

    logger.info("[6/7] Estimating completeness...")
    for genome in genomes.values():
        rec = {}
        rec["contig_id"] = genome.id
        rec["contig_length"] = genome.length

        if len(genome.aai) > 0:
            # fetch top hit
            rec["ref_name"] = genome.aai[0]["target"]
            rec["ref_aai"] = genome.aai[0]["identity"]
            rec["ref_af"] = genome.aai[0]["percent_length"]
            # get estimation error
            key1 = take_closest(error_keys["length"], genome.length)
            key2 = take_closest(error_keys["aai"], genome.aai[0]["identity"])
            key3 = take_closest(error_keys["cov"], genome.aai[0]["percent_length"])
            rec["est_error"] = error_rates[key1, key2, key3]
            rec["confidence"] = "low" if rec["est_error"] == "NA" else "high" if rec["est_error"] <= 5 else "medium" if rec["est_error"] <= 10 else "low"
            # compute expected genome length
            lengths = [refs[_["target"]].length for _ in genome.aai]
            weights = [refs[_["target"]].weight * _["score"] for _ in genome.aai]
            rec["expected_length"] = sum([l * w for l, w in zip(lengths, weights)]) / sum(weights)
            # estimate completness
            if genome.viral_length is not None:
                rec["viral_length"] = genome.viral_length
                rec["completeness"] = 100.0 * genome.viral_length / rec["expected_length"]
            else:
                rec["completeness"] = 100.0 * genome.length / rec["expected_length"]
            # summarize distribution of reference lengths within 50% of top hit
            len_75, len_25 = np.percentile(lengths, [75, 25])
            len_min, len_max = min(lengths), max(lengths)
            len_med = np.median(lengths)
            rec["ref_length_distribution"] = ",".join([str(len_min), str(len_25), str(len_med), str(len_75), str(len_max)])
            rec["ref_hits"] = len(genome.aai)
        genome.rec = rec

    logger.info("[7/7] Writing results...")
    p = os.path.join(args["output"], "completeness.tsv")
    fields = [
        "contig_id",
        "contig_length",
        "viral_length",
        "expected_length",
        "completeness",
        "confidence",
        "ref_name",
        "ref_aai",
        "ref_af",
        "ref_hits",
        "ref_length_distribution"
        ]
    with open(p, "w") as out:
        out.write("\t".join(fields) + "\n")
        for genome in genomes.values():
            row = [str(genome.rec[f]) if f in genome.rec else 'NA' for f in fields]
            out.write("\t".join(row) + "\n")

    # done!
    logger.info("\nDone!")
    logger.info("Run time: %s seconds" % round(time.time()-program_start,2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(),2))
