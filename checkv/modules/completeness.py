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
        "input", type=str, help="Input nucleotide sequences in FASTA format",
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
        "--percent_of_top_hit", type=float, default=50.0, help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--max_aai", type=float, default=None, help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--exclude_identical",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--exclude_list", type=str, default=None, metavar="PATH", help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Suppress logging messages",
    )


def set_defaults(args):
    key_values = [
        ("percent_of_top_hit", 50),
        ("max_aai", None),
        ("exclude_identical", False),
        ("exclude_list", None),
    ]
    for key, value in key_values:
        if key not in args:
            args[key] = value


def yield_query_alns(path):
    """
    Yields list of formatted blastp alignments per 'qname'
    """

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


def init_database(args):
    """
    Initialize genes, genomes, refs, and exclude data objects
    genes: dictionary of gene_id to gene object (from input)
    genomes: dictionary of genome_id to genome object (from input)
    refs: dictionary of checkv_id to genome object (from database)
    exclude: set of user-defined checkv_ids to exclude from AAI-based estimation
    """

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
        genome.hmms = set([])
        genome.hmm_name = None
        genome.hmm_completeness_lower = None
        genome.hmm_completeness_upper = None
        genome.expected_length = None
        genome.aai_completeness = None
        genome.aai_confidence = None
        genome.aai_error = None
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
            if "genome_id" in r:  ## temporary fix due to changing field format
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
        genome.weight = (
            1 if r["type"] == "circular" else 2 if r["type"] == "genbank" else None
        )
        refs[genome.id] = genome

    # exclude references
    exclude = set([])
    if args["exclude_list"]:
        for l in open(args["exclude_list"]):
            exclude.add(l.rstrip())

    # read hmm info
    hmms = {}
    p = os.path.join(args["db"], "hmm_db/genome_lengths.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        r["genomes"] = int(r["genomes"])
        r["cv"] = float(r["cv"])
        r["lengths"] = [int(_) for _ in r["lengths"].split(",")]
        hmms[r["hmm"]] = r

    return genes, genomes, refs, exclude, hmms


def compute_aai(blastp_path, out_path, genomes, genes, refs):
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
                lens = [genes[_].length for _ in alns]
                aligned_length = sum(lens)
                aligned_genes = len(alns)
                identity = round(
                    sum(x * y for x, y in zip(pids, lens)) / aligned_length, 2
                )
                percent_length = round(
                    100.0 * aligned_length / genomes[query].protlen, 2
                )
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


def store_aai(args, genomes, exclude, refs):
    """
    1. Loop over AAI records from disk
    2. Format numeric values
    3. Filter records based on exclusion criterea
    4. For each query, store all AAI records with alignment scores within 50% of top hit
    """

    for r in csv.DictReader(open(args["aai"]), delimiter="\t"):

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

        if len(genomes[r["query"]].aai) == 0:
            genomes[r["query"]].aai.append(r)
        else:
            top_score = genomes[r["query"]].aai[0]["score"]
            if (
                100.0 * (top_score - r["score"]) / top_score
                <= args["percent_of_top_hit"]
            ):
                genomes[r["query"]].aai.append(r)


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


def lookup_error_rate(genome, error_rates, error_keys):
    """
    Returns estimated error rate for genome based on AAI, AF, and length
    """

    key1 = take_closest(error_keys["length"], genome.length)
    key2 = take_closest(error_keys["aai"], genome.aai[0]["identity"])
    key3 = take_closest(error_keys["cov"], genome.aai[0]["percent_length"])
    error = error_rates[key1, key2, key3]
    confidence = (
        "low"
        if error == "NA"
        else "high"
        if error <= 5
        else "medium"
        if error <= 10
        else "low"
    )
    return error, confidence


def aai_based_completeness(args, genomes, exclude, refs):
    """
    0. Store error rates
    1. Fetch top database hit
    2. Get estimation error based on AAI and AF to top hit as well as contig length
    3. Compute expected genome length as weighted average of all hits with an alignment score within 50% of top hit
    4. Compute completeness as ratio of contig length (or viral length for proviruses) to expected length
    """

    # store aai
    store_aai(args, genomes, exclude, refs)

    # store error rates
    error_rates = {}
    p = os.path.join(args["db"], "genome_db/checkv_error.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        key = int(r["length"]), int(r["aai"]), int(r["cov"])
        error_rates[key] = float(r["error"]) if int(r["count"]) > 100 else "NA"

    # list keys for lookup
    error_keys = {
        "length": sorted(list({_[0] for _ in error_rates})),
        "aai": sorted(list({_[1] for _ in error_rates})),
        "cov": sorted(list({_[2] for _ in error_rates})),
    }

    # compute completeness
    for genome in genomes.values():
        if len(genome.aai) > 0:
            est_error, confidence = lookup_error_rate(genome, error_rates, error_keys)
            lengths = [refs[_["target"]].length for _ in genome.aai]
            weights = [refs[_["target"]].weight * _["score"] for _ in genome.aai]
            genome.expected_length = sum(l * w for l, w in zip(lengths, weights)) / sum(
                weights
            )
            query_length = (
                genome.viral_length
                if genome.viral_length is not None
                else genome.length
            )
            genome.aai_completeness = 100.0 * query_length / genome.expected_length
            genome.aai_confidence = confidence
            genome.aai_error = est_error


def fetch_percentile(num, lst):
    """
    Compute the percentile of num in lst. Example logic:
    if num is smaller than 1st element, rank = 0, percentile = 0
    if num is smaller than 2nd element, rank = 1, percentile = 1/len(list)
    if num is not smaller than last element, percentile = 100
    """
    for rank, length in enumerate(lst):
        if num < length:
            return 100.0 * (rank) / len(lst)
    return 100.0


def hmm_based_completeness(args, genomes, hmms):
    """
    1. Store HMM information (name, genome lengths, CV)
    2. List HMMs hitting each genome
    3. Identify the HMM whose database targets display the least genome size variation (i.e. smallest CV)
    4. Compare the contig length (or viral region on proviruses) to the database target lengths
    5. Estimate the 90% confidence interval of completeness
    """

    # list hmms hitting each genome
    p = os.path.join(args["tmp"], "gene_annotations.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["target_hmm"] in hmms and hmms[r["target_hmm"]]["genomes"] >= 10:
            genomes[r["contig_id"]].hmms.add(r["target_hmm"])

    # estimate minimum completeness
    for genome in genomes.values():

        # pick hmm for each query with lowest cv
        genome.hmms = list(genome.hmms)
        if len(genome.hmms) == 0:
            continue
        cvs = [hmms[_]["cv"] for _ in genome.hmms]
        hmm = [hmm for cv, hmm in sorted(zip(cvs, genome.hmms))][0]
        query_length = (
            genome.viral_length if genome.viral_length is not None else genome.length
        )

        # rank genome based on hypothetical complenetness
        # ranging from 1 to 100
        x = []
        for comp in range(0, 101, 1):
            size = query_length * 100.0 / comp if comp > 0 else float("Inf")
            perc = fetch_percentile(size, hmms[hmm]["lengths"])
            x.append([comp, perc])

        # completness_lower: 95% probability true value is higher
        # completness_upper: 95% probability true value is not higher
        genome.hmm_completeness_lower = [comp for comp, perc in x if perc >= 95][-1]
        genome.hmm_completeness_upper = [comp for comp, perc in x if perc >= 5][-1]
        genome.hmm_name = hmm


def main(args):

    program_start = time.time()
    set_defaults(args)
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

    logger.info(f"\nCheckV v{checkv.__version__}: completeness")

    args["faa"] = os.path.join(args["tmp"], "proteins.faa")
    if os.path.exists(args["faa"]):
        logger.info("[1/7] Skipping gene calling...")
    else:
        logger.info("[1/7] Calling genes with Prodigal...")
        utility.call_genes(args["input"], args["output"], args["threads"])

    logger.info("[2/7] Initializing queries and database...")
    genes, genomes, refs, exclude, hmms = init_database(args)

    args["blastp"] = os.path.join(args["tmp"], "diamond.tsv")
    if os.path.exists(args["blastp"]):
        logger.info("[3/7] Skipping DIAMOND blastp search...")
    else:
        logger.info("[3/7] Running DIAMOND blastp search...")
        utility.run_diamond(args["blastp"], args["db"], args["faa"], args["threads"])

    args["aai"] = os.path.join(args["tmp"], "aai.tsv")
    if os.path.exists(args["aai"]):
        logger.info("[4/7] Skipping AAI computation...")
    else:
        logger.info("[4/7] Computing AAI...")
        compute_aai(args["blastp"], args["aai"], genomes, genes, refs)

    logger.info("[5/7] Running AAI based completeness estimation...")
    aai_based_completeness(args, genomes, exclude, refs)

    annotation_path = os.path.join(args["tmp"], "gene_annotations.tsv")
    if os.path.exists(annotation_path):
        logger.info("[6/7] Running HMM based completeness estimation...")
        hmm_based_completeness(args, genomes, hmms)
    else:
        logger.info(
            "[6/7] Skipping HMM based completeness estimation: requires output from 'checkv contamination'"
        )

    logger.info("[7/7] Writing results...")
    p = os.path.join(args["output"], "completeness.tsv")
    fields = [
        "contig_id",
        "contig_length",
        "viral_length",
        "aai_expected_length",
        "aai_completeness",
        "aai_confidence",
        "aai_error",
        "aai_num_hits",
        "aai_top_hit",
        "aai_id",
        "aai_af",
        "hmm_completeness_lower",
        "hmm_completeness_upper",
        "hmm_name",
        "hmm_ref_genomes",
        "hmm_avg_length",
        "hmm_cv_length",
    ]
    with open(p, "w") as out:
        out.write("\t".join(fields) + "\n")
        for genome in genomes.values():
            rec = {
                "contig_id": genome.id,
                "contig_length": genome.length,
                "viral_length": genome.viral_length,
                "aai_expected_length": genome.expected_length,
                "aai_completeness": genome.aai_completeness,
                "aai_confidence": genome.aai_confidence,
                "aai_error": genome.aai_error,
                "aai_num_hits": len(genome.aai),
            }
            if len(genome.aai) > 0:
                rec["aai_top_hit"] = genome.aai[0]["target"]
                rec["aai_id"] = genome.aai[0]["identity"]
                rec["aai_af"] = genome.aai[0]["percent_length"]
            rec["hmm_completeness_lower"] = genome.hmm_completeness_lower
            rec["hmm_completeness_upper"] = genome.hmm_completeness_upper
            rec["hmm_name"] = genome.hmm_name
            if genome.hmm_name:
                rec["hmm_ref_genomes"] = hmms[genome.hmm_name]["genomes"]
                rec["hmm_avg_length"] = round(
                    np.mean(hmms[genome.hmm_name]["lengths"]), 1
                )
                rec["hmm_cv_length"] = hmms[genome.hmm_name]["cv"]
            row = [
                str(rec[f]) if f in rec and rec[f] is not None else "NA" for f in fields
            ]
            out.write("\t".join(row) + "\n")

    # done!
    logger.info("Run time: %s seconds" % round(time.time() - program_start, 2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
