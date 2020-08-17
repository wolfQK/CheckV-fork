import argparse
import csv
import logging
import os
import string
import sys
import time
import checkv
from checkv import utility


class Genome:
    def __init__(self):
        pass


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="quality_summary")
    parser.add_argument(
        "input", type=str, help="Input viral sequences in FASTA format",
    )
    parser.add_argument("output", type=str, help="Output directory")
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Suppress logging messages",
    )


def main(args):

    program_start = time.time()
    logger = utility.get_logger(args["quiet"])
    args["tmp"] = os.path.join(args["output"], "tmp")
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    for file in ["completeness.tsv", "complete_genomes.tsv"]:
        path = os.path.join(args["output"], file)
        if not os.path.exists(path):
            sys.exit(f"Error: input file does not exist: {path}\n")

    logger.info(f"\nCheckV v{checkv.__version__}: quality_summary")

    logger.info("[1/6] Reading input sequences...")
    genomes = {}
    for header, seq in utility.read_fasta(args["input"]):
        genome = Genome()
        genome.id = header.split()[0]
        genome.seq = seq
        genome.length = len(seq)
        genome.proviral_length = "NA"
        genome.kmer_freq = "NA"
        genome.contamination = "NA"
        genome.prophage = "NA"
        genome.total_genes = "NA"
        genome.viral_genes = "NA"
        genome.host_genes = "NA"
        genome.completeness = "NA"
        genome.complete = False
        genome.method = "NA"
        genome.warnings = []
        genomes[genome.id] = genome

    p = os.path.join(args["output"], "contamination.tsv")
    if os.path.exists(p):
        logger.info("[2/6] Reading results from contamination module...")
        for r in csv.DictReader(open(p), delimiter="\t"):
            genome = genomes[r["contig_id"]]
            # num genes
            genome.proviral_length = r["proviral_length"]
            genome.total_genes = r["total_genes"]
            genome.viral_genes = r["viral_genes"]
            genome.host_genes = r["host_genes"]
            # add warnings
            if int(genome.viral_genes) == 0:
                genome.warnings.append("no viral genes detected")
            if r["region_types"].count("viral") > 1:
                genome.warnings.append(">1 viral region detected")
            # prophage
            if r["provirus"] == "Yes":
                genome.prophage = "Yes"
                genome.contamination = round(
                    100.0 * int(r["host_length"]) / genome.length, 2
                )
            else:
                genome.prophage = "No"
                genome.contamination = 0.0
    else:
        logger.info("[2/6] Skipping contamination due to missing input file...")

    logger.info("[3/6] Reading results from completeness module...")
    p = os.path.join(args["output"], "completeness.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            genome = genomes[r["contig_id"]]
            genome.kmer_freq = r["kmer_freq"]
            # add warnings
            if float(genome.kmer_freq) >= 1.5:
                genome.warnings.append("high kmer_freq may indicate large duplication")
            if r["aai_confidence"] in ("high", "medium") and float(
                r["contig_length"]
            ) > 1.5 * float(r["aai_expected_length"]):
                genome.warnings.append(
                    "contig >1.5x longer than expected genome length"
                )
            # AAI-based esimtate (med/high-confidence)
            if r["aai_confidence"] in ("high", "medium"):
                genome.completeness = round(float(r["aai_completeness"]), 2)
                genome.method = "AAI-based (%s-confidence)" % r["aai_confidence"]
            # HMM-based estimate
            elif r["hmm_completeness_lower"] != "NA":
                genome.completeness = round(float(r["hmm_completeness_lower"]), 2)
                genome.method = "HMM-based (lower-bound)"
            # no completeness estimate whatsoever
            else:
                genome.completeness = "NA"
                genome.method = "NA"

    logger.info("[4/6] Reading results from complete genomes module...")
    p = os.path.join(args["output"], "complete_genomes.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        genome = genomes[r["contig_id"]]
        if r["confidence_level"] in ("medium", "high"):
            genome.complete = True
            genome.completeness = 100.0
            genome.method = "%s (%s-confidence)" % (r["prediction_type"], r["confidence_level"])
        else:
            genome.warnings.append(f"low-confidence {r['prediction_type']}")

    logger.info("[5/6] Classifying contigs into quality tiers...")
    for genome in genomes.values():
        if genome.complete:
            genome.quality = "Complete"
            genome.miuvig = "High-quality"
        elif genome.completeness == "NA":
            genome.quality = "Not-determined"
            genome.miuvig = "Genome-fragment"
        elif genome.completeness >= 90:
            genome.quality = "High-quality"
            genome.miuvig = "High-quality"
        elif genome.completeness >= 50:
            genome.quality = "Medium-quality"
            genome.miuvig = "Genome-fragment"
        elif genome.completeness < 50:
            genome.quality = "Low-quality"
            genome.miuvig = "Genome-fragment"

    logger.info("[6/6] Writing results...")
    header = [
        "contig_id",
        "contig_length",
        "provirus",
        "proviral_length",
        "gene_count",
        "viral_genes",
        "host_genes",
        "checkv_quality",
        "miuvig_quality",
        "completeness",
        "completeness_method",
        "contamination",
        "kmer_freq",
        "warnings",
    ]
    out = open(args["output"] + "/quality_summary.tsv", "w")
    out.write("\t".join(header) + "\n")
    for genome in genomes.values():
        row = [
            genome.id,
            genome.length,
            genome.prophage,
            genome.proviral_length,
            genome.total_genes,
            genome.viral_genes,
            genome.host_genes,
            genome.quality,
            genome.miuvig,
            genome.completeness,
            genome.method,
            genome.contamination,
            genome.kmer_freq,
            "; ".join(genome.warnings),
        ]
        out.write("\t".join([str(_) for _ in row]) + "\n")

    logger.info("Run time: %s seconds" % round(time.time() - program_start, 2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
