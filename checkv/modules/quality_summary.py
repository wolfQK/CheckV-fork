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
        pass


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="quality_summary")
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
        "--quiet",
        action="store_true",
        default=False,
        help="Suppress logging messages",
    )


def main(args):

    program_start = time.time()
    logger = utility.get_logger(args["quiet"])
    args["tmp"] = os.path.join(args["output"], "tmp")
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    logger.info("[1/6] Reading input sequences...")
    genomes = {}
    for r in Bio.SeqIO.parse(args["input"], "fasta"):
        genome = Genome()
        genome.id = r.id
        genome.seq = str(r.seq).upper()
        genome.length = len(genome.seq)
        genome.termini = "NA"
        genomes[genome.id] = genome

    logger.info("[2/6] Reading results from contamination module...")
    p = os.path.join(args["output"], "contamination.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            if "contig_id" not in r: r["contig_id"] = r["genome_id"]
            genome = genomes[r["contig_id"]]
            
            if r["region_types"] == "host,viral,host":
                genome.termini = "complete-prophage"
            
            if "viral" in r["region_types"]:
                if "host" in r["region_types"]:
                    genome.prophage = "Yes"
                    genome.contamination = round(100.0 * int(r["host_length"]) / genome.length,2)
                else:
                    genome.prophage = "No"
                    genome.contamination = 0.0
            else:
                genome.prophage = "No"
                genome.contamination = "NA"

    logger.info("[3/6] Reading results from completeness module...")
    p = os.path.join(args["output"], "completeness.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            if "contig_id" not in r: r["contig_id"] = r["genome_id"]
            genome = genomes[r["contig_id"]]
            if r["confidence"] in ["high", "medium"]:
                genome.completeness = round(float(r["completeness"]),2)
                if genome.completeness > 110:
                    genome.comment = "'Warning: completeness >110%. Estimate may be unreliable.'"
                else:
                    genome.comment = ""
            else:
                genome.completeness = "NA"
                genome.comment = ""

    logger.info("[4/6] Reading results from terminal repeats module...")
    p = os.path.join(args["output"], "terminal_repeats.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            if "contig_id" not in r: r["contig_id"] = r["genome_id"]
            genome = genomes[r["contig_id"]]
            genome.termini = r["repeat_length"] + "-bp-" + r["repeat_type"]

    logger.info("[5/6] Classifying contigs into quality tiers...")
    for genome in genomes.values():
        if "DTR" in genome.termini:
            if genome.completeness == "NA":
                genome.quality = "Complete"
                genome.miuvig = "High-quality"
            elif genome.completeness < 50:
                genome.quality = "Low-quality"
                genome.miuvig = "Genome-fragment"
            elif genome.completeness < 90:
                genome.quality = "Medium-quality"
                genome.miuvig = "Genome-fragment"
            else:
                genome.quality = "Complete"
                genome.miuvig = "High-quality"
        elif ("ITR" in genome.termini
                or "complete-prophage" in genome.termini):
            if genome.completeness == "NA":
                genome.quality = "Not-determined"
                genome.miuvig = "Genome-fragment"
            elif genome.completeness < 50:
                genome.quality = "Low-quality"
                genome.miuvig = "Genome-fragment"
            elif genome.completeness < 90:
                genome.quality = "Medium-quality"
                genome.miuvig = "Genome-fragment"
            else:
                genome.quality = "Complete"
                genome.miuvig = "High-quality"
        elif genome.completeness != "NA" and genome.completeness >= 90:
            genome.quality = "High-quality"
            genome.miuvig = "High-quality"
        elif genome.completeness != "NA" and genome.completeness >= 50:
            genome.quality = "Medium-quality"
            genome.miuvig = "Genome-fragment"
        elif genome.completeness != "NA":
            genome.quality = "Low-quality"
            genome.miuvig = "Genome-fragment"
        else:
            genome.quality = "Not-determined"
            genome.miuvig = "Genome-fragment"

    logger.info("[6/6] Writing results...")
    header = [
        "contig_id",
        "contig_length",
        "checkv_quality",
        "miuvig_quality",
        "completeness",
        "contamination",
        "prophage",
        "termini",
        "comments"
    ]
    out = open(args["output"] + "/quality_summary.tsv", "w")
    out.write("\t".join(header) + "\n")
    for genome in genomes.values():
        row = [
            genome.id,
            genome.length,
            genome.quality,
            genome.miuvig,
            genome.completeness,
            genome.contamination,
            genome.prophage,
            genome.termini,
            genome.comment
        ]
        out.write("\t".join([str(_) for _ in row]) + "\n")

    logger.info("\nDone!")
    logger.info("Run time: %s seconds" % round(time.time()-program_start,2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(),2))




