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

    path = os.path.join(args["output"], "completeness.tsv")
    if not os.path.exists(path): sys.exit(f"Error: input file does not exist: {path}\n")
    path = os.path.join(args["output"], "repeats.tsv")
    if not os.path.exists(path): sys.exit(f"Error: input file does not exist: {path}\n")

    logger.info(f"CheckV version: {checkv.__version__}")
    logger.info("")

    logger.info("[1/6] Reading input sequences...")
    genomes = {}
    for header, seq in utility.read_fasta(args["input"]):
        genome = Genome()
        genome.id = header.split()[0]
        genome.seq = seq
        genome.length = len(seq)
        genome.copies = "NA"
        genome.termini = "NA"
        genome.contamination = "NA"
        genome.prophage = "NA"
        genome.total_genes = "NA"
        genome.viral_genes = "NA"
        genome.host_genes = "NA"
        genome.completeness = "NA"
        genome.method = "NA"
        genomes[genome.id] = genome

    p = os.path.join(args["output"], "contamination.tsv")
    if os.path.exists(p):
        logger.info("[2/6] Reading results from contamination module...")
        for r in csv.DictReader(open(p), delimiter="\t"):
            genome = genomes[r["contig_id"]]
            
            # num genes
            genome.total_genes = r["total_genes"]
            genome.viral_genes = r["viral_genes"]
            genome.host_genes = r["host_genes"]
            
            # complete prophage
            if r["region_types"] == "host,viral,host":
                genome.termini = "complete-prophage"
        
            # partial prophage
            if "viral" in r["region_types"] and "host" in r["region_types"]:
                genome.prophage = "Yes"
                genome.contamination = round(
                    100.0 * int(r["host_length"]) / genome.length, 2
                )
            
            # entirely viral
            elif "viral" in r["region_types"] and "host" not in r["region_types"]:
                genome.prophage = "No"
                genome.contamination = 0.0

            # no viral region
            # contamination could be 100% or 0%
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
            if "hmm_completness" in r: r["hmm_completeness"] = r["hmm_completness"] ## temp fix
            
            # Med/High-confidence AAI-based estimate
            if r["aai_confidence"] in ["high", "medium"]:
                genome.completeness = min([round(float(r["aai_completeness"]), 2), 100.0])
                genome.method = "AAI-based"

            # HMM-based estimate
            elif r["hmm_completeness"] != "NA":
                genome.completeness = min([round(float(r["hmm_completeness"]), 2), 100.0])
                genome.method = "HMM-based"

    logger.info("[4/6] Reading results from repeats module...")
    p = os.path.join(args["output"], "repeats.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            if "contig_id" not in r: r["contig_id"] = r["genome_id"]
            genome = genomes[r["contig_id"]]
            genome.copies = r["genome_copies"]
            if r["repeat_type"] != "NA" and r["repeat_flagged"] == "No":
                genome.termini = r["repeat_length"] + "-bp-" + r["repeat_type"]

    logger.info("[5/6] Classifying contigs into quality tiers...")
    for genome in genomes.values():
        # DTR == complete if completeness > 90% or no AAI-based estimate
        if "DTR" in genome.termini:
            if genome.completeness == "NA" or genome.method == "HMM-based":
                genome.quality = "Complete"
                genome.miuvig = "High-quality"
            elif genome.completeness < 50 and genome.method == "AAI-based":
                genome.quality = "Low-quality"
                genome.miuvig = "Genome-fragment"
            elif genome.completeness < 90 and genome.method == "AAI-based":
                genome.quality = "Medium-quality"
                genome.miuvig = "Genome-fragment"
            else:
                genome.quality = "Complete"
                genome.miuvig = "High-quality"
        # ITR/prophage == complete if completeness > 90%
        elif "ITR" in genome.termini or "complete-prophage" in genome.termini:
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
        # high quality == completeness > 90%
        elif genome.completeness != "NA" and genome.completeness >= 90:
            genome.quality = "High-quality"
            genome.miuvig = "High-quality"
        # medium quality == completeness > 50%
        elif genome.completeness != "NA" and genome.completeness >= 50:
            genome.quality = "Medium-quality"
            genome.miuvig = "Genome-fragment"
        # low quality == completeness < 50%
        elif genome.completeness != "NA":
            genome.quality = "Low-quality"
            genome.miuvig = "Genome-fragment"
        # everything else is not determined
        else:
            genome.quality = "Not-determined"
            genome.miuvig = "Genome-fragment"

    logger.info("[6/6] Writing results...")
    header = [
        "contig_id",
        "contig_length",
        "genome_copies",
        "gene_count",
        "viral_genes",
        "host_genes",
        "checkv_quality",
        "miuvig_quality",
        "completeness",
        "completeness_method",
        "contamination",
        "prophage",
        "termini",
    ]
    out = open(args["output"] + "/quality_summary.tsv", "w")
    out.write("\t".join(header) + "\n")
    for genome in genomes.values():
        row = [
            genome.id,
            genome.length,
            genome.copies,
            genome.total_genes,
            genome.viral_genes,
            genome.host_genes,
            genome.quality,
            genome.miuvig,
            genome.completeness,
            genome.method,
            genome.contamination,
            genome.prophage,

            genome.termini,
        ]
        out.write("\t".join([str(_) for _ in row]) + "\n")

    logger.info("\nDone!")
    logger.info("Run time: %s seconds" % round(time.time() - program_start, 2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
