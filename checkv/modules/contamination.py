#!/usr/bin/env python

import argparse
import csv
import logging
import os
import subprocess as sp
import time
import shutil

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
    parser.set_defaults(program="contamination")
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
        help="Number of threads to use for Prodigal and hmmsearch",
    )
    parser.add_argument(
        "--restart",
        action="store_true",
        default=False,
        help="Overwrite existing intermediate files. By default CheckV continues where program left off",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Display logging messages",
    )
    parser.add_argument(
        "--exclude",
        type=str,
        required=False,
        metavar="PATH",
        help=argparse.SUPPRESS,
    )

def compute_gc(x):
    return 100.0 * (x.count("G") + x.count("C")) / len(x)


def annotate_genes(hmm_info, genomes, genes, args):
    """
    1. read hmm bitscore cutoffs for viral/microbial annotation
    2. assign genes to hmms according to their best hit using cutoffs from (1)
    3. Assign gene categories (1 for viral, -1 for microbial)
    4. Summarize the number of viral and microbial annotations per genome
    """

    # 1. identify best hit at cutoffs
    for r in utility.parse_hmmsearch(args["hmmout"]):
        r["tname"] = r["tname"] if r["tname"] in hmm_info else r["tacc"]
        gene = genes[r["qname"]]
        if r["tname"] not in hmm_info:
            continue
        elif r["score"] < float(hmm_info[r["tname"]]["score_cutoff"]):
            continue
        elif gene.hmm_hit is None or r["score"] > gene.hmm_hit["score"]:
            gene.hmm_hit = r

    # 2. annotate genes based on best hit at selected cutoff
    to_cat = {"viral": 1, "microbial": -1}
    for genome in genomes.values():
        for gene_id in genome.genes:
            hmm_hit = genes[gene_id].hmm_hit
            if hmm_hit is not None:
                hmm_type = hmm_info[hmm_hit["tname"]]["category"]
                genes[gene_id].cat = to_cat[hmm_type]

    # 3. summarize hits
    for genome in genomes.values():
        genome.count_viral = sum(1 for _ in genome.genes if genes[_].cat == 1)
        genome.count_host = sum(1 for _ in genome.genes if genes[_].cat == -1)


def define_regions(
    genome, genes, min_fract, min_genes=10, max_genes=50, gc_weight=0.02, delta_cutoff=1.2
):
    """
    1. Determine win size based on <min_genes>, <max_genes>, <min_fract>
    2. Score each possible breakpoint using combination of viral annotations and GC content
    3. Identify breakpoints. See code below for list of rules for this process
    4. Identify host/viral regions based on breakpoints
    5. Return list of breakpoints and regions
    """

    # 1. determine window size
    # window size adjusted based on contig length
    # at least min_genes or min_fract of genes (whichever is bigger)
    # but no more than max_genes
    count_fract = int(round(len(genome.genes) * min_fract))
    win_size = min([max([count_fract, min_genes]), max_genes])

    # 2. Score each possible breakpoint

    # get gene annotations and gc content
    # microbial annotation = -1
    # viral annotation = 1
    my_genes = [genes[_] for _ in genome.genes]

    # compute deltas
    deltas = []
    for i in range(1, len(my_genes)):

        # get gene indexes for 2 windows
        s1, e1 = max([i - win_size, 0]), i
        s2, e2 = i, min([i + win_size, len(my_genes)])

        # extend windows to ensure at least 1 annotated gene
        if all([g.cat == 0 for g in my_genes[s1:e1]]):
            for j in range(0, s1)[::-1]:
                s1 = j
                if my_genes[s1].cat != 0:
                    break
        if all([g.cat == 0 for g in my_genes[s2:e2]]):
            for j in range(e2 + 1, len(my_genes)):
                e2 = j
                if my_genes[e2].cat != 0:
                    break

        # get gene values for 2 windows
        win1 = my_genes[s1:e1]
        win2 = my_genes[s2:e2]
        v1 = [g.cat for g in win1 if g.cat != 0]
        v2 = [g.cat for g in win2 if g.cat != 0]
        g1 = [g.gc for g in win1]
        g2 = [g.gc for g in win2]

        # compute delta between windows
        if len(v1) > 0 and len(v2) > 0:
            delta_v = np.average(v1) - np.average(v2)
            delta_g = abs(np.average(g1) - np.average(g2))
            delta = (abs(delta_v) + delta_g * gc_weight) * np.sign(delta_v)
        else:
            delta = 0

        # store
        d = {
            "delta": delta,
            "coords": [s1, e1, s2, e2],
            "v1": v1,
            "v2": v2,
            "v1_len": len(v1),
            "v2_len": len(v2),
            "win1_len": len(win1),
            "win2_len": len(win2),
            "win1_fract_host": 1.0 * len([_ for _ in v1 if _ == -1]) / len(win1),
            "win2_fract_host": 1.0 * len([_ for _ in v2 if _ == -1]) / len(win2),
        }
        deltas.append(d)

    # 3. Identify breakpoints
    breaks = []
    for d in deltas:

        # filter breakpoint
        # delta not significant
        if abs(d["delta"]) < delta_cutoff:
            continue
        # at least 4 total annotated genes
        elif d["v1_len"] + d["v2_len"] < 4:
            continue
        # <20% host genes in both windows
        elif d["win1_fract_host"] < 0.20 and d["win2_fract_host"] < 0.20:
            continue

        # add first breakpoint
        if len(breaks) == 0:
            breaks.append(d)

        # deltas in opposite directions
        elif np.sign(d["delta"]) != np.sign(breaks[-1]["delta"]):

            # breakpoints both start at L edge
            # use breakpoint that is further R
            if d["coords"][0] == 0 and breaks[-1]["coords"][0] == 0:
                breaks[-1] = d
            # breakpoints both start at R edge
            # use breakpoint that is further L
            elif d["coords"][-1] == len(genome.genes) and breaks[-1]["coords"][-1] == len(
                genome.genes
            ):
                pass
            # does not overlap: add new breakpoint
            else:
                breaks.append(d)

        # deltas in same direction: update
        elif np.sign(d["delta"]) == np.sign(breaks[-1]["delta"]):
            # larger delta
            if abs(d["delta"]) > abs(breaks[-1]["delta"]):
                breaks[-1] = d
            # same delta, negative --> host-virus --> use left-most boundary
            elif abs(d["delta"]) == abs(breaks[-1]["delta"]) and d["delta"] < 0:
                pass
            # same delta, positive, virus-host --> use right-most boundary
            elif abs(d["delta"]) == abs(breaks[-1]["delta"]) and d["delta"] > 0:
                breaks[-1] = d

    # 4. identify viral/host regions
    regions = []
    for d in breaks:
        s1, e1, s2, e2 = d["coords"]
        region = {
            "type": "host" if d["delta"] < 0 else "viral",
            "delta": d["delta"],
            "start_pos": regions[-1]["end_pos"] + 1
            if len(regions) > 0
            else 1,  # 1-indexed
            "end_pos": my_genes[e1 - 1].end,  # 1-indexed
            "start_gene": regions[-1]["end_gene"] if len(regions) > 0 else 0,  # 0-indexed
            "end_gene": e1,  # 0-indexed
        }
        region["size"] = region["end_gene"] - region["start_gene"]
        region["length"] = region["end_pos"] - region["start_pos"] + 1
        regions.append(region)

    # handle last region
    region = {
        "start_pos": regions[-1]["end_pos"] + 1 if len(regions) > 0 else 1,
        "end_pos": genome.length,
        "start_gene": regions[-1]["end_gene"] if len(regions) > 0 else 0,
        "end_gene": len(genome.genes),
    }
    region["size"] = region["end_gene"] - region["start_gene"]
    region["length"] = region["end_pos"] - region["start_pos"] + 1

    # determine if host/viral/unclassified
    # no breaks detected; determine if region is viral, host, or unclassified
    if len(regions) == 0:
        if genome.count_viral > 0 and genome.count_viral >= genome.count_host:
            region["type"] = "viral"
        elif genome.count_host > genome.count_viral:
            region["type"] = "host"
        elif genome.count_viral == 0:
            region["type"] = "unclassified"
    # last delta was postive, indicating a viral-host breakpoint
    elif regions[-1]["delta"] > 0:
        region["type"] = "host"
    # last delta was negative, indicating a host-virus breakpoint
    else:
        region["type"] = "viral"
    regions.append(region)

    return regions


def main(args):

    program_start = time.time()
    logger = utility.get_logger(args["verbose"])
    utility.check_executables(["prodigal", "hmmsearch"])
    args["db"] = utility.check_database(args["db"])
    args["tmp"] = os.path.join(args["output"], "tmp")
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
    if args["restart"] and os.path.exists(args["tmp"]):
        shutil.rmtree(args["tmp"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    logger.info("[1/8] Reading database info...")
    hmm_info = {}
    p = os.path.join(args["db"], "checkv_hmms.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        r["score_cutoff"] = max([float(r["score_cutoff"]), 25.0]) # min cutoff == 25 bits
        hmm_info[r["hmm"]] = r
    if args["exclude"]:
        for l in open(args["exclude"]):
            del hmm_info[l.rstrip()]

    logger.info("[2/8] Reading genome info...")
    genomes = {}
    for r in utility.read_fasta(args["input"]):
        header, seq = r
        genome = Genome()
        genome.id = header.split()[0]
        genome.length = len(seq)
        genome.genes = []
        genome.seq = seq
        genome.viral_hits = {}
        genomes[genome.id] = genome

    logger.info("[3/8] Calling genes with prodigal...")
    args["faa"] = os.path.join(args["tmp"], "proteins.faa")
    if not os.path.exists(args["faa"]):
        utility.call_genes(args["input"], args["output"], args["threads"])

    logger.info("[4/8] Reading gene info...")  # assumes PRODIGAL V2.6.3 FORMAT
    genes = {}
    for header, seq in utility.read_fasta(args["faa"]):
        gene = Gene()
        gene.id = header.split()[0]
        gene.start = int(header.split()[2])
        gene.end = int(header.split()[4])
        gene.strand = int(header.split()[6])
        gene.genome_id = gene.id.rsplit("_", 1)[0]
        gene.gc = compute_gc(genomes[gene.genome_id].seq[gene.start - 1 : gene.end])
        gene.cat = 0
        gene.hmm_hit = None
        genes[gene.id] = gene
        genomes[gene.genome_id].genes.append(gene.id)

    logger.info("[5/8] Running hmmsearch...")
    args["hmmout"] = os.path.join(args["tmp"], "hmmsearch.txt")
    if not os.path.exists(args["hmmout"]):
        db_dir = os.path.join(args["db"], "checkv_hmms")
        utility.search_hmms(args["tmp"], args["threads"], db_dir)

    logger.info("[6/8] Annotating genes...")
    annotate_genes(hmm_info, genomes, genes, args)

    logger.info("[7/8] Identifying host regions...")
    for genome in genomes.values():
        genome.regions = define_regions(genome, genes, min_fract=1.0, max_genes=35)

    logger.info("[8/8] Writing results...")
    out = open(args["output"] + "/contamination.tsv", "w")
    header = ["contig_id", "contig_length", "viral_length", "host_length"]
    header += ["total_genes", "viral_genes", "host_genes"]
    header += ["region_types", "region_lengths", "region_coords", "region_genes"]
    out.write("\t".join(header) + "\n")
    for genome in genomes.values():
        viral_length = sum(r["length"] for r in genome.regions if r["type"] == "viral")
        host_length = sum(r["length"] for r in genome.regions if r["type"] == "host")
        region_types = ",".join([r["type"] for r in genome.regions])
        region_lengths = ",".join([str(r["length"]) for r in genome.regions])
        region_coords = ",".join([str(r["start_pos"])+"-"+str(r["end_pos"]) for r in genome.regions])
        region_genes = ",".join([str(r["start_gene"]+1)+"-"+str(r["end_gene"]) for r in genome.regions])
        row = [genome.id, genome.length, viral_length, host_length]
        row += [len(genome.genes), genome.count_viral, genome.count_host]
        row += [region_types, region_lengths, region_coords, region_genes]
        out.write("\t".join([str(_) for _ in row]) + "\n")

    with open(args["tmp"]+"/gene_features.tsv", "w") as out:
        header = ["contig_id", "gene_num", "start", "end", "strand", "hmm_cat", "gc"]
        out.write("\t".join(header)+"\n")
        for genome in genomes.values():
            for gene_id in genome.genes:
                gene = genes[gene_id]
                num = gene_id.split("_")[-1]
                row = [genome.id, num, gene.start, gene.end, gene.strand, gene.cat, round(gene.gc,1)]
                out.write("\t".join([str(_) for _ in row])+"\n")

    with open(args["tmp"]+"/gene_annotations.tsv", "w") as out:
        header = ["contig_id", "gene_num", "target_hmm", "hmm_db", "hmm_cat", "target_score", "target_evalue"]
        out.write("\t".join(header)+"\n")
        for genome in genomes.values():
            for gene_id in genome.genes:
                gene = genes[gene_id]
                if gene.hmm_hit:
                    num = gene_id.split("_")[-1]
                    hmm = gene.hmm_hit["tname"]
                    db = hmm_info[hmm]["database"]
                    cat = hmm_info[hmm]["category"]
                    score = gene.hmm_hit["score"]
                    eval = gene.hmm_hit["eval"]
                    row = [genome.id, num, hmm, db, cat, score, eval]
                    out.write("\t".join([str(_) for _ in row])+"\n")

    logger.info("\nDone!")
    logger.info("Run time: %s seconds" % round(time.time()-program_start,2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(),2))



