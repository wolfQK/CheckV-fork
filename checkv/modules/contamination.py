import argparse
import csv
import os
import shutil
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
    parser.set_defaults(program="contamination")
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
        help="Number of threads to use for Prodigal and hmmsearch",
    )
    parser.add_argument(
        "--restart",
        action="store_true",
        default=False,
        help="Overwrite existing intermediate files. By default CheckV continues where program left off",
    )
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Suppress logging messages",
    )
    parser.add_argument(
        "--exclude", type=str, required=False, metavar="PATH", help=argparse.SUPPRESS,
    )


def set_defaults(args):
    key_values = [("exclude", None)]
    for key, value in key_values:
        if key not in args:
            args[key] = value


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
        genome.count_viral = sum(genes[_].cat == 1 for _ in genome.genes)
        genome.count_host = sum(genes[_].cat == -1 for _ in genome.genes)


def compute_delta(my_genes, s1, e1, s2, e2, gc_weight):

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

    return d


def define_regions(
    genome,
    genes,
    gc_weight,
    delta_cutoff,
    min_host_fract,
    min_host_genes,
    min_viral_genes,
):
    """
    1. Score each possible breakpoint using combination of viral annotations and GC content
    2. Identify breakpoints. See code below for list of rules for this process
    3. Identify host/viral regions based on breakpoints
    4. Return list of breakpoints and regions
    """

    # fetch genes
    my_genes = [genes[_] for _ in genome.genes]

    # determine window size
    # 30% of genes with 15 gene min and 50 gene max
    win_size = min([max([15, int(round(0.30 * len(my_genes)))]), 50])

    # identify breakpoints
    breaks = []
    while True:

        # determine s1
        if len(breaks) == 0:
            s1 = 0
        else:
            s1 = breaks[-1]["coords"][-2]

        # determine window coords (gene indexes)
        coords = []
        for i in range(1, len(my_genes)):
            e1 = i
            s2, e2 = i, min([i + win_size, len(my_genes)])
            size1 = e1 - s1
            size2 = e2 - s2
            edge1 = True if s1 == 0 else False
            edge2 = True if e2 == len(my_genes) else False
            # windows must be 40 genes unless they start/end at contig boundary
            if size1 < win_size and not edge1 or size2 < win_size and not edge2:
                continue
            coords.append([s1, e1, s2, e2])

        # score each possible breakpoint
        deltas = []
        for s1, e1, s2, e2 in coords:
            d = compute_delta(my_genes, s1, e1, s2, e2, gc_weight)
            deltas.append(d)

        # filter each possible breakpoint
        filtered = []
        for d in deltas:
            # delta not significant
            if abs(d["delta"]) < delta_cutoff:
                continue
            # always require at least 1 annotated gene per region
            if d["delta"] < 0 and (d["v1"].count(-1) == 0 or d["v2"].count(1) == 0):
                continue
            if d["delta"] > 0 and (d["v1"].count(1) == 0 or d["v2"].count(-1) == 0):
                continue
            # require at least N host genes per region for contigs with > 10 genes
            if (
                d["delta"] < 0
                and d["v1"].count(-1) < min_host_genes
                and len(my_genes) > 10
            ):
                continue
            if (
                d["delta"] > 0
                and d["v2"].count(-1) < min_host_genes
                and len(my_genes) > 10
            ):
                continue
            # require at least N viral genes per region for contigs with > 10 genes
            if (
                d["delta"] > 0
                and d["v1"].count(1) < min_viral_genes
                and len(my_genes) > 10
            ):
                continue
            if (
                d["delta"] < 0
                and d["v2"].count(1) < min_viral_genes
                and len(my_genes) > 10
            ):
                continue
            # require minimum % of host genes in either window
            if (
                d["win1_fract_host"] < min_host_fract
                and d["win2_fract_host"] < min_host_fract
            ):
                continue
            # store
            filtered.append(d)

        # select breakpoint
        selected = None
        for d in filtered:
            # add first breakpoint
            if selected is None:
                selected = d
            # breakpoint in same orientation --> update
            elif np.sign(d["delta"]) == np.sign(selected["delta"]):
                if abs(d["delta"]) > abs(selected["delta"]):
                    selected = d
            # breakpoint in opposite orientation
            else:
                break

        # if breakpoint selected, add it to list, look for next breakpoint
        if selected is None:
            break
        else:
            breaks.append(selected)

    # update last break so end coord is contig end
    if len(breaks) > 0:
        breaks[-1]["coords"][-1] = len(my_genes)

    # define regions based on breakpoints
    regions = []
    for b in breaks:
        s1, e1, s2, e2 = b["coords"]
        d = compute_delta(my_genes, s1, e1, s2, e2, gc_weight)
        region = {}
        region["delta"] = d["delta"]
        region["type"] = "host" if d["delta"] < 0 else "viral"
        region["start_gene"] = s1  # 0-indexed
        region["end_gene"] = e1  # 0-indexed
        region["start_pos"] = (
            regions[-1]["end_pos"] + 1 if len(regions) > 0 else 1
        )  # 1-indexed
        region["end_pos"] = my_genes[e1 - 1].end  # 1-indexed
        region["size"] = region["end_gene"] - region["start_gene"]
        region["length"] = region["end_pos"] - region["start_pos"] + 1
        region["host_genes"] = len([_ for _ in d["v1"] if _ == -1])
        region["viral_genes"] = len([_ for _ in d["v1"] if _ == 1])
        regions.append(region)

    # handle last region
    if len(regions) > 0:
        s1, e1 = regions[-1]["start_gene"], regions[-1]["end_gene"]
        s2, e2 = e1, len(genome.genes)
        d = compute_delta(my_genes, s1, e1, s2, e2, gc_weight)
        region = {}
        region["type"] = "viral" if regions[-1]["type"] == "host" else "host"
        region["start_gene"] = s2  # 0-indexed
        region["end_gene"] = e2  # 0-indexed
        region["start_pos"] = regions[-1]["end_pos"] + 1  # 1-indexed
        region["end_pos"] = genome.length  # 1-indexed
        region["size"] = region["end_gene"] - region["start_gene"]
        region["length"] = region["end_pos"] - region["start_pos"] + 1
        region["host_genes"] = len([_ for _ in d["v2"] if _ == -1])
        region["viral_genes"] = len([_ for _ in d["v2"] if _ == 1])
        regions.append(region)

    return regions


def main(args):

    program_start = time.time()
    set_defaults(args)
    logger = utility.get_logger(args["quiet"])
    utility.check_executables(["prodigal", "hmmsearch"])
    args["db"] = utility.check_database(args["db"])
    args["tmp"] = os.path.join(args["output"], "tmp")
    if not os.path.exists(args["output"]):
        os.makedirs(args["output"])
    if args["restart"] and os.path.exists(args["tmp"]):
        shutil.rmtree(args["tmp"])
    if not os.path.exists(args["tmp"]):
        os.makedirs(args["tmp"])

    logger.info(f"\nCheckV v{checkv.__version__}: contamination")

    logger.info("[1/8] Reading database info...")
    hmm_info = {}
    p = os.path.join(args["db"], "hmm_db/checkv_hmms.tsv")
    for r in csv.DictReader(open(p), delimiter="\t"):
        r["score_cutoff"] = max(
            [float(r["score_cutoff"]), 25.0]
        )  # min cutoff == 25 bits
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
        genome.regions = None
        genomes[genome.id] = genome

    args["faa"] = os.path.join(args["tmp"], "proteins.faa")
    if not os.path.exists(args["faa"]):
        logger.info("[3/8] Calling genes with Prodigal...")
        utility.call_genes(args["input"], args["output"], args["threads"])
    else:
        logger.info("[3/8] Skipping gene calling...")

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

    args["hmmout"] = os.path.join(args["tmp"], "hmmsearch.txt")
    if not os.path.exists(args["hmmout"]):
        logger.info("[5/8] Running hmmsearch...")
        db_dir = os.path.join(args["db"], "hmm_db/checkv_hmms")
        utility.search_hmms(args["tmp"], args["threads"], db_dir)
    else:
        logger.info("[5/8] Skipping hmmsearch...")

    logger.info("[6/8] Annotating genes...")
    annotate_genes(hmm_info, genomes, genes, args)

    logger.info("[7/8] Identifying host regions...")
    for genome in genomes.values():
        genome.regions = define_regions(
            genome,
            genes,
            min_host_genes=2,
            min_viral_genes=2,
            min_host_fract=0.30,
            gc_weight=0.02,
            delta_cutoff=1.2,
        )

    logger.info("[8/8] Writing results...")

    with open(os.path.join(args["output"], "viruses.fna"), "w") as out:
        for genome in genomes.values():
            if len(genome.regions) == 0:
                out.write(">" + genome.id + "\n" + genome.seq + "\n")

    with open(os.path.join(args["output"], "proviruses.fna"), "w") as out:
        for genome in genomes.values():
            if len(genome.regions) > 0:
                viral_regions = [r for r in genome.regions if r["type"] == "viral"]
                if len(viral_regions) == 0:
                    viral_regions = genome.regions
                for i, r in enumerate(viral_regions):
                    header = (
                        genome.id
                        + "_"
                        + str(i + 1)
                        + " "
                        + str(r["start_pos"])
                        + "-"
                        + str(r["end_pos"])
                        + "/"
                        + str(genome.length)
                    )
                    seq = genome.seq[r["start_pos"] - 1 : r["end_pos"]]
                    out.write(">" + header + "\n" + seq + "\n")

    with open(os.path.join(args["output"], "contamination.tsv"), "w") as out:
        header = ["contig_id", "contig_length"]
        header += ["total_genes", "viral_genes", "host_genes"]
        header += ["provirus", "proviral_length", "host_length"]
        header += [
            "region_types",
            "region_lengths",
            "region_coords_bp",
            "region_coords_genes",
        ]
        header += ["region_viral_genes", "region_host_genes"]
        out.write("\t".join(header) + "\n")
        for genome in genomes.values():
            row = [
                genome.id,
                genome.length,
                len(genome.genes),
                genome.count_viral,
                genome.count_host,
            ]
            if len(genome.regions) > 0:
                viral_length = sum(
                    r["length"] for r in genome.regions if r["type"] == "viral"
                )
                host_length = sum(
                    r["length"] for r in genome.regions if r["type"] == "host"
                )
                region_types = ",".join([r["type"] for r in genome.regions])
                region_lengths = ",".join([str(r["length"]) for r in genome.regions])
                region_coords_bp = ",".join(
                    [
                        str(r["start_pos"]) + "-" + str(r["end_pos"])
                        for r in genome.regions
                    ]
                )
                region_coords_genes = ",".join(
                    [
                        str(r["start_gene"] + 1) + "-" + str(r["end_gene"])
                        for r in genome.regions
                    ]
                )
                region_viral_genes = ",".join(
                    [str(r["viral_genes"]) for r in genome.regions]
                )
                region_host_genes = ",".join(
                    [str(r["host_genes"]) for r in genome.regions]
                )
                row += ["Yes", viral_length, host_length]
                row += [
                    region_types,
                    region_lengths,
                    region_coords_bp,
                    region_coords_genes,
                ]
                row += [region_viral_genes, region_host_genes]
            else:
                row += ["No"] + ["NA"] * 8
            out.write("\t".join([str(_) for _ in row]) + "\n")

    with open(os.path.join(args["tmp"], "gene_features.tsv"), "w") as out:
        header = [
            "contig_id",
            "gene_num",
            "start",
            "end",
            "strand",
            "gc",
            "hmm_cat",
            "hmm_db",
            "hmm_name",
            "evalue",
            "score",
        ]
        out.write("\t".join(header) + "\n")
        for genome in genomes.values():
            for gene_id in genome.genes:
                gene = genes[gene_id]
                row = [
                    genome.id,
                    gene_id.split("_")[-1],
                    gene.start,
                    gene.end,
                    gene.strand,
                    round(gene.gc, 1),
                    gene.cat,
                ]
                if gene.hmm_hit:
                    row.append(hmm_info[gene.hmm_hit["tname"]]["database"])
                    row.append(gene.hmm_hit["tname"])
                    row.append(gene.hmm_hit["eval"])
                    row.append(gene.hmm_hit["score"])
                else:
                    row.append("NA")
                    row.append("NA")
                    row.append("NA")
                    row.append("NA")
                out.write("\t".join([str(_) for _ in row]) + "\n")

    logger.info("Run time: %s seconds" % round(time.time() - program_start, 2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
