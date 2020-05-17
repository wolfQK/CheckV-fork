import logging
import os
import shutil
import sys
import time
import checkv
import subprocess as sp
import Bio.SeqIO
from checkv import utility

def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="update_database")
    parser.add_argument(
        "source_db",
        type=str,
        help="Path to current CheckV database.",
    )
    parser.add_argument(
        "dest_db",
        type=str,
        help="Path to updated CheckV database.",
    )
    parser.add_argument(
        "genomes",
        type=str,
        help="FASTA file of complete genomes to add to database, where each nucleotide sequence represents one genome.",
    )
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Suppress logging messages",
    )
    parser.add_argument(
        "--restart", action="store_true", default=False, help="Overwrite existing database",
    )
    parser.add_argument(
        "--threads", metavar="INT", type=int, default=1, help="Number of threads for Prodigal and DIAMOND",
    )


def main(args):
    program_start = time.time()
    logger = utility.get_logger(args["quiet"])

    if os.path.exists(args["dest_db"]):
        if args["restart"]:
            shutil.rmtree(args["dest_db"])
        else:
            sys.exit("Error: database already exists")

    logger.info(f"\nCheckV v{checkv.__version__}: update_database")

    logger.info("[1/5] Copying database to new destination...")
    shutil.copytree(args["source_db"], args["dest_db"])

    logger.info("[2/5] Calling genes on new genomes with Prodigal...")
    tmp = os.path.join(args["dest_db"], "tmp")
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    utility.call_genes(args["genomes"], args["dest_db"], args["threads"])

    logger.info("[3/5] Updating FASTA databases...")
    faa = os.path.join(args["dest_db"], "genome_db/checkv_reps.faa")
    with open(faa, "a") as out:
        for line in open(os.path.join(args["dest_db"], "tmp/proteins.faa")):
            out.write(line)
    fna = os.path.join(args["dest_db"], "genome_db/checkv_reps.fna")
    with open(fna, "a") as out:
        for line in open(args["genomes"]):
            out.write(line)

    logger.info("[4/5] Updating DIAMOND database...")
    dmnd = os.path.join(args["dest_db"], "genome_db/checkv_reps.dmnd")
    os.remove(dmnd)
    cmd = f"diamond makedb --in {faa} --db {dmnd} --threads {args['threads']} &> {dmnd}.log"
    p = sp.Popen(cmd, shell=True)
    p.wait()

    logger.info("[5/5] Updating genome info...")
    tsv = os.path.join(args["dest_db"], "genome_db/checkv_reps.tsv")
    with open(tsv, "a") as f:
        for r in Bio.SeqIO.parse(fna, "fasta"):
            f.write("\t".join([r.id, "circular", str(len(r.seq))])+"\n")

    logger.info("Run time: %s seconds" % round(time.time() - program_start, 2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
    shutil.rmtree(tmp)

