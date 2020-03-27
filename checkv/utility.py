import logging
import math
import multiprocessing as mp
import os
import shutil
import signal
import subprocess as sp
import sys
import time

import psutil
from Bio import SeqIO


def get_logger(verbosity):
    if verbosity:
        logging.basicConfig(level=logging.INFO, format='%(message)s')
    else:
        logging.basicConfig(level=logging.WARNING, format='%(message)s')
    return logging.getLogger()


def check_executables(requirements):
    fails = 0
    for program in requirements:
        found = False
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path.strip('"'), program)
            if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
                found = True
                break
        if not found:
            msg = f"Error: required program '{program}' not executable or not found on $PATH\n"
            sys.stderr.write(msg)
            fails += 1
    if fails > 0:
        sys.exit()


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def killtree(pid, including_parent=True):
    parent = psutil.Process(pid)
    for child in parent.children(recursive=True):
        child.kill()
    if including_parent:
        parent.kill()


def parallel(function, argument_list, threads):
    """ Based on: https://gist.github.com/admackin/003dd646e5fadee8b8d6 """
    threads = len(argument_list)
    pool = mp.Pool(threads, init_worker)
    try:
        results = []
        for arguments in argument_list:
            p = pool.apply_async(function, args=arguments)
            results.append(p)
        pool.close()
        while True:
            if all(r.ready() for r in results):
                return [r.get() for r in results]
            time.sleep(1)
    except KeyboardInterrupt:
        # when you want to kill everything, including this program
        # https://www.reddit.com/r/learnpython/comments/7vwyez/how_to_kill_child_processes_when_using/dtw3oh4/
        pid=os.getpid()
        killtree(pid)
        sys.exit("\nKeyboardInterrupt")


def check_database(dbdir):
    """check existence of database files"""
    if dbdir is None:
        if "CHECKVDB" not in os.environ:
            msg = "Error: database dir not specified\nUse -d or set CHECKDB environmental variable"
            sys.exit(msg)
        else:
            dbdir = os.environ["CHECKVDB"]
    dbdir = os.path.abspath(dbdir)
    if not os.path.exists(dbdir):
        msg = f"Error: database dir not found '{dbdir}'"
        sys.exit(msg)
    files = ["checkv_refs.dmnd", "checkv_refs.tsv"]
    for f in files:
        path = os.path.join(dbdir, f)
        if not os.path.exists(path):
            msg = f"Error: database file not found '{path}'"
            sys.exit(msg)
    return dbdir


def read_fasta(path):
    """ read fasta file and yield (header, sequence)"""
    with open(path) as f:
        for record in SeqIO.parse(f, "fasta"):
            name = record.description
            seq = str(record.seq).upper()
            yield name, seq


def run_prodigal(out):
    cmd = "prodigal "
    cmd += "-p meta "
    cmd += f"-i {out}.fna "
    cmd += f"-a {out}.faa "
    cmd += "1> /dev/null "
    cmd += f"2> {out}.log"
    p = sp.Popen(cmd, shell=True)
    p.wait()


def run_dustmasker(input, out):
    cmd = "dustmasker "
    cmd += f"-in {input} "
    cmd += f"-out {out} "
    cmd += "-outfmt acclist"
    p = sp.Popen(cmd, shell=True)
    p.wait()


def run_diamond(out, db, faa, threads):
    cmd = "diamond blastp "
    cmd += "--outfmt 6 "
    cmd += "--evalue 1e-5 "
    cmd += "--query-cover 50 "
    cmd += "--subject-cover 50 "
    cmd += f"--query {faa} "
    cmd += f"--db {db}/checkv_refs.dmnd "
    cmd += f"--threads {threads} "
    cmd += "-k 10000 "
    cmd += f"> {out} "
    cmd += "2> /dev/null"
    p = sp.Popen(cmd, shell=True)
    p.wait()


def run_hmmsearch(out, db, faa, evalue=10, threads=1):
    cmd = "hmmsearch "
    cmd += "--noali "
    cmd += "-o /dev/null "
    cmd += f"-E {evalue} "
    cmd += f"--tblout {out} "
    cmd += f"--cpu {threads} "
    cmd += f"{db}/checkv_hmms.hmm "
    cmd += f"{faa} "
    p = sp.Popen(cmd, shell=True)
    p.wait()


def search_hmms(out_dir, threads, db_dir):
    # make tmp
    tmp = f"{out_dir}/tmp/hmmsearch"
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    # list faa files
    faa = [
        f
        for f in os.listdir(out_dir + "/tmp/proteins")
        if f.split(".")[-1] == "faa"
    ]
    # run hmmer
    args_list = []
    for f in faa:
        out = f"{tmp}/{f.split('.')[0]}.hmmout"
        args_list.append([out, db_dir, out_dir + "/tmp/proteins/" + f])
    parallel(run_hmmsearch, args_list, threads)
    # cat output
    with open(f"{tmp}.txt", "w") as f:
        for f in os.listdir(tmp):
            with open(f"{tmp}/{f}") as subf:
                for line in subf:
                    f.write(line)


def call_genes(in_fna, out_dir, threads):
    # make tmp
    tmp = f"{out_dir}/tmp/proteins"
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    # count seqs
    num_seqs = sum(1 for _ in read_fasta(in_fna))
    # split fna
    split_size = int(math.ceil(1.0 * num_seqs / threads))
    iteration = 1
    count = 0
    out = open(f"{tmp}/{iteration}.fna", "w")
    for id, seq in read_fasta(in_fna):
        out.write(">" + id + "\n" + seq + "\n")
        count += 1
        if count == split_size:
            count = 0
            iteration += 1
            # avoid creating empty fasta file
            if iteration <= threads:
                out = open(f"{tmp}/{iteration}.fna", "w")
    out.close()
    # call genes
    args_list = []
    for i in range(1, threads + 1):
        out = f"{tmp}/{i}"
        args_list.append([out])
    parallel(run_prodigal, args_list, threads)
    # cat output
    with open(f"{tmp}.faa", "w") as f:
        for i in range(1, iteration + 1):
            # avoid trying to read empty fasta file
            if i <= threads:
                with open(f"{tmp}/{i}.faa") as subf:
                    for line in subf:
                        f.write(line)


def parse_blastp(path):
    with open(path) as f:
        names = [
            "qname",
            "tname",
            "pid",
            "aln",
            "mis",
            "gap",
            "qstart",
            "qstop",
            "tstart",
            "tstop",
            "eval",
            "score",
        ]
        formats = [str, str, float, int, int, int, int, int, int, int, float, float]
        for line in f:
            values = line.split()
            yield dict([(names[i], formats[i](values[i])) for i in range(12)])


def parse_hmmsearch(path):
    with open(path) as f:
        names = [
            "qname",
            "qacc",
            "tname",
            "tacc",
            "eval",
            "score",
            "bias",
            "beval",
            "bscore",
            "bbias",
        ]
        formats = [str, str, str, str, float, float, float, float, float, float]
        for line in f:
            if not line.startswith("#"):
                values = line.split()
                yield dict([(names[i], formats[i](values[i])) for i in range(10)])
