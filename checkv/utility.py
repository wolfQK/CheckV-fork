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
import resource
import platform

def max_mem_usage():
    """ Return max mem usage (Gb) of self and child processes """
    max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if platform.system() == 'Linux':
        return (max_mem_self + max_mem_child)/float(1e6)
    else:
        return (max_mem_self + max_mem_child)/float(1e9)


def get_logger(quiet):
    if not quiet:
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


def terminate_tree(pid, including_parent=True):
    parent = psutil.Process(pid)
    for child in parent.children(recursive=True):
        child.terminate()
    if including_parent:
        parent.terminate()


def parallel(function, argument_list, threads):
    """ Based on: https://gist.github.com/admackin/003dd646e5fadee8b8d6 """
    #threads = len(argument_list) ## why is this being defined again here?
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
        terminate_tree(pid)


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
    with open(f"{out}.cmd", "w") as file:
        file.write(cmd+"\n")
    p = sp.Popen(cmd, shell=True)
    p.wait()


def run_dustmasker(input, out):
    cmd = "dustmasker "
    cmd += "-outfmt acclist "
    cmd += f"-in {input} "
    cmd += f"-out {out} "
    cmd += f"2> {out}.log"
    with open(f"{out}.cmd", "w") as file:
        file.write(cmd+"\n")
    p = sp.Popen(cmd, shell=True)
    p.wait()


def run_diamond(out, db, faa, threads):
    cmd = "diamond blastp "
    cmd += "--outfmt 6 "
    cmd += "--evalue 1e-5 "
    cmd += "--query-cover 50 "
    cmd += "--subject-cover 50 "
    cmd += "-k 10000 "
    cmd += f"--query {faa} "
    cmd += f"--db {db}/checkv_refs.dmnd "
    cmd += f"--threads {threads} "
    cmd += f"> {out} "
    cmd += f"2> {out}.log"
    with open(f"{out}.cmd", "w") as file:
        file.write(cmd+"\n")
    p = sp.Popen(cmd, shell=True)
    p.wait()


def run_hmmsearch(out, db, faa, threads=1, evalue=10):
    cmd = "hmmsearch "
    cmd += "--noali "
    cmd += "-o /dev/null "
    cmd += f"-E {evalue} "
    cmd += f"--tblout {out} "
    cmd += f"--cpu {threads} "
    cmd += f"{db} "
    cmd += f"{faa} "
    cmd += f"2> {out}.log "
    with open(f"{out}.cmd", "w") as file:
        file.write(cmd+"\n")
    p = sp.Popen(cmd, shell=True)
    p.wait()

def search_hmms(tmp_dir, threads, db_dir):
    # make tmp
    hmm_dir = os.path.join(tmp_dir, "hmmsearch")
    if not os.path.exists(hmm_dir):
        os.makedirs(hmm_dir)
    # list faa files
    faa = [
        file
        for file in os.listdir(os.path.join(tmp_dir, "proteins"))
        if file.split(".")[-1] == "faa"
    ]
    # list splits to process
    splits = []
    for file in os.listdir(db_dir):
        if file.split(".")[0] != "hmmout": continue
        split = file.split(".")[0]
        out = os.path.join(hmm_dir, f"{split}.hmmout")
        # file doesn't exist; add to list for processing
        if not os.path.exists(out):
            splits.append(split)
        # check if file is complete
        else:
            x = False
            with open(out) as subf:
                for line in subf:
                    if line == "# [ok]\n":
                        x = True
            if not x:
                splits.append(split)
    # run hmmer
    args_list = []
    for split in splits:
        out = os.path.join(hmm_dir, f"{split}.hmmout")
        hmmdb = os.path.join(db_dir, f"{split}.hmm")
        faa = os.path.join(tmp_dir, "proteins.faa")
        args_list.append([out, hmmdb, faa])
    parallel(run_hmmsearch, args_list, threads)
    # check outputs are complete
    complete = []
    for file in os.listdir(hmm_dir):
        if file.split(".")[-1] != "log":
            x = False
            with open(os.path.join(hmm_dir, file)) as subf:
                for line in subf:
                    if line == "# [ok]\n":
                        x = True
            complete.append(x)
    num_fails = complete.count(False)
    if num_fails > 0:
        sys.exit(f"\nError: {num_fails}/80 hmmsearch tasks failed. Program should be rerun")
    # cat output
    with open(os.path.join(tmp_dir, "hmmsearch.txt"), "w") as f:
        for file in os.listdir(hmm_dir):
            if file.split(".")[-1] != "log":
                with open(os.path.join(hmm_dir, file)) as subf:
                    for line in subf:
                        f.write(line)


def call_genes(in_fna, out_dir, threads):
    # make tmp dir
    tmp = f"{out_dir}/tmp/proteins"
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    # count seqs in fasta
    num_seqs = sum(1 for _ in read_fasta(in_fna))
    # split fna into equal sized chunks
    split_size = int(math.ceil(1.0 * num_seqs / threads))
    iteration = 1
    count = 0
    out = open(os.path.join(tmp, f"{iteration}.fna"), "w")
    for id, seq in read_fasta(in_fna):
        # check if new file should be opened
        if count == split_size:
            count = 0
            iteration += 1
            out = open(os.path.join(tmp, f"{iteration}.fna"), "w")
        # write seq to file
        out.write(">" + id + "\n" + seq + "\n")
        count += 1
    out.close()
    # call genes
    args_list = []
    for i in range(1, threads + 1):
        out = os.path.join(tmp, str(i))
        args_list.append([out])
    parallel(run_prodigal, args_list, threads)
    # cat output
    with open(f"{tmp}.faa", "w") as f:
        for i in range(1, iteration + 1):
            # avoid trying to read empty fasta file
            if i <= threads:
                with open(os.path.join(tmp, f"{i}.faa")) as subf:
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

