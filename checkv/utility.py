
import os, sys

def check_database(dbdir):
	""" check existence of database files """
	if dbdir is None:
		if 'CHECKVDB' not in os.environ:
			msg = "Error: database dir not specified\nUse -d or set CHECKDB environmental variable"
			sys.exit(msg)
		else:
			dbdir = os.environ['CHECKVDB']
	dbdir = os.path.abspath(dbdir)
	if not os.path.exists(dbdir):
		msg = "Error: database dir not found '%s'" % dbdir
		sys.exit(msg)
	files = ['checkv_refs.dmnd', 'checkv_refs.tsv']
	for file in files:
		path = os.path.join(dbdir, file)
		if not os.path.exists(path):
			msg = "Error: database file not found '%s'" % path
			sys.exit(msg)
	return dbdir

def read_fasta(p):
	""" read fasta file and yield (header, sequence)"""
	with open(p) as fp:
		last = None
		while True:
			if not last:
				for l in fp:
					if l[0] == '>':
						last = l[:-1]
						break
			if not last: break
			name, seqs, last = last[1:], [], None
			for l in fp:
				if l[0] == '>':
					last = l[:-1]
					break
				seqs.append(l[:-1])
			if not last or last[0] != '+':
				yield name, ''.join(seqs)
				if not last: break

def parse_blastp(p):
	names = ['qname', 'tname', 'pid', 'aln', 'mis', 'gap',
	         'qstart', 'qstop', 'tstart', 'tstop', 'eval', 'score']
	formats = [str, str, float, int, int, int, int, int, int, int, float, float]
	with open(p) as f:
		for l in f:
			values = l.split()
			yield dict([(names[i], formats[i](values[i])) for i in range(12)])

def parse_hmmsearch(p):
	names = ['qname', 'qacc', 'tname', 'tacc', 'eval', 'score', 'bias', 'beval',
	         'bscore', 'bbias']
	formats = [str, str, str, str, float, float, float, float, float, float]
	with open(p) as f:
		for l in f:
			if not l.startswith( '#' ):
				values = l.split()
				yield dict([(names[i], formats[i](values[i])) for i in range(10)])
