#!/usr/bin/env python

import csv, os, sys, time, argparse, subprocess as sp, numpy as np
from operator import itemgetter
from bisect import bisect_left
from checkv import utility, prodigal

class Genome:
	def __init__(self):
		pass

class Gene:
	def __init__(self):
		pass

def parse_arguments():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS
	)
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('-i', dest='fna', type=str, required=True, metavar='PATH',
		help="""Input nucleotide sequences in FASTA format""")
	parser.add_argument('-o', dest='out', type=str, required=True, metavar='PATH',
		help="""Output directory""")
	parser.add_argument('-d', dest='db', type=str, required=False, metavar='PATH',
		help="""Reference database path; by default CHECKDB environmental variable used""")
	parser.add_argument('-t', dest='threads', type=int, default=1, metavar='INT',
		help="""Number of threads to use for DIAMOND (1)""")
	parser.add_argument('--restart', action='store_true', default=False,
		help="""Restart program from start overwriting existing files""")
	parser.add_argument('--percent_of_top_hit', type=float, default=50, metavar='FLOAT',
		help=argparse.SUPPRESS)
	parser.add_argument('--exclude_identical', action='store_true', default=False,
		help=argparse.SUPPRESS)
	parser.add_argument('--exclude_list', type=str, default=None, metavar='PATH',
		help=argparse.SUPPRESS)

	args = vars(parser.parse_args())
	return args

def yield_query_alns(p):
	handle = utility.parse_blastp(p)
	alns = [next(handle)]
	query = alns[0]['qname'].rsplit('_', 1)[0]
	for r in handle:
		if r['qname'].rsplit('_', 1)[0] != query:
			yield query, alns
			alns = []
			query = r['qname'].rsplit('_', 1)[0]
		alns.append(r)
	yield query, alns

def take_closest(myList, myNumber):
    """ Assumes myList is sorted. Returns closest value to myNumber.
    	If two numbers are equally close, return the smallest number.
    """
    myNumber = float(myNumber)
    pos = bisect_left(myList, float(myNumber))
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

def compute_aai(blastp_path, out_path, genomes, refs):
	""" Compute AAI to reference genomes
	1. loop over alignment blocks from <blastp_path> (1 per query)
	2. compute average aa identity to each reference genome
	3. write all hits to disk
	"""
	with open(out_path, 'w') as out:
	
		header = ['query', 'target', 'aligned_genes', 'percent_genes', 'aligned_length', 'percent_length', 'identity', 'score']
		out.write('\t'.join(header)+'\n')

		for query, alns in yield_query_alns(blastp_path):
			
			# get best hits for each gene to each target genome
			hits = {}
			for r in alns:
				target = r['tname'].rsplit('_', 1)[0]
				# store gene aln
				if target not in hits:
					hits[target] = {}
				if r['qname'] not in hits[target]:
					hits[target][r['qname']] = r
				elif r['score'] > hits[target][r['qname']]['score']:
					hits[target][r['qname']] = r
			
			# compute aai
			aai = []
			for target, alns in hits.items():
				pids = [r['pid'] for r in alns.values()]
				lens = [r['aln'] - r['gap'] for r in alns.values()]
				aligned_length = sum(lens)
				aligned_genes = len(alns)
				identity = round(sum([x * y for x, y in zip(pids, lens)]) / aligned_length, 2)
				percent_length = round(100.0 * aligned_length / genomes[query].protlen, 2)
				percent_genes = round(100.0 * aligned_genes / genomes[query].genes, 2)
				score = round(identity * aligned_length / 100,2)
				row = [query, target, aligned_genes, percent_genes, aligned_length, percent_length, identity, score]
				aai.append(row)
			if len(aai) == 0: continue
		
			# write all hits to disk
			aai = sorted(aai, key=itemgetter(-1), reverse=True)
			top_score = max([_[-1] for _ in aai])
			for row in aai:
				score = row[-1]
				out.write('\t'.join([str(_) for _ in row])+'\n')

def main():
	
	program_start = time.time()
	args = parse_arguments()
	args['db'] = utility.check_database(args['db'])
	args['tmp'] = os.path.join(args['out'], 'tmp')
	if not os.path.exists(args['out']):
		os.makedirs(args['out'])
	if not os.path.exists(args['tmp']):
		os.makedirs(args['tmp'])

	print("calling genes with prodigal...")
	args['faa'] = os.path.join(args['tmp'], 'proteins.faa')
	if args['restart'] or not os.path.exists(args['faa']):
		prodigal.call_genes(args['fna'], args['out'], args['threads'])

	print("running blastp search...")
	args['blastp'] = os.path.join(args['tmp'], 'diamond.tsv')
	if args['restart'] or not os.path.exists(args['blastp']):
		cmd = "diamond blastp "
		cmd += "--outfmt 6 "
		cmd += "--evalue 1e-5 "
		cmd += "--query-cover 50 "
		cmd += "--subject-cover 50 "
		cmd += "--query %s " % args['faa']
		cmd += "--db %s/checkv_refs.dmnd " % args['db']
		cmd += "--threads %s " % args['threads']
		cmd += "-k 10000 "
		cmd += "> %s " % args['blastp']
		cmd += "2> /dev/null"
		p = sp.Popen(cmd, shell=True)
		p.wait()

	print("initializing queries and database...")
	genomes = {}
	for r in utility.read_fasta(args['fna']):
		header, seq = r
		genome = Genome()
		genome.id = header.split()[0]
		genome.length = len(seq)
		genome.genes = 0
		genome.protlen = 0
		genome.aai = []
		genomes[genome.id] = genome
	# read input genes
	genes = {}
	for r in utility.read_fasta(args['faa']):
		header, seq = r
		gene = Gene()
		gene.id = header.split()[0]
		gene.length = len(seq)
		gene.genome_id = gene.id.rsplit('_', 1)[0]
		genes[gene.id] = gene
		genomes[gene.genome_id].genes += 1
		genomes[gene.genome_id].protlen += gene.length
	# read reference genomes
	refs = {}
	p = os.path.join(args['db'], 'checkv_refs.tsv')
	for r in csv.DictReader(open(p), delimiter='\t'):
		genome = Genome()
		genome.id = r['checkv_id']
		genome.length = int(r['length'])
		refs[genome.id] = genome

	# exclude references
	exclude = set([])
	for l in open(args['db']+'/exclude_genomes.list'):
		exclude.add(l.rstrip())
	if args['exclude_list']:
		for l in open(args['exclude_list']):
			exclude.add(l.rstrip())

	# estimated error rates for alignment cutoffs
	error_rates = {}
	p = os.path.join(args['db'], 'error_rates.tsv')
	for r in csv.DictReader(open(p), delimiter='\t'):
		key = int(r['length']), int(r['aai']), int(r['cov'])
		error_rates[key] = float(r['error']) if int(r['count']) > 100 else 'NA'
	error_keys = {}
	error_keys['length'] = sorted(list(set([_[0] for _ in error_rates.keys()])))
	error_keys['aai'] = sorted(list(set([_[1] for _ in error_rates.keys()])))
	error_keys['cov'] = sorted(list(set([_[2] for _ in error_rates.keys()])))

	print("computing aai...")
	args['aai'] = os.path.join(args['tmp'], 'aai.tsv')
	if args['restart'] or not os.path.exists(args['aai']):
		compute_aai(args['blastp'], args['aai'], genomes, refs)

	print("store aai...")
	for r in csv.DictReader(open(args['aai']), delimiter='\t'):
		
		# better formatting here to floats
		r['identity'] = float(r['identity'])
		r['aligned_genes'] = float(r['aligned_genes'])
		r['aligned_length'] = int(r['aligned_length'])
		r['score'] = float(r['score'])
		r['percent_length'] = float(r['percent_length'])

		if r['target'] in exclude:
			continue
		elif args['exclude_identical'] and r['identity'] == 100 and r['percent_length'] == 100:
			continue
		elif len(genomes[r['query']].aai) == 0:
			genomes[r['query']].aai.append(r)
		else:
			top_score = genomes[r['query']].aai[0]['score']
			if 100.0*(top_score-r['score'])/top_score <= args['percent_of_top_hit']:
				genomes[r['query']].aai.append(r)

	print("estimating completeness...")
	module_start = time.time()
	p = os.path.join(args['out'], 'completeness.tsv')
	with open(p, 'w') as out:
		header = ['genome_id', 'genome_length',
		          'num_hits', 'ref_length_min', 'ref_length_q1', 'ref_length_mean', 'ref_length_q2', 'ref_length_max',
		          'top_hit', 'top_hit_length', 'top_hit_aai', 'top_hit_cov', 'error_rate']
		out.write('\t'.join(header)+'\n')
		
		for genome in genomes.values():
			
			# at least 1 hit to a reference
			if len(genome.aai) > 0:
				
				# fetch top hit
				top = genome.aai[0]
				key1 = take_closest(error_keys['length'], genome.length)
				key2 = take_closest(error_keys['aai'], top['identity'])
				key3 = take_closest(error_keys['cov'], top['percent_length'])
				error_rate = error_rates[key1, key2, key3]
				
				# estimate genome size
				scores = [_['score'] for _ in genome.aai]
				lengths = [refs[_['target']].length for _ in genome.aai]
				avg_len = sum([l*s for l,s in zip(lengths, scores)])/sum(scores)
				len_75, len_25 = np.percentile(lengths, [75 ,25])
				min_len, max_len = min(lengths), max(lengths)
				num_hits = len(lengths)
			
				# write
				row = [genome.id, genome.length]
				row += [num_hits, min_len, len_25, avg_len, len_75, max_len]
				row += [top['target'], refs[top['target']].length, top['identity'], top['percent_length'], error_rate]
				out.write('\t'.join([str(_) for _ in row])+'\n')
				
			# no hits to any reference
			else:
				row = [genome.id, genome.length] + ['NA'] * 11
				out.write('\t'.join([str(_) for _ in row])+'\n')

	# done!
	print("done!")

