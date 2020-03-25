
import os
import shutil
import subprocess as sp
from checkv import utility
from math import ceil

def parallel(function, argument_list, threads):
	""" Based on: https://gist.github.com/admackin/003dd646e5fadee8b8d6 """
	import multiprocessing as mp, signal, time
	def init_worker():
		signal.signal(signal.SIGINT, signal.SIG_IGN)
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
		pool.terminate()
		pool.join()
		sys.exit("\nKeyboardInterrupt")

def run_prodigal(out):
	cmd = "prodigal -i %s.fna -a %s.faa -p meta 1> /dev/null 2> %s.log" % (out, out, out)
	p = sp.Popen(cmd, shell=True)
	p.wait()

def run_hmmer(out, db, faa):
	cmd = "hmmsearch "
	cmd += "--noali "
	cmd += "--tblout %s " % out
	cmd += "--cpu 1 "
	cmd += "%s/checkv_hmms.hmm " % db
	cmd += "%s " % faa
	cmd += "&> /dev/null"
	p = sp.Popen(cmd, shell=True)
	p.wait()

def search_hmms(out_dir, threads, db_dir):
	
	# make tmp
	tmp = '%s/tmp/hmmsearch' % out_dir
	if not os.path.exists(tmp):
		os.makedirs(tmp)
	
	# list faa files
	faa = []
	for file in os.listdir(out_dir+'/tmp/proteins'):
		if file.split('.')[-1] == 'faa':
			faa.append(file)

	# run hmmer
	args_list = []
	for file in faa:
		out = '%s/%s.hmmout' % (tmp, file.split('.')[0])
		args_list.append([out, db_dir, out_dir+'/tmp/proteins/'+file])
	parallel(run_hmmer, args_list, threads)

	# cat output
	with open('%s.txt' % tmp, 'w') as f:
		for file in os.listdir(tmp):
			for l in open('%s/%s' % (tmp, file)):
				f.write(l)


def call_genes(in_fna, out_dir, threads):
	
	# make tmp
	tmp = '%s/tmp/proteins' % out_dir
	if not os.path.exists(tmp):
		os.makedirs(tmp)
	
	# count seqs
	num_seqs = 0
	for seq in utility.read_fasta(in_fna):
		num_seqs += 1
	
	# split fna
	split_size = int(ceil(1.0*num_seqs/threads))
	iter = 1
	count = 0
	out = open('%s/%s.fna' % (tmp, iter), 'w')
	for id, seq in utility.read_fasta(in_fna):
		out.write('>'+id+'\n'+seq+'\n')
		count += 1
		if count == split_size:
			count = 0
			iter += 1
			out = open('%s/%s.fna' % (tmp, iter), 'w')
	out.close()
	
	# call genes
	args_list = []
	for i in range(1,iter+1):
		out = '%s/%s' % (tmp, i)
		args_list.append([out])
	parallel(run_prodigal, args_list, threads)

	# cat output
	with open('%s.faa' % tmp, 'w') as f:
		for i in range(1,iter+1):
			for l in open('%s/%s.faa' % (tmp, i)):
				f.write(l)


