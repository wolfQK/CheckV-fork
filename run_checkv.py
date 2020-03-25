#!/usr/bin/env python

__version__ = '0.0.1'

import os, sys, subprocess as sp, numpy as np, time, shutil

def get_program():
	""" Get program specified by user """
	if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
		print('CheckV v%s: estimating the quality of viral genomes from metagenomes' % __version__)
		print('')
		print('Usage: checkv.py <command> [options]')
		print('')
		print('Commands:')
		print('   circularity: identify closed, circular genomes')
		print('   completeness: estimate the % genome completeness')
		print('   contamination: estimate the % host contamination for prophages')
		print('')
		print('Note: use checkv.py <command> -h to view usage for a specific command')
		quit()
	else:
		return sys.argv[1]

if __name__ == "__main__":

	# check executables
	fails = 0
	for program in ['dustmasker', 'blastn', 'prodigal', 'diamond', 'hmmsearch']:
		found = False
		for path in os.environ["PATH"].split(os.pathsep):
			exe_file = os.path.join(path.strip('"'), program)
			if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
				found = True
				break
		if not found:
			msg = "Error: required program '%s' not executable or not found on $PATH\n"
			sys.stderr.write(msg % program)
			fails += 1
	if fails > 0: sys.exit()

	program = get_program()
	if program == 'circularity':
		from checkv import dtr
		dtr.main()
	elif program == 'completeness':
		from checkv import completeness
		completeness.main()
	elif program == 'contamination':
		from checkv import contamination
		contamination.main()
	elif program == 'marker':
		from checkv import marker_genes
		marker_genes.main()
	elif program not in ['circular', 'completeness', 'contamination', 'marker']:
		sys.exit("\nError: Unrecognized command: '%s'\n" % program)
