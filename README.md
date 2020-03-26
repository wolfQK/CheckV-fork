# CheckV
Assessing the quality of viral genomes from metagenomes & viromes

## Installation

Clone repository:
`git clone bitbucket.org/snayfachlbl/checkv`

Install dependencies:
`conda install -c conda-forge -c bioconda numpy biopython blast diamond hmmer prodigal`

Download database:
`wget https://www.dropbox.com/s/ab1hj499qeypn6z/checkv-db-v0.1.tar.gz`
`tar -zxvf checkv-db-v0.1.tar.gz`

Setup environment:
`export PATH=$PATH:/path/to/CheckV`
`export PYTHONPATH=$PYTHONPATH:/path/to/CheckV`
`export CHECKVDB=/path/to/checkv-db-v0.1`

## Quick start

Navigate to CheckV test directory:
`cd /path/to/CheckV/test`

Identify circular genomes:
`run_checkv.py circularity -i test.fna -o checkv_out`

Estimate completeness for genome fragments:
`run_checkv.py completeness -i test.fna -t 16 -o checkv_out`

Estimate contamination for integrated prophages:
`run_checkv.py contamination -i test.fna -t 16 -o checkv_out`

Generate final output:
`run_checkv.py summary checkv_out`




