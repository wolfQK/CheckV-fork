# CheckV
Assessing the quality of viral genomes from metagenomes & viromes

## Installation

Download repository:  
```bash
git clone bitbucket.org/berkeleylab/checkv
```

Install dependencies:  
```bash
conda install -c conda-forge -c importlib_metadata bioconda biopython numpy psutil blast diamond=0.9.30 hmmer prodigal
```

Install CheckV:  
```bash
cd /path/to/checkv
python setup.py develop
```

Download database:  
```bash
wget https://www.dropbox.com/s/ysv382w01ee7y3t/checkv-db-v0.2.tar.gz
tar -zxvf checkv-db-v0.1.tar.gz
```

Setup environment:  
```bash
export CHECKVDB=/path/to/checkv-db-v0.1
```

## Quick start

Navigate to CheckV test directory:  
```bash
cd /path/to/checkv/test
```

Estimate host contamination on integrated prophages:
```bash
checkv contamination test.fna checkv_out -t 16
```

Estimate completeness for genome fragments:
```bash
checkv completeness test.fna checkv_out -t 16
```

Identify (possible) complete genomes with terminal repeats:
```bash
checkv terminal_repeats test.fna checkv_out
```

Summarize CheckV output & classify contigs into quality tiers:
```bash
checkv quality_summary test.fna checkv_out
```




