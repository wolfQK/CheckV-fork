# CheckV
Assessing the quality of viral genomes from metagenomes & viromes

## Installation

Clone repository:
```bash
git clone bitbucket.org/berkeleylab/checkv
```

Install dependencies:
```bash
conda install -c conda-forge -c importlib_metadata bioconda biopython numpy psutil blast diamond hmmer prodigal
```

Install CheckV:
```bash
cd /path/to/checkv  
python setup.py develop
```

Download database:
```bash
wget https://www.dropbox.com/s/s2pymkcdyujp47g/checkv-db-v0.1.tar.gz
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

Identify circular genomes:
```bash
checkv circularity test.fna checkv_out --verbose
```

Estimate completeness for genome fragments:
```bash
checkv completeness test.fna checkv_out -t 16 --verbose
```

Estimate contamination for integrated prophages:

* using full database (slower):
```bash
checkv contamination test.fna checkv_out -t 16 --hmm-db full --verbose
```

* or using reduced database (faster):
```bash
checkv contamination test.fna checkv_out -t 16 --hmm-db reduced --verbose
```

Generate final output (**not ready**):
```bash
checkv summary checkv_out
```




