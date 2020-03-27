# CheckV
Assessing the quality of viral genomes from metagenomes & viromes

## Installation

Clone repository:
```bash
git clone bitbucket.org/snayfachlbl/checkv
```

Install dependencies:
```bash
conda install -c conda-forge -c bioconda biopython numpy psutil blast diamond hmmer prodigal
```

Download database:
```bash
wget https://www.dropbox.com/s/ab1hj499qeypn6z/checkv-db-v0.1.tar.gz
tar -zxvf checkv-db-v0.1.tar.gz
```

Setup environment:
```bash
export CHECKVDB=/path/to/checkv-db-v0.1
```

## Quick start

Navigate to CheckV test directory:
```bash
cd /path/to/CheckV/test
```

Identify circular genomes:
```bash
checkv circularity test.fna checkv_out
```

Estimate completeness for genome fragments:
```bash
checkv completeness test.fna checkv_out -t 16
```

Estimate contamination for integrated prophages:
```bash
checkv contamination test.fna checkv_out -t 16
```

Generate final output (**not ready**):
```bash
checkv summary checkv_out
```




