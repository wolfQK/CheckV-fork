# CheckV
Assessing the quality of viral genomes from metagenomes & viromes

## Installation

There are two methods to install CheckV in your computer:

- Using `conda`:  

```bash
conda install -c conda-forge -c bioconda checkv
```

- Using `pip`:

```bash
pip install checkv
```

If you decide to install CheckV via `pip`, make sure you also have the following external dependencies installed:

- BLAST
- DIAMOND
- HMMER
- Prodigal

### CheckV database

Whichever method you choose to install CheckV you will need to download and extract database in order to use it:

```bash
wget https://www.dropbox.com/s/ysv382w01ee7y3t/checkv-db-v0.2.tar.gz
tar -zxvf checkv-db-v0.1.tar.gz
```

Update your environment:

```bash
export CHECKVDB=/path/to/checkv-db-v0.1
```

If you don't want to put the database into your PATH, you can still use it through the `-d` parameter.

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




