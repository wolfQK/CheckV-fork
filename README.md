<p align="center"><a href="https://bitbucket.org/berkeleylab/checkv/src/master/"><img src="https://bitbucket.org/berkeleylab/checkv/raw/f739cd5c9622b9b1f799b52115c68d8387f17d25/logo.png" width="310rem"></a></p>

<p align="center">Assessing the quality of metagenome-assembled viral genomes</p>

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

- BLAST+ (v2.5.0)
- DIAMOND (v0.9.30)
- HMMER (v3.3)
- Prodigal (v2.6.3)

### CheckV database

Whichever method you choose to install CheckV you will need to download and extract database in order to use it:

```bash
wget https://www.dropbox.com/s/ysv382w01ee7y3t/checkv-db-v0.2.tar.gz
tar -zxvf checkv-db-v0.1.tar.gz
```

Update your environment:

```bash
export CHECKVDB=/path/to/checkv-db-v0.2
```

If you don't want to set the environmet variable, you can still use the database through the `-d` parameter.

## Quick start

Navigate to CheckV test directory:  

```bash
cd /path/to/checkv/test
```

Identify flanking host regions on integrated prophages:

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




