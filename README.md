![](https://bitbucket.org/berkeleylab/checkv/raw/758a99a857ee874f273c7de326679dfdf7e38847/logo.png)

Assessing the quality of metagenome-assembled viral genomes

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

The versions listed above were the ones that were properly tested. Different versions may also work.

### CheckV database

Whichever method you choose to install CheckV you will need to download and extract database in order to use it:

```bash
wget https://www.dropbox.com/s/xz8h7d1ycrf4fjf/checkv-db-v0.4.tar.gz
tar -zxvf checkv-db-v0.4.tar.gz
```

Update your environment:

```bash
export CHECKVDB=/path/to/checkv-db-v0.4
```

If you don't want to set the environmet variable, you can still use the database through the `-d` parameter of the `contamination` and `completeness` modules.

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




