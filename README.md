![](https://bitbucket.org/berkeleylab/checkv/raw/6d4448f738ac8549551c8ef9511afb05bc394813/logo.png)

[![PyPI](https://img.shields.io/pypi/v/checkv.svg?label=PyPI&color=green)](https://pypi.python.org/pypi/checkv)
[![Conda](https://img.shields.io/conda/vn/bioconda/checkv.svg?label=Conda&color=green)](https://anaconda.org/bioconda/checkv)
[![PyPI downloads](https://img.shields.io/pypi/dm/checkv?label=PyPI%20downloads&color=blue)](https://pypi.python.org/pypi/checkv)
[![Conda downloads](https://img.shields.io/conda/dn/bioconda/checkv.svg?label=Conda%20downloads&color=blue)](https://anaconda.org/bioconda/checkv)
[![Docker pulls](https://img.shields.io/docker/pulls/snayfach/checkv-0.7.0.svg?label=Docker%20pulls&color=blue)](https://hub.docker.com/r/snayfach/checkv-0.7.0)

CheckV is a fully automated command-line pipeline for assessing the quality of single-contig viral genomes, including identification of host contamination for integrated proviruses, estimating completeness for genome fragments, and identification of closed genomes.

The pipeline can be broken down into 4 main steps:

![](https://bitbucket.org/berkeleylab/checkv/raw/657fde9b1c696185a399456fbcbb4ca82066abb6/pipeline.png)

**A: Remove host contamination on proviruses**  

* Genes are first annotated as viral or microbial based on comparison to a custom database of HMMs  
* CheckV scans over the contig (5' to 3') comparing gene annotations and GC content between a pair of adjacent gene windows
* This information is used to compute a score at each intergenic position and identify host-virus breakpoints
* Host-virus breakpoints are identified with:
	* High scores (>1.2)
	* A minimum of 2 host-specific genes in the putative host region (for contigs with >=10 genes)
	* A minimum of 2 virus-specific genes in the putative viral region (for contigs with >=10 genes)
	* A minimum of 30% genes annotated as microbial in the putative host region

**B: Estimate genome completeness (2 algorithms)** 

* AAI-based approach (accurate point estimate for genome completeness; pictured above)
	* First, proteins are compared to the CheckV genome database using AAI (average amino acid identity)
	* After identifying the top hits, completeness is computed as a ratio between the contig length (or viral region length for proviruses) and the length of matched reference genomes
	* A confidence level is reported based on the strength of the alignment and the length of the contig
	* Generally, high- and medium-confidence estimates are quite accurate and can be trusted
* HMM-based approach (estimated range for genome completeness)
	* Highly novel viruses may not match a CheckV genome with sufficient AAI (i.e. low-confidence estimate)
	* In these cases CheckV identifies the viral HMMs on the contig (see panel A) and compares the contig length with reference genomes sharing the same HMMs
	* CheckV then returns the estimated range for genome completeness (e.g. 35% to 60% completeness), which represents the 90% confidence interval based on the distribution of lengths of reference genomes with the same viral HMMs

**C: Predict closed genomes (3 signatures)**

* Direct terminal repeats (DTRs) 
	* Repeated sequence of >20-bp at start/end of contig
	* Most trusted signature in our experience
	* May indicate circular genome or linear genome replicated from a circular template (i.e. concatamer)
* Proviruses 
	* Viral region with predicted host boundaries at 5' and 3' ends (see panel A)
	* Note: CheckV will not detect proviruses if host regions have already been removed (e.g. using VIBRANT or VirSorter)
* Inverted terminal repeats (ITRs)
 	* Repeated sequence of >20-bp at start/end of contig (3' repeat is inverted)
	* Least trusted signature in our experience 

CheckV will also report a confidence level based on comparison to completeness estimates (panel B):  

* High-confidence: >90% estimated completeness
* Medium-confidence: 80-90% estimated completeness
* Low-confidence: <80% estimated completeness

For DTRs and ITRs, CheckV performs some additional filtering/checks:

* Ambiguous bases in repeat (e.g. "NNNNN"): <= 20% of repeat sequence with Ns
* Mode base frequency in repeat (e.g. "AAAAA"): <= 75% of repeat sequence composed of single base 
* Maximum occurences of repeat sequence: <= 8 times per contig (removes highly repetetive sequences)
* Maximum kmer-frquency: <= 1.5 (removes contigs with the same genome repeated back-to-back)


**D: Summarize quality.** 

Based on the results of A-C, CheckV generates a report file and assigns query contigs to one of five quality tiers (consistent with and expand upon the MIUViG quality tiers):
 
* Complete (high- or medium-confidence predictions; see panel C)
* High-quality (>90% completeness)
* Medium-quality (50-90% completeness)
* Low-quality (<50% completeness)
* Undetermined quality


## Installation

There are three methods to install CheckV:

- Using `conda`:

```bash
conda install -c conda-forge -c bioconda checkv
```

- Using `docker`: see section below  

- Using `pip`:

```bash
pip install checkv
```

If you decide to install CheckV via `pip`, make sure you also have the following external dependencies installed:

- DIAMOND (v2.0.4)
- HMMER (v3.3)
- Prodigal (v2.6.3)

The versions listed above were the ones that were properly tested. Different versions may also work.

### CheckV database

If you install using `conda` or `pip` you will need to download the database:

```bash
checkv download_database ./
```

Or (to download and extract manually):

```bash
wget https://portal.nersc.gov/CheckV/checkv-db-v1.0.tar.gz
tar -zxvf checkv-db-v1.0.tar.gz
```

And update your environment (optional):

```bash
export CHECKVDB=/path/to/checkv-db-v1.0
```

Some users may wish to update the database using their own complete genomes:
```bash
checkv update_database /path/to/checkv-db /path/to/updated-checkv-db genomes.fna
```

## Quick start

There are two ways to run CheckV:

- Using individual commands for each step in the pipeline in the following order:

```bash
checkv contamination input_file.fna output_directory -t 16
checkv completeness input_file.fna output_directory -t 16
checkv complete_genomes input_file.fna output_directory
checkv quality_summary input_file.fna output_directory
```

- Using a single command to run the full pipeline:

```bash
checkv end_to_end input_file.fna output_directory -t 16
```

- For a full listing of checkv programs and options, use: `checkv -h` and `checkv <program> -h`

## Docker usage

- Pull the image (including the database)

```bash
docker pull snayfach/checkv-0.7.0
```

- Run the image

```bash
docker run -ti --rm -v "$(pwd):/app" snayfach/checkv-0.7.0 end_to_end input_file.fna output_directory -t 16
```


## Automated testing

CheckV uses the [pytest](https://docs.pytest.org/en/latest/) framework for checking the accuracy of the outputs of its main functions. The testing pipeline executes the `contamination`, `completeness`, `complete_genomes` and `quality_summary` modules sequentially and then compares the output files with the ground truth.

To avoid introducing unwanted breaking changes, we recommend developers and users to execute the tests before committing any changes or submitting new PRs to the repository.

```bash
pytest --database=/path/to/checkv-db --threads=16
```

## Known Issues

* DIAMOND v0.9.30 does not work on OSX
* Excel formatting issue that corrupts data when imported to excel
* (SOLVED) Bug in completeness module (v0.5.1) that suppresses output of HMM-based results


## Output files

#### quality_summary.tsv

This contains integrated results from the three main modules and should be the main output referred to. Below is an example to demonstrate the type of results you can expect in your data:

| contig_id | contig\_length | 	provirus | 	proviral\_length | 	gene_count | 	viral_genes | 	host_genes | 	checkv_quality | 	miuvig_quality | 	completeness | 	completeness\_method | 	complete\_genome_type | 	contamination | 	kmer_freq | 	warnings |     
|---------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|  
| 1 | 	5325 | 	No | 	NA | 	11 | 	0 | 	2 | 	Not-determined | 	Genome-fragment | 	NA | 	NA | 	NA | 	0 | 	1 | 	no viral genes detected | 
| 2 | 	41803 | 	No | 	NA | 	72 | 	27 | 	1 | 	Low-quality | 	Genome-fragment | 	21.99 | 	AAI-based (medium-confidence) | 	NA | 	0 | 	1 | flagged DTR	 | 
| 3 | 	38254 | 	Yes | 	36072 | 	54 | 	23 | 	2 | 	Medium-quality | 	Genome-fragment | 	80.3 | 	HMM-based (lower-bound) | 	NA | 	5.7 | 	1 | 	 | 
| 4 | 	67622 | 	No | 	NA | 	143 | 	25 | 	0 | 	High-quality | 	High-quality | 	100 | 	AAI-based (high-confidence) | 	NA | 	0 | 	1.76 | 	high kmer_freq | 
| 5 | 	98051 | 	No | 	NA | 	158 | 	27 | 	1 | 	Complete | 	High-quality | 	100 | 	AAI-based (high-confidence) | 	DTR | 	0 | 	1 | 	 | 

In the example, above there are results for 6 viral contigs:

* The first 5325 bp contig has no completeness prediction, which is indicated by 'Not-determined' for the 'checkv_quality' field. This contig also has no viral genes identified, so there's a chance it may not even be a virus.   
* The second 41803 bp contig is classified as 'Low-quality' since its completeness is <50%. This is estimate is based on the 'AAI' method. Note that only either high- or medium-confidence estimates are reported in the quality_summary.tsv file. You can see 'completeness.tsv' for more details. This contig had a DTR, but it was flagged for some reason (see complete_genomes.tsv for details)
* The third contig is considered 'Medium-quality' since its completeness is estimated to be 80%, which is based on the 'HMM' method. This means that it was too novel to estimate completeness based on AAI, but shared an HMM with CheckV reference genomes. Note that this value represents a lower bound (meaning the true completeness may be higher but not lower than this value). Note that this contig is also classified as a provirus.
* The fourth contig is classified as High-quality based on a completness of >90%. However, note that value of 'kmer_freq' is 1.7. This indicates that the viral genome is represented multiple times in the contig. These cases are quite rare, but something to watch out for.
* The fifth contig is classified as Complete based on the presence of a direct terminal repeat (DTR) and has 100% completeness based on the AAI method. This sequence can condifently treated as a complete genome.


#### completeness.tsv

A detailed overview of how completeness was estimated:

| contig_id  | contig_length  | proviral_length  | aai\_expected_length  | aai_completeness  | aai_confidence  | aai_error  | aai\_num_hits  | aai\_top_hit  | aai_id  | aai_af  | hmm\_completeness_lower  | hmm\_completeness_upper  | hmm_hits  | 
|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
1  | 9837  | 5713  | 53242.8  | 10.7  | high  | 3.7  | 10  | DTR_517157  | 78.5  | 34.6  | 5  | 15  | 4  | 
2  | 39498  | NA  | 37309  | 100.0  | medium  | 7.7  | 11  | DTR_357456  | 45.18  | 30.46  | 75  | 100 | 22 |
3  | 29224  | NA  | 44960.1  | 65.8  | low  | 15.2  | 17  | DTR_091230  | 39.74  | 19.54  | 52  | 70  | 10  | 
4  | 23404  | NA  | NA  | NA  | NA  | NA  | 0  | NA  | NA  | NA  | NA  | NA  | 0 | 

In the example, above there are results for 4 viral contigs:  

* The first proviral contig has an estimated completeness of 10.7% using on the AAI-based method (100 x 5713 / 53242.8). The confidence for this estimate is high, based on a relative estimated error of 3.7%, which is in turn based on the aai\_id (average amino acid identity) and aai\_af (alignment fraction of contig) to the CheckV reference DTR_517517  
* The second contig has a completeness of 100% using the AAI-based method and a completeness range of 75 - 100% using the HMM-based method. Note that the contig length is a bit longer than the expected genome length of 37,309 bp.
* The third contig is estimated to be 65.8% complete based on the AAI approach. However we can't trust this all that much since the aai_confidence is low (meaning the top hit based on AAI was fairly weak). To be conservative, we may wish to report the range of completeness (52-70%) based on the HMM approach  
* The last contig doesn't have any hits based on AAI and doesn't have any viral HMMs, so there's nothing we can say about this sequence

#### contamination.tsv

A detailed overview of how contamination was estimated:

| contig_id | 	contig_length | 	total_genes | 	viral_genes | 	host_genes | 	provirus | 	proviral_length | 	host_length | 	region_types | 	region_lengths | 	region\_coords_bp | 	region\_coords_genes | 	region\_viral_genes | 	region\_host_genes |  
|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
| 1 | 	98051 | 	158 | 	27 | 	1 | 	No | 	NA | 	NA | 	NA | 	NA | 	NA | 	NA | 	NA | 	NA | 
| 2 | 	38254 | 	54 | 	23 | 	2 | 	Yes | 	36072 | 	2182 | 	host,viral | 	1-2182,2183-38254 | 	1-2182,2183-38254 | 	1-4,5-54 | 	0,23 | 	2,0 | 
| 3 | 	6930 | 	9 | 	1 | 	2 | 	Yes | 	3023 | 	3907 | 	viral,host | 	3023,3907 | 	1-3023,3024-6930 | 	1-5,6-9 | 	1,0 | 	0,2 | 
| 4 | 	101630 | 	103 | 	7 | 	24 | 	Yes | 	28170 | 	73460 | 	host,viral,host | 	46804,28170,26656 | 	1-46804,46805-74974,74975-101630 | 	1-43,44-85,86-103 | 	0,7,0 | 	13,0,11 | 

In the example, above there are results for 4 viral contigs:  

* The first contig is not a predicted provirus
* The second contig has a predicted host region covering 2182 bp
* The third 6930 bp contig has a host region identified on the left side
* The fourth 101630 bp contig has 103 genes, including 7 viral and 24 host genes. CheckV identified two host-virus boundaries

#### complete_genomes.tsv

A detailed overview of putative complete genomes identified:

| contig_id | 	contig_length | 	prediction_type | 	confidence_level | 	confidence_reason | 	repeat_length | 	repeat_count | 
|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
| 1 | 	44824 | 	DTR | 	high | 	AAI-based completeness > 90% | 	253 | 	2 | 
| 2 | 	38147 | 	DTR | 	low | 	Low complexity TR; Repetetive TR | 	20 | 	10 | 
| 3 | 	67622 | 	DTR | 	low | 	Multiple genome copies detected | 	26857 | 	2 | 
| 4 | 	5477 | 	ITR | 	medium | 	AAI-based completeness > 80% | 	91 | 	2 | 
| 5 | 	101630 | 	Provirus | 	not-determined | 	NA | 	NA | 	NA | 

In the example, above there are results for 5 viral contigs:  

* The first viral contig has a direct terminal repeat of 253 bp. It is classified as high-confidence based on having an estimated completeness > 90%
* The second viral contig has a 20-bp DTR, but is the DTR is low complexity and can't be trusted, resulting in a low confidence level. The DTR also occurs 10x and is considered repetetive.
* The third viral contig has DTR of 26857 bp! This indicates that a very large fraction of the genome is repeated. CheckV classifies these as low confidence, but users may which to manually resolve these duplications
* The fourth viral contig has ITR of 91 bp. This is considered medium-confidence based on having AAI-based completeness > 80%
* The fifth viral contig is flanked by host on both sides (provirus). However CheckV was unable to assess completeness, so the confidence is left as not-determined

## Frequently asked questions


**Q: What is the difference between AAI- and HMM-based completeness?**  
A: AAI-based completeness produces a single estimated completeness value, which was designed to be very accurate and can be trusted when the reported confidence level is medium or high. 
HMM-based completeness gives the 90% confidence interval of completeness (e.g. 30-75%) in cases where AAI-based completeness is not reliable. In this example, we can be 90% sure (in theory) that the completeness is between 30% to 75%.

**Q: What is the meaning of the kmer_freq field?**  
A: This is a measure of how many times the viral genome is represented in the contig. Most times this is 1.0 (or very close to 1.0). In rare cases assembly errors may occur in which the contig sequence represents multiple concatenated copies of the viral genome. In these cases genome_copies will exceed 1.0.

**Q: Why does my DTR contig have <100% estimated completeness?**  
A: If the estimated completeness is close to 100% (e.g. 90-110%) then the query is likely complete. However sometimes incomplete genome fragments may contain a direct terminal repeat (DTR), in which case we should expect their estimated completeness to be <90%, and sometimes much less. In other cases, the contig will truly be circular, but the estimated completeness is incorrect. This may also happen if the query a complete segment of a multipartite genome (common for RNA viruses). By default, CheckV uses the 90% completeness cutoff for verification, but a user may wish to make their own judgement in these ambiguous cases.

**Q: Why is my sequence considered "high-quality" when it has high contamination?**  
A: CheckV determines sequence quality solely based on completeness. Host contamination is easily removed, so is not factored into these quality tiers.

**Q: I performed binning and generated viral MAGs. Can I use CheckV on these?**  
A: CheckV can estimate completeness but not contamination for these. You'll need to concatentate the contigs from each MAG into a single sequence prior to running CheckV.

**Q: Can I use CheckV to predict (pro)viruses from whole (meta)genomes?**  
A: Possibly, though this has not been tested.

**Q: Can I use CheckV to remove false positive viral predictions?**  
A: Probably, though this has not been tested.

**Q: How should I handle putative "closed genomes" with no completeless estimate?**  
A: In some cases, you won't be able to verify the completeness of a sequence with terminal repeats or provirus integration sites. DTRs are a fairly reliable indicator (>90% of the time) and can likely be trusted with no completeness estimate. However, complete proviruses and ITRs are much less reliable indicators, and therefore require >90% estimated completeness.

**Q: Why is my contig classified as "undetermined quality"?**  
A: This happens when the sequence doesn't match any CheckV reference genome with high enough similarity to confidently estimate completeness and doesn't have any viral HMMs. There are a few explanations for this, in order of likely frequency: 1) your contig is very short, and by chance it does not share any genes with a CheckV reference, 2) your contig is from a very novel virus that is distantly related to all genomes in the CheckV database, 3) your contig is not a virus at all and so doesn't match any of the references.

**Q: How should I handle sequences with "undetermined quality"?**  
A: While it is not possible to estimate completeness for these, you may choose to still analyze sequences above a certain length (e.g. >30 kb). 

**Q: Does CheckV taxonomically annotate my sequences?**  
A: CheckV does not perform taxonomic annotation of viral contigs. However, some taxonomic information is available for users to find on their own.
You can implement this yourself by looking at the taxonomy (and also the source environment too) for the top hit to the genome database. 
First, you can find the top hit info from the 'aai_top_hit' field in the 'completeness.tsv' output file. 
Then you can look up the taxonomy of the top hit.
If the top hit starts with 'DTR' look here at the 'lineage' field in 'checkv_circular.tsv' database file. 
You can also look at 'habitat' field here as well. 
For GenBank references starting with 'GCA' look at the 'vog_clade' or 'lineage' field in 'checkv_genbank.tsv' database file.

## Supporting code

**Rapid genome clustering based on pairwise ANI**

First, create a blast+ database:  
`makeblastdb -in <my_seqs.fna> -dbtype nucl -out <my_db>`

Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:   
`blastn -query <my_seqs.fna> -db <my_db> -outfmt '6 std qlen slen' -max_target_seqs 10000 -perc_identity 90 -o <my_blast.tsv> -num_threads 32`

Next, calculate pairwise ANI by combining local alignments between sequence pairs:  
`anicalc.py -i <my_blast.tsv> -o <my_ani.tsv>` 

Finally, perform UCLUST-like clustering:  
`aniclust.py --fna <my_seqs.fna> --ani <my_ani.tsv> --out <my_clusters.tsv> --min_ani 95 --min_qcov 0 --min_tcov 70`
