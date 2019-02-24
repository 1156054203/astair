# asTair

## A little tool for analysis of bisulfite-free and base-resolution sequencing data generated with TET Assisted Pic-borane Sequencing (TAPS) or other modified cytosine to thymine conversion methods.
_These scripts transform aligned reads to bed-like file containing cytosine positions and their modification level._

### Installation of dependencies with _pip_
For python 2.7:
```bash
pip install -r requirements.txt
```
  
For python 3 (>=3.5):
```bash
pip3 install -r requirements.txt
```

### Add asTair to the path and begin the analysis

export PATH=$PATH:/dir/to/asTair

### Analysis of TAPS data 
_Let us have an example paired-end sequencing data from step 1 called `lambda.phage_test_sample_R1.fastq.gz` and `lambda.phage_test_sample_R2.fastq.gz` to guide us through the process._

1 - Get your data hot from the sequencing center.  
2 - Choose genomic aligner of preference. Ours was  `bwa mem`, which we wrapped in astair_aligner_v3.py, for which we currently require samtools.

```bash
astair_aligner_v3.py -bp /dir/to/bwa -sp /dir/to/samtools -f lambda.phage.fa -1 lambda.phage_test_sample_R1.fastq.gz -2 lambda.phage_test_sample_R2.fastq.gz -d /output/directory/
```
3 - Sort and index your cram/bam files.

4 - Summon astair_mod_caller_v3.py for the modifcation calling. The current version does not require splitting by chromosomes, and can give both at least 1x covered positions only (default mode), or 0x covered and above ('--zero coverage' option). Splitting by cytosine context and/or chromosome prior to the modification calling might be desirable for larger files.  

```bash
astair_mod_caller_v3.py -i lambda.phage_test_sample.bam (cram) -f lambda_phage.fa -co CpG -sc -bq 13 -d /output/directory/
```
  
### Recommendations and extras

1 - Do QC check of the sequencing reads and do quality trimming before mapping and dispose of very short reads. 

2 - If the DNA fragmentation during the library construction was random, it will be best to remove PCR duplicates after step 3 before calling samtools mpileup.

3 - Check the fragment (insert) size distribution and decide on overlap removal method for paired-end reads. The simplest option is the default removal of overlaps handled by astair_mod_caller_v3.py, which can be disabled by the `-sc` option.  

```bash
clipOverlap --in lambda.phage_test_sample.bam --out lambda.phage_test_sample_clipped.bam
```

4 - Some more scripts can give your fuller information about the sequencing calling. Currently, in this category are `astair_phred_values_v3.py` that outputs the average quality per read (if available) for each of the four bases T, C, A, G, and `astair_mbias_v3.py` that gives modification bias along the reads in tabular format and as an image. Another script that might be helpful is `astair_taps_modification_simulation_v3.py` that enables modified data simulation.  

4.1

```bash
astair_phred_values_v3.py -1 lambda.phage_test_sample_R1.fastq.gz -2 lambda.phage_test_sample_R2.fastq.gz -d /output/directory/ -s 1000 -p
```

4.2

```bash
astair_mbias_v3.py -i lambda.phage_test_sample.bam -l 75 -d /output/directory/ -p
```

4.3

```bash
astair_taps_modification_simulation_v3.py -i lambda.phage_test_sample.bam  -f lambda.phage.fa -l 75 -si bam -ml 40 -co CHG -s 10 -d /output/directory/
```
