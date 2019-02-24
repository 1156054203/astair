# asTair

A little tool for analysis of bisulfite-free and base-resolution sequencing data generated with TET Assisted Pic-borane Sequencing (TAPS) or other modified cytosine to thymine conversion methods (CmtoT).
_These scripts transform aligned reads to bed-like file containing cytosine positions and their modification level._

Contents:

1. [Creating a python virtual environment](#pyenv)
2. [Installation of dependencies with _pip_](#pip)
3. [Creating a conda virtual environment and installing requirements](#conda)
4. [Add asTair to the path and begin the analysis](#path)
5. [Analysis of TAPS  or other modified cytosine to thymine conversion methods (CmtoT) data ](#taps)
6. [Recommendations and extras](#extra)

<a name="pyenv"></a>
### Creating a python virtual environment 

Using python virtualenv and creating an environment folder for your project:

```bash
virtualenv -p /path/to/(python2.7 or python3) your_project_name
```
Activate the environment before proceeding with pip installation:

```bash
source your_project_name/bin/activate
```

To deactivate the environment:
```bash
deactivate
```
To use the python environment you can add it to the path or use it directly:

```bash
/path/to/environment/your_project_name/bin/python(2 or 3) command --options
  
export PATH=$PATH:/path/to/environment/your_project_name/bin/python(2 or 3)
```
<a name="pip"></a>
### Installation of dependencies with _pip_
For python 2.7:
```bash
pip install -r requirements.txt
```
  
For python 3 (>=3.5):
```bash
pip3 install -r requirements.txt
```

<a name="conda"></a>
### Creating a conda virtual environment and installing requirements

For python 2.7:
```bash
conda create -n your_project_name python=2.7
```

For python 3 (>=3.5):
```bash
conda create -n your_project_name python=3.5
```

```bash
conda install -n your_project_name --file requirements.txt
```

Finally, add the the directory of the conda virual environment to the path:

```bash
export PATH=$PATH:/path/to/conda/environment/your_project_name/bin/
```

<a name="path"></a>
### Add asTair to the path and begin the analysis

```bash
export PATH=$PATH:/dir/to/asTair
```

<a name="taps"></a>
### Analysis of TAPS  or other modified cytosine to thymine conversion methods (CmtoT) data 
_Let us have an example paired-end sequencing data from step 1 called `lambda.phage_test_sample_R1.fastq.gz` and `lambda.phage_test_sample_R2.fastq.gz` to guide us through the process._

1 - Get your data hot from the sequencing center.  
2 - Choose genomic aligner of preference. Ours was  `bwa mem`, which we wrapped in `astair_aligner_v3.py`, for which we currently require samtools.

```bash
astair_aligner_v3.py -bp /dir/to/bwa -sp /dir/to/samtools -f lambda.phage.fa -1 lambda.phage_test_sample_R1.fastq.gz -2 lambda.phage_test_sample_R2.fastq.gz -d /output/directory/
```
3 - Sort and index your cram/bam files.

4 - Summon `astair_mod_caller_v3.py` for the modifcation calling. The current version does not require splitting by chromosomes, and can give both at least 1x covered positions only (default mode), or 0x covered and above (`--zero coverage` option). Splitting by cytosine context and/or chromosome prior to the modification calling might be desirable for larger files.  

```bash
astair_mod_caller_v3.py -i lambda.phage_test_sample.bam (cram) -f lambda_phage.fa -co CpG -sc -bq 13 -d /output/directory/
```
  
<a name="extra"></a>
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

### Help with asTair

```bash
Usage: astair_aligner_v3.py [OPTIONS]
Options:
  -1, --fq1 TEXT                  First in pair (R1) sequencing reads file in
                                  fastq.gz format  [required]
  -2, --fq2 TEXT                  Second in pair (R2) sequencing reads file in
                                  fastq.gz format  [required]
  -f, --fasta_file TEXT           DNA sequence in fasta format used for
                                  aligning the sequencing reads.  [required]
  -bp, --bwa_path TEXT            The path to BWA.
  -sp, --samtools_path TEXT       The path to Samtools.
  -d, --directory TEXT            Output directory to save files.  [required]
  -m, --method [CtoT|CmtoT]       Specify sequencing method, possible options
                                  are CtoT (unmodified cytosines are converted
                                  to thymines, bisulfite sequencing-like) and
                                  CmtoT (modified cytosines are converted to
                                  thymines, TAPS-like).
  -mq, --minimum_mapping_quality INTEGER
                                  Set the minimum mapping quality for a read
                                  to be output to file (Default >=1).
  -u, --keep_unmapped             Outputs the unmapped reads (Default false).
  -t, --N_threads INTEGER         The number of threads to spawn (the default
                                  value is 10).  [required]
  --help                          Show this message and exit.
```


```bash
Usage: astair_mod_caller_v3.py [OPTIONS]

Options:
  -i, --input_file TEXT           BAM format file containing sequencing reads.
                                  [required]
  -f, --fasta_file TEXT           DNA sequence in fasta format used for
                                  aligning the sequencing reads and mpileup.
                                  [required]
  -z, --zero_coverage             When set to True, outputs positions not
                                  covered in the bam file. Uncovering zero
                                  coverage positions takes longer time than
                                  using the default option.
  -co, --context [all|CpG|CHG|CHH]
                                  Explains which cytosine sequence contexts
                                  are to be expected in the output file.
                                  Default behaviour is all, which includes
                                  CpG, CHG, CHH contexts and their sub-
                                  contexts for downstream filtering and
                                  analysis.
  -uc, --user_defined_context TEXT
                                  At least two-letter contexts other than CG,
                                  CHH and CHG to be evaluated, will return the
                                  genomic coordinates for the first cytosine
                                  in the string.
  -sc, --skip_clip_overlap        Skipping the random removal of overlapping
                                  bases between paired-end reads. Not
                                  recommended for paired-end libraries, unless
                                  the overlaps are removed prior to calling.
  -bq, --minimum_base_quality INTEGER
                                  Set the minimum base quality for a read base
                                  to be used in the pileup (Default 20).
  -chr, --per_chromosome TEXT     When used, it calculates the modification
                                  rates only per the chromosome given.
                                  (Default None
  -d, --directory TEXT            Output directory to save files.  [required]
  --help                          Show this message and exit.
```

```bash
Usage: astair_phred_values_v3.py [OPTIONS]

Options:
  -1, --fq1 TEXT                  First in pair (R1) sequencing reads file in
                                  fastq.gz format  [required]
  -2, --fq2 TEXT                  Second in pair (R2) sequencing reads file in
                                  fastq.gz format  [required]
  -cm, --calculation_mode [means|absolute]
                                  Gives the mode of computation used for the
                                  Phred scores summary, where means runs
                                  faster. (Default is means)
  -d, --directory TEXT            Output directory to save files.  [required]
  -s, --sample_size INTEGER       The number of reads to sample for the
                                  analysis. (Default 10 000 000)
  -p, --plot                      Phred scores will be visualised and output
                                  as a pdf file. Requires installed
                                  matplotlib.
  -q, --minimum_score INTEGER     Minimum Phred score used for visualisation
                                  only. (Default 15)
  -c, --colors LIST               List of color values used for visualistion
                                  of A, C, G, T, they are given as
                                  color1,color2,color3,color4. Accepts valid
                                  matplotlib color names, RGB and RGBA hex
                                  strings and  single letters denoting color
                                  {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}.
                                  (Default
                                  skyblue,mediumaquamarine,khaki,lightcoral)
  --help                          Show this message and exit.
```

```bash
Usage: astair_mbias_v3.py [OPTIONS]

Options:
  -i, --input_file TEXT      BAM format file containing sequencing reads.
                             [required]
  -d, --directory TEXT       Output directory to save files.  [required]
  -l, --read_length INTEGER  The read length is needed to calculate the
                             M-bias.  [required]
  -m, --method [CtoT|CmtoT]  Specify sequencing method, possible options are
                             CtoT (unmodified cytosines are converted to
                             thymines, bisulfite sequencing-like) and CmtoT
                             (modified cytosines are converted to thymines,
                             TAPS-like).
  -p, --plot                 Phred scores will be visualised and output as a
                             pdf file. Requires installed matplotlib.
  -c, --colors LIST          List of color values used for visualistion of
                             CpG, CHG and CHH modification levels per read,
                             which are given as color1,color2,color3. Accepts
                             valid matplotlib color names, RGB and RGBA hex
                             strings and  single letters denoting color {'b',
                             'g', 'r', 'c', 'm', 'y', 'k', 'w'}. (Default
                             'teal','gray','maroon')
  --help                     Show this message and exit.

```

```bash
Usage: astair_taps_modification_simulation_v3.py [OPTIONS]

Options:
  -f, --fasta_file TEXT           DNA sequence in fasta format used for
                                  aligning the sequencing reads and mpileup.
                                  [required]
  -l, --read_length INTEGER       Desired length of pair-end sequencing reads.
                                  [required]
  -i, --input_file TEXT           Sequencing reads as a Bam file or fasta
                                  sequence to generate reads.  [required]
  -si, --simulation_input [bam|fasta|Ns]
                                  Input file format according to the desired
                                  outcome. Bam files can be generated with
                                  other WGS simulators allowing for sequencing
                                  errors and read distributions or can be
                                  real-life sequencing data; fasta can be used
                                  to generate simple simulated data where only
                                  read positions and modification positions
                                  are known; Ns is of use to generate
                                  statistically possible sequences with
                                  changes in N nucleotides at certain
                                  positions.
  -m, --method [CtoT|CmtoT]       Specify sequencing method, possible options
                                  are CtoT (unmodified cytosines are converted
                                  to thymines, bisulfite sequencing-like) and
                                  CmtoT (modified cytosines are converted to
                                  thymines, TAPS-like).
  -ml, --modification_level INTEGER
                                  Desired modification level; can take any
                                  value between 0 and 100.
  -lb, --library [directional|ptat]
                                  Provide the correct library construction
                                  method.
  -mp, --modified_positions TEXT  Provide a tab-delimited list of positions to
                                  be modified. By default the simulator
                                  randomly modifies certain positions. Please
                                  use seed for replication if no list is
                                  given.
  -co, --context [all|CpG|CHG|CHH]
                                  Explains which cytosine sequence contexts
                                  are to be modified in the output file.
                                  Default behaviour is all, which modifies
                                  positions in CpG, CHG, CHH contexts.
  -uc, --user_defined_context TEXT
                                  At least two-letter contexts other than CG,
                                  CHH and CHG to be evaluated, will return the
                                  genomic coordinates for the first cytosine
                                  in the string.
  -cv, --coverage TEXT            Desired depth of sequencing coverage.
  -r, --region <TEXT INTEGER INTEGER>...
                                  The one-based genomic coordinates of the
                                  specific region of interest given in the
                                  form chromosome, start position, end
                                  position, e.g. chr1 100 2000.
  -d, --directory TEXT            Output directory to save files.  [required]
  -s, --seed INTEGER              An integer number to be used as a seed for
                                  the random generators to ensure replication.
  --help                          Show this message and exit.
```
