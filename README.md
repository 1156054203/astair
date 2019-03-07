
_`asTair` is a toolchain to process DNA modification sequencing data. `asTair` was designed primarily for handling [TET-Assisted Pyridine Borane (TAPS) sequencing](https://www.nature.com/articles/s41587-019-0041-2) output, but also contains functions that are useful for Bisulfite Sequencing (BS) data._

Contents:

[TOC]

## Installation

### asTair installation through pip

```bash
pip install astair
```


### asTair installation through direct download of the repository

```bash
wget https://bitbucket.org/bsblabludwig/astair/get/v3.0.0.tar.gz
# or if you don't have wget, try
# curl -O https://bitbucket.org/bsblabludwig/astair/get/v3.0.0.tar.gz
tar -xzf v3.0.0.tar.gz
```

#### Creating a python virtual environment 

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

#### Installation with python install

For python 2.7:
```bash
python setup.py install --user
```
For python 3 (>=3.5):
```bash
python3 setup.py install --user
```
#### Using asTair without installation

##### Installation of dependencies with _pip_
For python 2.7:
```bash
pip install -r requirements.txt
```
  
For python 3 (>=3.5):
```bash
pip3 install -r requirements.txt
```

##### Creating a conda virtual environment and installing requirements

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

##### Add asTair to the path and begin the analysis

```bash
export PATH=$PATH:/dir/to/asTair
```

### Downloading the test data

All files needed to starting testing asTair can be downloaded from Zenodo through our bash wrapper `download_test_data.sh`. 

Another classic option is:

```bash
wget -nc -np -nd -A bam,bam.bai,fa,fa.fai,fq.gz  -P path/to/test_data -r https://zenodo.org/record/2582855/
```

## Analysis of TAPS  or other modified cytosine to thymine conversion methods (mCtoT) data 

_Let us have an example paired-end sequencing data from step 1 called `lambda.phage_test_sample_R1.fastq.gz` and `lambda.phage_test_sample_R2.fastq.gz` to guide us through the process._

1 - Get your data hot from the sequencing center.  
2 - Choose genomic aligner of preference. Ours was  `__bwa mem__`, which we wrapped in `aligner.py`, for which we currently require `__samtools__`.

```bash
aligner.py -bp /dir/to/bwa -sp /dir/to/samtools -f lambda.phage.fa -1 lambda.phage_test_sample_R1.fastq.gz -2 lambda.phage_test_sample_R2.fastq.gz -d /output/directory/
```
3 - Sort and index your cram/bam files (done on step two if the aligner.py is used).

4 - Summon `caller.py` `--method mCtoT` (default) for the modifcation calling. The current version does not require splitting by chromosomes, and can give both at least 1x covered positions only (default mode), or 0x covered and above (`--zero coverage` option). Splitting by cytosine context and/or chromosome prior to the modification calling might be desirable for larger files.  

```bash
caller.py -i lambda.phage_test_sample.bam (cram) -f lambda_phage.fa -co CpG -sc False -bq 13 -d /output/directory/
```
  
## Analysis of WGBS  or other unmodified cytosine to thymine conversion methods (CtoT) data 

_The analysis pipeline for bisulfie sequencing data does follows the same steps as TAPS data analysis, but requires different options._

1 - Get your data hot from the sequencing center.  
2 - Choose genomic aligner of preference. Ours was  `__bwa-meth.py__`, which we wrapped in `aligner.py`, for which we currently require `__samtools__` and setting `--method CtoT`.

```bash
aligner.py -bp /dir/to/bwa-meth.py -sp /dir/to/samtools -f lambda.phage.fa -1 lambda.phage_test_sample_R1.fastq.gz -2 lambda.phage_test_sample_R2.fastq.gz --method CtoT -d /output/directory/
```
3 - Sort and index your cram/bam files (done on step two if the aligner.py is used).

4 - Summon `caller.py` `--method CtoT` for the modifcation calling. The current version does not require splitting by chromosomes, and can give both at least 1x covered positions only (default mode), or 0x covered and above (`--zero coverage` option). Splitting by cytosine context and/or chromosome prior to the modification calling might be desirable for larger files.  

```bash
caller.py -i lambda.phage_test_sample.bam (cram) -f lambda_phage.fa --method CtoT -co CpG -sc False -bq 13 -d /output/directory/
```
  


## Recommendations and extras

1 - Do QC check of the sequencing reads and do quality trimming before mapping and dispose of very short reads. 

2 - If the DNA fragmentation during the library construction was random, it will be best to remove PCR duplicates after step 3 before calling samtools mpileup.

3 - Check the fragment (insert) size distribution and decide on overlap removal method for paired-end reads. The simplest option is the default removal of overlaps handled by caller.py, which can be disabled by the `-sc` option.  

```bash
clipOverlap --in lambda.phage_test_sample.bam --out lambda.phage_test_sample_clipped.bam
```

4 - Some more scripts can give your fuller information about the sequencing calling. Currently, in this category are `phred.py` that outputs the average quality per read (if available) for each of the four bases T, C, A, G, and `mbias.py` that gives modification bias along the reads in tabular format and as an image. Another script that might be helpful is `simulator.py` that enables modified data simulation. For visual output from mbias.py and phred.py, installation of __matplotlib__ i required.  

4.1

```bash
phred.py -1 lambda.phage_test_sample_R1.fastq.gz -2 lambda.phage_test_sample_R2.fastq.gz -d /output/directory/ -s 1000 -p
```

4.2

```bash
mbias.py -i lambda.phage_test_sample.bam -l 75 -d /output/directory/ -p
```

4.3

```bash
simulator.py -i lambda.phage_test_sample.bam  -f lambda.phage.fa -l 75 -si bam -ml 40 -co CHG -s 10 -d /output/directory/
```

## asTair general usage

After pip installation:

```bash
Program: asTair (tools for processing cytosine modification sequencing data)
    Version: 3.0

    Usage: Command [command specific options/arguments]

    Commands:
  
    astair_align         Align reads
    
    astair_call          Call methylation
    
    astair_simulate      Simulate TAPS or WGBS data
    
    astair_mbias         Visualise modification bias
    
    astair_phred         Visualise Phred scores
    

    __________________________________About__________________________________
    
    asTair was written by Gergana V. Velikova and Benjamin  Schuster-Boeckler. 
    This code is made available under the GNU General Public License, see 
    LICENSE.txt for more details.
```

After download of the repository scripts need to be added to the path and called directly.

### Help with asTair

```bash
Usage: aligner.py [OPTIONS] | astair_align [OPTIONS] 

Options:
  -1, --fq1 TEXT                  First in pair (R1) sequencing reads file in
                                  fastq.gz format  [required]
  -2, --fq2 TEXT                  Second in pair (R2) sequencing reads file in
                                  fastq.gz format  [required]
  -f, --reference TEXT            Reference DNA sequence in FASTA format used
                                  for aligning of the sequencing reads.
                                  [required]
  -bp, --bwa_path TEXT            The path to BWA for TAPS-like data and to
                                  bwameth.py for bisulfite sequencing.
  -sp, --samtools_path TEXT       The path to Samtools.
  -d, --directory TEXT            Output directory to save files.  [required]
  -m, --method [CtoT|mCtoT]       Specify sequencing method, possible options
                                  are CtoT (unmodified cytosines are converted
                                  to thymines, bisulfite sequencing-like) and
                                  mCtoT (modified cytosines are converted to
                                  thymines, TAPS-like).
  -O, --output_format [BAM|CRAM]  Specify output format, possible options are
                                  BAM and CRAM. The default is CRAM.
  -mq, --minimum_mapping_quality INTEGER
                                  Set the minimum mapping quality for a read
                                  to be output to file (Default >=1).
  -u, --keep_unmapped             Outputs the unmapped reads (Default false).
  -t, --N_threads INTEGER         The number of threads to spawn (the default
                                  value is 1).  [required]
  -k, --minimum_seed_length INTEGER
                                  The minimum seed length used for alignment,
                                  see BWA manual (the default value is 19).
  -w, --band_width INTEGER        The band width for banded alignment, see BWA
                                  manual (the default value is 100).
  -D, --dropoff INTEGER           The off-diagonal X-dropoff, see BWA manual
                                  (the default value is 100).
  -r, --internal_seeds FLOAT      Looks for internal seeds inside a seed
                                  longer than minimum_seed_length *
                                  internal_seeds, see BWA manual (the default
                                  value is 1.5).
  -y, --reseeding_occurence INTEGER
                                  The seed occurrence for the 3rd round
                                  seeding, see BWA manual (the default value
                                  is 20).
  -c, --N_skip_seeds INTEGER      Skips seeds with more than the given seed
                                  occurrences, see BWA manual (the default
                                  value is 500).
  -dc, --drop_chains FLOAT        Drops chains shorter than the specified
                                  fraction of the longest overlapping chain,
                                  see BWA manual (the default value is 0.5).
  -W, --discard_chains INTEGER    Discards a chain if seeded bases shorter
                                  than the specified value, see BWA manual
                                  (the default value is 0).
  -mr, --N_mate_rescues INTEGER   Performs at most the specified rounds of
                                  mate rescues for each read, see BWA manual
                                  (the default value is 50).
  -s, --skip_mate_rescue          NB: Does not recommend unless necessary:
                                  skips mate rescue, in mCtoT mode 
                                  see BWA manual. If set, orphan reads 
                                  (paired reads that are not in a
                                  proper pair) will be generated. Ensure
                                  ignore_orphans in the caller is set to
                                  False.
  -P, --skip_pairing              NB: Does not recommend unless necessary:
                                  skips read pairing in mCtoT mode, 
                                  but does rescue mates unless
                                  mate_skipping is also performed, see
                                  BWA manual. If set, orphan reads (paired
                                  reads that are not in a proper pair) will be
                                  generated. Ensure ignore_orphans in the
                                  caller is set to False.
  -A, --match_score INTEGER       The score for a sequence match, which scales
                                  the remaing scoring options, see BWA
                                  manual(the default value is 1).
  -B, --mismatch_penalty INTEGER  The penalty for a mismatch, see BWA manual
                                  (the default value is 4).
  -o, --gap_open_penalty TEXT     The gap open penalties for deletions and
                                  insertions, see BWA manual (the default
                                  value is 6,6).
  -E, --gap_extension_penalty TEXT
                                  The gap extension penalty with a cost size
                                  calculated as {-O} + {-E}*k, see BWA manual
                                  (the default value is 1,1).
  -L, --end_clipping_penalty TEXT
                                  The penalty for 5-prime- and 3-prime-end
                                  clipping, see BWA manual (the default value
                                  is 5,5).
  -U, --unpaired_penalty INTEGER  The penalty for an unpaired read pair, see
                                  BWA manual (the default value is 17).
  -x, --read_type [null|pacbio|ont2d|intractg]
                                  Changes multiple parameters unless
                                  overridden, see BWA manual (the default
                                  value is None).
  --help                          Show this message and exit.
  
```


```bash
Usage: caller.py [OPTIONS] | astair_call [OPTIONS]

Options:
  -i, --input_file TEXT           BAM|CRAM format file containing sequencing reads.
                                  [required]
  -f, --reference TEXT            Reference DNA sequence in FASTA format used
                                  for aligning of the sequencing reads and for
                                  pileup.  [required]
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
  -m, --method [CtoT|mCtoT]       Specify sequencing method, possible options
                                  are CtoT (unmodified cytosines are converted
                                  to thymines, bisulfite sequencing-like) and
                                  mCtoT (modified cytosines are converted to
                                  thymines, TAPS-like).
  -sc, --skip_clip_overlap BOOLEAN
                                  Skipping the random removal of overlapping
                                  bases between paired-end reads. Not
                                  recommended for paired-end libraries, unless
                                  the overlaps are removed prior to calling.
                                  (Default False)
  -bq, --minimum_base_quality INTEGER
                                  Set the minimum base quality for a read base
                                  to be used in the pileup (Default 20).
  -mq, --minimum_mapping_quality INTEGER
                                  Set the minimum mapping quality for a read
                                  to be used in the pileup (Default 0).
  -amq, --adjust_capq_threshold INTEGER
                                  Used to adjust the mapping quality with
                                  default 0 for no adjustment and a
                                  recommended value for adjustment 50.
                                  (Default 0).
  -mm, --mark_matches BOOLEAN     Output bases matching the reference per
                                  strand (Default True).
  -me, --mark_ends BOOLEAN        Marks head and tail bases in the read
                                  (Default True).
  -ai, --add_indels BOOLEAN       Adds inserted bases and Ns for base skipped
                                  from the reference (Default True).
  -rbq, --redo_baq BOOLEAN        Re-calculates per-Base Alignment Qualities
                                  ignoring existing base qualities (Default
                                  False).
  -cbq, --compute_baq BOOLEAN     Performs re-alignment computing of per-Base
                                  Alignment Qualities (Default True).
  -io, --ignore_orphans BOOLEAN   Ignore reads not in proper pairs (Default
                                  True).
  -md, --max_depth INTEGER        Set the maximum read depth for the pileup,
                                  maximum value 8000 (Default 250).
  -chr, --per_chromosome TEXT     When used, it calculates the modification
                                  rates only per the chromosome given.
                                  (Default None
  -t, --N_threads INTEGER         The number of threads to spawn (the default
                                  value is 1).  [required]
  -d, --directory TEXT            Output directory to save files.  [required]
  --help                          Show this message and exit.

```

```bash
Usage: simulator.py [OPTIONS] | astair_simulate [OPTIONS] 

Options:
  -f, --reference TEXT            Reference DNA sequence in FASTA format used
                                  for generation and modification of the
                                  sequencing reads at desired contexts.
                                  [required]
  -l, --read_length INTEGER       Desired length of pair-end sequencing reads.
                                  [required]
  -i, --input_file TEXT           Sequencing reads as a BAM|CRAM file or fasta
                                  sequence to generate reads.  [required]
  -si, --simulation_input [bam]   Input file format according to the desired
                                  outcome. BAM|CRAM files can be generated with
                                  other WGS simulators allowing for sequencing
                                  errors and read distributions or can be
                                  real-life sequencing data.
  -m, --method [CtoT|mCtoT]       Specify sequencing method, possible options
                                  are CtoT (unmodified cytosines are converted
                                  to thymines, bisulfite sequencing-like) and
                                  mCtoT (modified cytosines are converted to
                                  thymines, TAPS-like). (Default mCtoT)
  -ml, --modification_level INTEGER
                                  Desired modification level; can take any
                                  value between 0 and 100.
  -lb, --library [directional]    Provide the correct library construction
                                  method. NB: Non-directional methods under
                                  development.
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
                                  (Default all)
  -uc, --user_defined_context TEXT
                                  At least two-letter contexts other than CG,
                                  CHH and CHG to be evaluated, will return the
                                  genomic coordinates for the first cytosine
                                  in the string.
  -cv, --coverage INTEGER         Desired depth of sequencing coverage.
  -r, --region <TEXT INTEGER INTEGER>...
                                  The one-based genomic coordinates of the
                                  specific region of interest given in the
                                  form chromosome, start position, end
                                  position, e.g. chr1 100 2000.
  -uc, --user_defined_context TEXT
                                  At least two-letter contexts other than CG,
                                  CHH and CHG to be evaluated, will return the
                                  genomic coordinates for the first cytosine
                                  in the string.
  -ov, --overwrite BOOLEAN        Indicates whether existing output files with
                                  matching names will be overwritten. (Default
                                  False)
  -gc, --GC_bias FLOAT            The value of total GC levels in the read
                                  above which lower coverage will be observed
                                  in Ns and fasta modes. (Default 0.5)
                                  [required]
  -sb, --sequence_bias FLOAT      The proportion of lower-case letters in the
                                  read string for the Ns and fasta modes that
                                  will decrease the chance of the read being
                                  output. (Default 0.1)  [required]
  -t, --N_threads INTEGER         The number of threads to spawn for the bam
                                  mode (the default value is 1).  [required]
  -d, --directory TEXT            Output directory to save files.  [required]
  -s, --seed INTEGER              An integer number to be used as a seed for
                                  the random generators to ensure replication.
  --help                          Show this message and exit.

```

```bash
Usage: mbias.py [OPTIONS] | astair_mbias [OPTIONS]

Options:
  -i, --input_file TEXT      BAM|CRAM format file containing sequencing reads.
                             [required]
  -d, --directory TEXT       Output directory to save files.  [required]
  -l, --read_length INTEGER  The read length is needed to calculate the
                             M-bias.  [required]
  -m, --method [CtoT|mCtoT]  Specify sequencing method, possible options are
                             CtoT (unmodified cytosines are converted to
                             thymines, bisulfite sequencing-like) and mCtoT
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
  -t, --N_threads INTEGER    The number of threads to spawn (the default value
                             is 1).  [required]
  --help                     Show this message and exit.
  
```

```bash
Usage: phred.py [OPTIONS] | astair_phred [OPTIONS]

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

#License

This software is made available under the terms of the [GNU General Public License v3](http://www.gnu.org/licenses/gpl-3.0.html).

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
