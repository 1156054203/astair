
_`asTair` is a toolchain to process DNA modification sequencing data. `asTair` was designed primarily for handling [TET-Assisted Pyridine Borane (TAPS) sequencing](https://www.nature.com/articles/s41587-019-0041-2) output, but also contains functions that are useful for Bisulfite Sequencing (BS) data._

![Build status](https://img.shields.io/bitbucket/pipelines/bsblabludwig/astair.svg "Build Status")

# Basic usage
## 0. Installation 

Installation through `pip` is the easiest way to get `asTair`, and it works in python2 and 3:

```bash
pip install astair
```

You should now be able to call `astair`:

```bash
astair --help
```
```text
Usage: astair [OPTIONS] COMMAND [ARGS]...

  asTair (tools for processing cytosine modification sequencing data)

Options:
  --help  Show this message and exit.

Commands:
  align     Align raw reads in fastq format to a reference genome.
  call      Call modified cytosines from a bam or cram file.
  mbias     Generate modification per read length information (Mbias).
  phred     Calculate per base (A, C, T, G) Phred scores for each strand.
  simulate  Simulate TAPS/BS conversion on top of an existing bam/cram file.

  __________________________________About__________________________________
  asTair was written by Gergana V. Velikova and Benjamin Schuster-Boeckler.
  This code is made available under the GNU General Public License, see
  LICENSE.txt for more details.
                                                           Version: 3.0.0
```

In general, you can use `--help` on all `astair` sub-commands to get detailed instructions on the available options.

(If for some reason `pip` is not an option, [see below](#markdown-header-installing-astair-without-pip) for further ways to install `asTair`.)

## 1. Align reads

We will assume that you have generated paired-end sequencing data, which is stored in two fastq files. For this brief tutorial, we assume the files are called `lambda.phage_test_sample_R1.fq.gz` and `lambda.phage_test_sample_R2.fq.gz`. If you want to follow this tutorial, you can download the files here:

```bash
# Or use curl -O if wget is not available
wget https://zenodo.org/record/2582855/files/lambda.phage_test_sample_1.fq.gz
wget https://zenodo.org/record/2582855/files/lambda.phage_test_sample_2.fq.gz
```

The raw reads need to be aligned. asTair contains a command to help with this. It assumes that [`bwa`](https://github.com/lh3/bwa) and [`samtools`](http://www.htslib.org/) are available on your system. (If you prefer to use a different aligner, [skip to step 2](#markdown-header-2-call-methylation).)

You will also need an indexed reference genome to align to. For this example we are using the lambda phage genome, which you can download with

```bash
wget https://zenodo.org/record/2582855/files/lambda_phage.fa
wget https://zenodo.org/record/2582855/files/lambda_phage.fa.fai
```

Now, you are ready to align:
```bash
mkdir -p output_dir
astair align -f lambda_phage.fa -1 lambda.phage_test_sample_1.fq.gz -2 lambda.phage_test_sample_2.fq.gz -d output_dir
```

## 2. Call methylation

Once your fastq files are aligned and sorted (done automatically by `astair align`), you can run `astair call` to create a list of putative modified positions:

```bash
astair call -i output_dir/lambda.phage_test_sample_mCtoT.cram -f lambda_phage.fa --context CpG --minimum_base_quality 13 -d output_dir/
```

## 3. Interpret results
After calling methylation, you will find two additional files in `output_dir`:

1. `lambda.phage_test_sample_mCtoT_mCtoT_CpG.stats`
2. `lambda.phage_test_sample_mCtoT_mCtoT_CpG.mods`

The `.stats` file contains global statistics on the modification rate in different sequence contexts. You can use this to get an idea of the overall level of modification in your sample. In our example here, we used a 1:1 mixture of in-vitro modified and unmodified lambda phage, so the results show a methylation rate of approx. 50% :


| CONTEXT | SPECIFIC_CONTEXT | MEAN_MODIFICATION_RATE_PERCENT | TOTAL_POSITIONS | COVERED_POSITIONS |
| ------- | ---------------- | ------------------------------ | --------------- | ----------------- |
| CpG     |                  | **48.225**                     | 6225            | 6225              |
|         | CGA              | 44.647                         | 1210            | 1210              |
|         | CGC              | 48.595                         | 1730            | 1730              |
|         | CGG              | 48.936                         | 1847            | 1847              |
|         | CGT              | 49.862                         | 1438            | 1438              |

The `.mods` file contains per-position information on your sample:

| CHROM | START | END   | MOD_LEVEL | MOD    | UNMOD   | REF   | ALT  | CONTEXT | SPECIFIC_CONTEXT  | SNV     | TOTAL_DEPTH |
| ----- | ----- | ----- | --------- | ------ | ------- | ----- | ---- | --------| ----------------- | ------- | ----------- |
| lambda |3 | 4 | 1.0 | 23 | 0 | C | T | CGG | CpG | No | 57 |
| lambda |4 | 5 | 0.0 | 0 | 34 | G | A | CGC | CpG | No | 71 |
| lambda  |6 | 7 | 1.0 | 38 | 0 | C | T | CGA | CpG | No | 104 |
| lambda  |7 | 8 | 1.0 | 58 | 0 | G | A | CGC | CpG | No | 127 |
| lambda | 12 | 13 | 1.0 | 88 | 0 | C | T | CGC | CpG | No | 240 |
| lambda  |13 | 14 | 0.0 | 0 | 139 | G | A | CGA | CpG | No | 250 |

The header should be mostly self-explanatory. `MOD` and `UNMOD` refer to the number of reads covering that base that showed evidence of modification/no modification, and were of the right orientation to be meaningful for modification calling. The total coverage, including reads that were oriented in a way that no modification information can be extracted, is shown in `TOTAL_DEPTH`. `SNV` gives a heuristic indication whether the position is indeed a modified base, or a genetic C to T variant in the genome of the sample.

# Other useful information

## Recommendations for data pre-processing

1. Do quality control of the sequencing reads and do quality trimming before mapping and dispose of very short reads, using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) or similar tools.
2. In most cases, it will be best to remove PCR duplicates before running the modification caller, unless your reads are non-randomly fragmented (e.g. enzymatically digested).
3. Check the fragment (insert) size distribution and decide on an overlap removal method for paired-end reads. The simplest option is the default removal of overlaps handled by `astair call`, which will randomly select one of two overlapping reads. This behaviour can be disabled by the `-sc` option, in case you are using a more sophisticated overlap-clipping tool.

## Installing asTair without `pip`

If pip is not available for some reason, you can simply download the package directly and install it manually:

```bash
wget https://bitbucket.org/bsblabludwig/astair/get/v3.0.2.tar.gz
# or if you don't have wget, try
# curl -O https://bitbucket.org/bsblabludwig/astair/get/v3.0.2.tar.gz

tar -xzf v3.0.2.tar.gz -C astair_3.0.2 --strip-components=1

cd astair_3.0.2
```

We would recommend using a virtual environment to avoid issues with globally installed packages, or if you are on a shared system and need to install locally:

```bash
virtualenv venv

# Activate the environment
source venv/bin/activate
```

To install, execute:

```bash
python setup.py install
# or, if you want to install as a user, but without using a virtual environment:
# python setup.py install --user
```

Now, you can perform all analysis and processing steps described [above](#markdown-header-basic-usage).

If you were using a virtual environment, you can deactivate the environment after you're finished using astair by calling:

```bash
deactivate
```

## Analysis of WGBS data (or other unmodified cytosine to thymine conversion methods)

The analysis pipeline for bisulfie sequencing data does follows the same steps as TAPS data analysis, but requires different options. We again start from fastq files. To avoid Bismark-style double-alignments, we prefer to use `bwa meth`, which can be used directly through `astair align` when you choose the `--method CtoT` option.

```bash
mkdir -p output_dir
astair align -f lambda_phage.fa -1 lambda.phage_test_sample_BS_1.fastq.gz -2 lambda.phage_test_sample_BS_2.fastq.gz --method CtoT -d output_dir/
```

You can now use `astair call` with `--method CtoT` for the modifcation calling:
```bash
astair call -i output_dir/lambda.phage_test_sample_BS_CtoT.cram -f lambda_phage.fa --method CtoT --context CpG --minimum_base_quality 13 -d output_dir/
```

# License

This software is made available under the terms of the [GNU General Public License v3](http://www.gnu.org/licenses/gpl-3.0.html).

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
