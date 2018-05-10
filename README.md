# asTair

## A little tool for analysis of bisulfite-free and base-resolution sequencing data generated with TET Assisted Pic-borane Sequencing (TAPS).
_These scripts transform aligned reads to bed-like file containing cytosine positions and their modification level._

### Analysis of TAPS data 
_Let us have an example paired-end sequencing data from step 1 called TAPS_mESC_L001_R1.fastq.gz and TAPS_mESC_L001_R2.fastq.gz to guide us through the process._  

1. Get your data hot from the sequencing center.  
2. Choose genomic aligner of preference. Ours was *bwa mem*.  
  ```
  bwa mem /dir/mm9.fa TAPS_mESC_L001_R1.fastq.gz TAPS_mESC_L001_R2.fastq.gz  | samtools view -bS > TAPS_mESC_L001.bam
  ```
3. Use MF_split.py to separate the sequencing reads coming from the original top (OT) and original bottom (OB) strand.  
  ```
 python3 MF_split.py --input_file TAPS_mESC_L001.bam
  ```
4. Call variants on the OT and OB strand separately.  
  ```
  samtools mplieup -t "DP,ADF,ADR,AD" -v -f /dir/mm9.fa  TAPS_mESC_LOO1_OT/OB.bam > TAPS_mESC_LOO1_OT/OB.vcf.gz
  ```
5. Summon MF_callerMOD.py for the modifcation calling. The current version does not require splitting by chromosomes, and can give both at least 1x covered positions only (default mode), or 0x covered and above ('--zero coverage option'). The default mode has some speed advantage, so if running with 0x option, splitting by chromosome prior to the modification calling might be useful.  
  ```
  python3 MF_callerMOD.py --input_file TAPS_mESC_LOO1_OT/OB.vcf.gz --fasta_file /dir/mm9.fa (--zero_coverage)
  ```
  
### Recommendations and extras

1. Do QC check of the sequencing reads and do quality trimming before mapping and dispose of very short reads. 

2. If the DNA fragmentation during the library construction was random, it will be best to remove PCR duplicates after step 3 before calling samtools mpileup.

3. Check the fragment (insert) size distribution and decide on overlap removal method for paired-end reads. The simplest option is *-x* added to the samtools mpileup command, but we prefered bamUtil tool clipOverlap before calling variants.  
  ```
  clipOverlap -in TAPS_mESC_LOO1_OT/OB.bam --out TAPS_mESC_LOO1_OT/OB_clipped.bam
  ```
4. A script enabling you to split your bam files by chromosome is MF_chrom_split.py, a handy option that can enable you to run scipts in parallel. We split by chromosome before using samtools mpileup followed by MF_callerMOD.py:  
  4.1.  ```
        samtools index TAPS_mESC_LOO1_OT/OB.bam
        ```  
  4.2. ```
      python3 MF_chrom_split.py --input_file TAPS_mESC_LOO1_OT/OB.bam --fasta_file /dir/mm9.fa
       ```
5. Some more scripts can give your fuller information about the sequencing calling. Currently in this category is MF_phredder.py that outputs the average quality per read (if available) for each of the four bases T, C, A, G.  
  ```
  python3 MF_phredder.py --input_file TAPS_mESC_L001_R1.fastq.gz
  ```




