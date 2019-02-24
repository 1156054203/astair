#!/usr/bin/env python3

import click
import re
from datetime import datetime
import logging
import warnings
from distutils.spawn import find_executable
import os
import subprocess


@click.command()
@click.option('fq1', '--fq1', '-1', required=True, help='First in pair (R1) sequencing reads file in fastq.gz format')
@click.option('fq2', '--fq2', '-2', required=True, help='Second in pair (R2) sequencing reads file in fastq.gz format')
@click.option('fasta_file', '--fasta_file', '-f', required=True, help='DNA sequence in fasta format used for aligning the sequencing reads.')
@click.option('bwa_path', '--bwa_path', '-bp', required=False, help='The path to BWA.')
@click.option('samtools_path', '--samtools_path', '-sp', required=False, help='The path to Samtools.')
@click.option('directory', '--directory', '-d', required=True, help='Output directory to save files.')
@click.option('method', '--method', '-m', required=False, default='CmtoT', type=click.Choice(['CtoT', 'CmtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and CmtoT (modified cytosines are converted to thymines, TAPS-like).')
@click.option('minimum_mapping_quality', '--minimum_mapping_quality', '-mq', required=False, type=int, default=1, help='Set the minimum mapping quality for a read to be output to file (Default >=1).')
@click.option('keep_unmapped', '--keep_unmapped', '-u', default=False, is_flag=True, help='Outputs the unmapped reads (Default false).')
@click.option('N_threads', '--N_threads', '-t', default=10, required=True, help='The number of threads to spawn (the default value is 10).')
def aligner_exec(fq1, fq2, fasta_file, bwa_path, samtools_path, directory, method, minimum_mapping_quality, keep_unmapped, N_threads):
    run_alignment(fq1, fq2, fasta_file, bwa_path, samtools_path, directory, method, minimum_mapping_quality, keep_unmapped, N_threads)


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def which_path(bwa_path, samtools_path):
    if bwa_path:
        use_bwa = bwa_path
    else:
        use_bwa = find_executable('bwa')

    if samtools_path:
        use_samtools = samtools_path
    else:
        use_samtools = find_executable('samtools')

    return use_bwa, use_samtools

def check_index(use_bwa, fasta_file):
    if os.path.isfile(fasta_file + '.bwt') == False:
        build_command = '{} index {}'.format(use_bwa, fasta_file)
        subprocess.Popen(build_command, shell=True)


def run_alignment(fq1, fq2, fasta_file, bwa_path, samtools_path, directory, method, minimum_mapping_quality, keep_unmapped, N_threads):
    name = os.path.splitext(os.path.basename(fq1))[0]
    name = re.sub('_(R1|1).fq', '', name)
    directory = os.path.abspath(directory)
    if list(directory)[-1] != "/":
        directory = directory + "/"
    use_bwa, use_samtools = which_path(bwa_path, samtools_path)
    if method == 'CmtoT':
        align = 'mem'
    else:
        align = 'meth'
    check_index(use_bwa, fasta_file)
    if keep_unmapped:
        aligned_string = ''
    else:
        aligned_string = '-F 4 '
    alignment_command = '{} {} -t {} {} {} {} | {} view -hC -T {} -q {} {} | {} sort -@ {} > {}'.format(use_bwa, align, N_threads,
    fasta_file, fq1, fq2, use_samtools, fasta_file, minimum_mapping_quality, aligned_string, use_samtools, N_threads, os.path.join(directory + name + ".cram"))
    align = subprocess.Popen(alignment_command, shell=True)
    exit_code = align.wait()
    if exit_code:
        indexing_command = '{} index {}'.format( use_samtools,os.path.join(directory + name + ".cram"))
        index = subprocess.Popen(indexing_command, shell=True)
        index.wait()
    time_e = datetime.now()
    logs.info("asTair genome aligner finished running. {} seconds".format((time_e - time_b).total_seconds()))

if __name__ == '__main__':
    aligner_exec()
