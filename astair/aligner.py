#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import re
import sys
import pdb
import click
import logging
import warnings
import subprocess
from datetime import datetime
from distutils.spawn import find_executable


@click.command()
@click.option('fq1', '--fq1', '-1', required=True, help='First in pair (R1) sequencing reads file in fastq.gz format')
@click.option('fq2', '--fq2', '-2', required=True, help='Second in pair (R2) sequencing reads file in fastq.gz format')
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format used for aligning of the sequencing reads.')
@click.option('bwa_path', '--bwa_path', '-bp', required=False, help='The path to BWA for TAPS-like data and to bwameth.py for bisulfite sequencing.')
@click.option('samtools_path', '--samtools_path', '-sp', required=False, help='The path to Samtools.')
@click.option('directory', '--directory', '-d', required=True, help='Output directory to save files.')
@click.option('method', '--method', '-m', required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like).')
@click.option('output_format', '--output_format', '-O', required=False, default='CRAM', type=click.Choice(['CRAM', 'BAM']), help='Specify output format, possible options are BAM and CRAM. The default is CRAM.')
@click.option('minimum_mapping_quality', '--minimum_mapping_quality', '-mq', required=False, type=int, default=1, help='Set the minimum mapping quality for a read to be output to file (Default >=1).')
@click.option('keep_unmapped', '--keep_unmapped', '-u', default=False, is_flag=True, help='Outputs the unmapped reads (Default false).')
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('minimum_seed_length', '--minimum_seed_length', '-k', default=19, type=int, required=False, help='The minimum seed length used for alignment, see BWA manual (the default value is 19).')
@click.option('band_width', '--band_width', '-w', default=100, type=int, required=False, help='The band width for banded alignment, see BWA manual (the default value is 100).')
@click.option('dropoff', '--dropoff', '-D', default=100, type=int, required=False, help='The off-diagonal X-dropoff, see BWA manual (the default value is 100).')
@click.option('internal_seeds', '--internal_seeds', '-r', default=1.5, type=float, required=False, help='Looks for internal seeds inside a seed longer than minimum_seed_length * internal_seeds, see BWA manual (the default value is 1.5).')
@click.option('reseeding_occurence', '--reseeding_occurence', '-y', default=20, type=int, required=False, help='The seed occurrence for the 3rd round seeding, see BWA manual (the default value is 20).')
@click.option('N_skip_seeds', '--N_skip_seeds', '-c', default=500, type=int, required=False, help='Skips seeds with more than the given seed occurrences, see BWA manual (the default value is 500).')
@click.option('drop_chains', '--drop_chains', '-dc', default=0.5, type=float, required=False, help='Drops chains shorter than the specified fraction of the longest overlapping chain, see BWA manual (the default value is 0.5).')
@click.option('discard_chains', '--discard_chains', '-W', default=0, type=int, required=False, help='Discards a chain if seeded bases shorter than the specified value, see BWA manual (the default value is 0).')
@click.option('N_mate_rescues', '--N_mate_rescues', '-mr', default=50, type=int, required=False, help='Performs at most the specified rounds of mate rescues for each read, see BWA manual (the default value is 50).')
@click.option('skip_mate_rescue', '--skip_mate_rescue', '-s', is_flag=True, required=False, help='NB: Does not recommend unless necessary: skips mate rescue in mCtoT mode, see BWA manual. If set, orphan reads (paired reads that are not in a proper pair) will be generated. Ensure ignore_orphans in the caller is set to False.')
@click.option('skip_pairing', '--skip_pairing', '-P', is_flag=True, required=False, help='NB: Does not recommend unless necessary: skips read pairing in mCtoT mode, but does rescue mates unless mate_skipping is also performed, see BWA manual. If set, orphan reads (paired reads that are not in a proper pair) will be generated. Ensure ignore_orphans in the caller is set to False.')
@click.option('match_score', '--match_score', '-A', default=1, type=int, required=False, help='The score for a sequence match, which scales the remaing scoring options, see BWA manual(the default value is 1).')
@click.option('mismatch_penalty', '--mismatch_penalty', '-B', default=4, type=int, required=False, help='The penalty for a mismatch, see BWA manual (the default value is 4).')
@click.option('gap_open_penalty', '--gap_open_penalty', '-o', default='6,6', type=str, required=False, help='The gap open penalties for deletions and insertions, see BWA manual (the default value is 6,6).')
@click.option('gap_extension_penalty', '--gap_extension_penalty', '-E', default='1,1', type=str, required=False, help='The gap extension penalty with a cost size calculated as {-O} + {-E}*k, see BWA manual (the default value is 1,1).')
@click.option('end_clipping_penalty', '--end_clipping_penalty', '-L', default='5,5', type=str, required=False, help='The penalty for 5-prime- and 3-prime-end clipping, see BWA manual (the default value is 5,5).')
@click.option('unpaired_penalty', '--unpaired_penalty', '-U', default=17, type=int, required=False, help='The penalty for an unpaired read pair, see BWA manual (the default value is 17).')
@click.option('read_type', '--read_type', '-x', default='null', type=click.Choice(['null', 'pacbio', 'ont2d', 'intractg']), required=False, help='Changes multiple parameters unless overridden, see BWA manual (the default value is None).')
def align(fq1, fq2, reference, bwa_path, samtools_path, directory, method, output_format, minimum_mapping_quality, keep_unmapped, N_threads, minimum_seed_length, band_width, dropoff,
                 internal_seeds, reseeding_occurence, N_skip_seeds, drop_chains, discard_chains, N_mate_rescues, skip_mate_rescue, skip_pairing, match_score, mismatch_penalty,
                 gap_open_penalty, gap_extension_penalty, end_clipping_penalty, unpaired_penalty, read_type):
    """Align raw reads in fastq format to a reference genome. bwa is required to align TAPS reads, and bwa-meth fif you plan to process BS-seq data."""
    run_alignment(fq1, fq2, reference, bwa_path, samtools_path, directory, method, output_format, minimum_mapping_quality, keep_unmapped, N_threads, minimum_seed_length, band_width, dropoff,
                 internal_seeds, reseeding_occurence, N_skip_seeds, drop_chains, discard_chains, N_mate_rescues, skip_mate_rescue, skip_pairing, match_score, mismatch_penalty,
                 gap_open_penalty, gap_extension_penalty, end_clipping_penalty, unpaired_penalty, read_type)


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

#logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def which_path(bwa_path, samtools_path, method):
    """Discovers the paths to BWA and Samtools."""
    if bwa_path:
        use_bwa = bwa_path
    elif method == 'mCtoT':
        use_bwa = find_executable('bwa')
    elif method == 'CtoT':
        use_bwa = find_executable('bwameth.py')

    if samtools_path:
        use_samtools = samtools_path
    else:
        use_samtools = find_executable('samtools')

    return use_bwa, use_samtools

def check_index(use_bwa, reference, method):
    """Checks if the provided reference is indexed, and creates an index if one is not found."""
    if (os.path.isfile(reference + '.bwt') == False and method == 'mCtoT') \
            or (os.path.isfile(reference + '.bwameth.c2t') == False and method == 'CtoT'):
        build_command = '{} index {}'.format(use_bwa, reference)
        index_fasta = subprocess.Popen(build_command, shell=True)
        index_fasta.wait()


def run_alignment(fq1, fq2, reference, bwa_path, samtools_path, directory, method, output_format, minimum_mapping_quality, keep_unmapped, N_threads, minimum_seed_length, band_width, dropoff,
                 internal_seeds, reseeding_occurence, N_skip_seeds, drop_chains, discard_chains, N_mate_rescues, skip_mate_rescue, skip_pairing, match_score, mismatch_penalty,
                 gap_open_penalty, gap_extension_penalty, end_clipping_penalty, unpaired_penalty, read_type):
    """Aligns the provided pair-end reads to the reference according to the method specified.
    Outputs a sorted and indexed file."""
    name = os.path.splitext(os.path.basename(fq1))[0]
    name = re.sub('_(R1|1).fq', '', name)
    directory = os.path.abspath(directory)
    if list(directory)[-1] != "/":
        directory = directory + "/"
    if keep_unmapped:
        aligned_string = ''
    else:
        aligned_string = '-F 4 '
    use_bwa, use_samtools = which_path(bwa_path, samtools_path, method)
    if output_format == 'CRAM':
        output_f = '-hC'
    else:
        output_f = '-hb'
    check_index(use_bwa, reference, method)
    if skip_mate_rescue:
        skip_mate_rescue = '-S'
    else:
        skip_mate_rescue = ''
    if skip_pairing:
        skip_pairing = '-P'
    else:
        skip_pairing = ''
    if read_type == 'null':
        read_type = ''
    elif read_type == 'pacbio':
        read_type = '-x pacbio'
    elif read_type == 'ont2d':
        read_type = '-x ont2d'
    elif read_type == 'intractg':
        read_type = '-x intractg'
    if method == 'mCtoT':
        alignment_command = '{} {} -t {} -k {} -w {} -d {} -r {} -y {} -c {} -D {} -W {} -m {} {} ' \
                            ' {} -A {} -B {} -O {} -E {} -L {} -U {} {} {} {} {} | {} view {} -T {} -q {} {} -O {} ' \
                            '| {} sort -@ {} -O {} > {}'.format(use_bwa, 'mem', N_threads, minimum_seed_length, band_width, dropoff,
                 internal_seeds, reseeding_occurence, N_skip_seeds, drop_chains, discard_chains, N_mate_rescues, skip_mate_rescue, skip_pairing, match_score, mismatch_penalty,
                 gap_open_penalty, gap_extension_penalty, end_clipping_penalty, unpaired_penalty, read_type,
                   reference, fq1, fq2, use_samtools, output_f, reference, minimum_mapping_quality,
                   aligned_string, output_format, use_samtools, N_threads, output_format, os.path.join(directory + name + '_' + method + "." + output_format.lower()))
    else:
        alignment_command = 'python2 {} -t {} --reference {} {} {} | {} view {} -T {} -q {} {} -O {} | {} sort -@ {} -O {} > {}'.\
            format(use_bwa, N_threads, reference, fq1, fq2, use_samtools, output_f, reference, minimum_mapping_quality,
                   aligned_string, output_format, use_samtools, N_threads, output_format, os.path.join(directory + name + '_' + method + "." + output_format.lower()))
    try:
        if os.path.isfile(os.path.join(directory + name + '_' + method + "." + output_format.lower())):
            logs.error('The output files will not be overwritten. Please rename the input or the existing output files before rerunning if the input is different.')
            sys.exit(1)
        else:
            align = subprocess.Popen(alignment_command, shell=True)
            exit_code = align.wait()
            if exit_code == 0:
                indexing_command = '{} index {}'.format(use_samtools, os.path.join(directory + name + '_' + method + "." + output_format.lower()))
                index = subprocess.Popen(indexing_command, shell=True)
                index.wait()
        time_e = datetime.now()
        logs.info("asTair genome aligner finished running. {} seconds".format((time_e - time_b).total_seconds()))
    except IOError:
        logs.error('asTair cannot write to alignment file and its index.', exc_info=True)


if __name__ == '__main__':
    align()
