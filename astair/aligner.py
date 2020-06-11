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

from astair.simple_fasta_parser import fasta_splitting_by_sequence


@click.command()
@click.option('fq1', '--fq1', '-1', required=True, help='First in pair (R1) sequencing reads file in fastq.gz format.')
@click.option('fq2', '--fq2', '-2', required=False, help='Second in pair (R2) sequencing reads file in fastq.gz format.')
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format used for aligning of the sequencing reads.')
@click.option('bwa_path', '--bwa_path', '-bp', required=False, help='The path to BWA for TAPS-like data and to bwameth.py for bisulfite sequencing.')
@click.option('samtools_path', '--samtools_path', '-sp', required=False, help='The path to Samtools.')
@click.option('directory', '--directory', '-d', required=True, help='Output directory to save files.')
@click.option('method', '--method', '-m', required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT)')
@click.option('output_format', '--output_format', '-O', required=False, default='BAM', type=click.Choice(['CRAM', 'BAM']), help='Specify output format, possible options are BAM and CRAM. The Default BAM.')
@click.option('minimum_mapping_quality', '--minimum_mapping_quality', '-mq', required=False, type=int, default=1, help='Set the minimum mapping quality for a read to be output to file (Default >=1).')
@click.option('keep_unmapped', '--keep_unmapped', '-u', default=False, is_flag=True, help='Outputs the unmapped reads (Default false).')
@click.option('single_end', '--se', '-se', default=False, is_flag=True, required=False, help='Indicates single-end sequencing reads (Default False).')
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('minimum_seed_length', '--minimum_seed_length', '-k', default=19, type=int, required=False, help='The minimum seed length used for alignment, see BWA manual. (Default 19).')
@click.option('band_width', '--band_width', '-w', default=100, type=int, required=False, help='The band width for banded alignment, see BWA manual. (Default 100).')
@click.option('dropoff', '--dropoff', '-D', default=100, type=int, required=False, help='The off-diagonal X-dropoff, see BWA manual. (Default 100).')
@click.option('internal_seeds', '--internal_seeds', '-r', default=1.5, type=float, required=False, help='Looks for internal seeds inside a seed longer than minimum_seed_length * internal_seeds, see BWA manual. (Default 1.5).')
@click.option('reseeding_occurence', '--reseeding_occurence', '-y', default=20, type=int, required=False, help='The seed occurrence for the 3rd round seeding, see BWA manual. (Default 20).')
@click.option('N_skip_seeds', '--N_skip_seeds', '-c', default=500, type=int, required=False, help='Skips seeds with more than the given seed occurrences, see BWA manual. (Default 500).')
@click.option('drop_chains', '--drop_chains', '-dc', default=0.5, type=float, required=False, help='Drops chains shorter than the specified fraction of the longest overlapping chain, see BWA manual. (Default 0.5).')
@click.option('discard_chains', '--discard_chains', '-W', default=0, type=int, required=False, help='Discards a chain if seeded bases shorter than the specified value, see BWA manual. (Default 0).')
@click.option('N_mate_rescues', '--N_mate_rescues', '-mr', default=50, type=int, required=False, help='Performs at most the specified rounds of mate rescues for each read, see BWA manual. (Default 50).')
@click.option('skip_mate_rescue', '--skip_mate_rescue', '-s', is_flag=True, required=False, help='NB: Does not recommend unless necessary: skips mate rescue in mCtoT mode, see BWA manual. If set, orphan reads (paired reads that are not in a proper pair) will be generated. Ensure ignore_orphans in the caller is set to False.')
@click.option('skip_pairing', '--skip_pairing', '-P', is_flag=True, required=False, help='NB: Does not recommend unless necessary: skips read pairing in mCtoT mode, but does rescue mates unless mate_skipping is also performed, see BWA manual. If set, orphan reads (paired reads that are not in a proper pair) will be generated. Ensure ignore_orphans in the caller is set to False.')
@click.option('match_score', '--match_score', '-A', default=1, type=int, required=False, help='The score for a sequence match, which scales the remaing scoring options, see BWA manual. (Default 1).')
@click.option('mismatch_penalty', '--mismatch_penalty', '-B', default=4, type=int, required=False, help='The penalty for a mismatch, see BWA manual. (Default 4).')
@click.option('gap_open_penalty', '--gap_open_penalty', '-o', default='6,6', type=str, required=False, help='The gap open penalties for deletions and insertions, see BWA manual. (Default 6,6).')
@click.option('gap_extension_penalty', '--gap_extension_penalty', '-E', default='1,1', type=str, required=False, help='The gap extension penalty with a cost size calculated as {-O} + {-E}*k, see BWA manual. (Default 1,1).')
@click.option('end_clipping_penalty', '--end_clipping_penalty', '-L', default='5,5', type=str, required=False, help='The penalty for 5-prime- and 3-prime-end clipping, see BWA manual. (Default 5,5).')
@click.option('unpaired_penalty', '--unpaired_penalty', '-U', default=17, type=int, required=False, help='The penalty for an unpaired read pair, see BWA manual. (Default 17).')
@click.option('read_type', '--read_type', '-x', default='null', type=click.Choice(['null', 'pacbio', 'ont2d', 'intractg']), required=False, help='Changes multiple parameters unless overridden, see BWA manual. (Default null).')
@click.option('smart_pairing', '--smart_pairing', '-smp', default=False, is_flag=True,  required=False, help='Ignores read2, see BWA manual.')
@click.option('read_group', '--read_group', '-rg', default=None, type=str,  required=False, help='Adds the given read group line "@RG\\tID:ids\\tSM:name", see BWA manual. (Default None).')
@click.option('header_string', '--header_string', '-hs', default=None, type=str,  required=False, help='Adds the given string to header if it starts with @, or to file if it does not, see BWA manual.')
@click.option('include_alt', '--include_alt', '-al', default=False, is_flag=True,  required=False, help='Treats ALT contigs as part of the primary assembly, see BWA manual.')
@click.option('split_alignment', '--split_alignment', '-spa', default=False, is_flag=True,  required=False, help='Takes the alignment with the smallest coordinate as primary in case of split alignment, see BWA manual.')
@click.option('supplementary_mapq', '--supplementary_mapq', '-smq', default=False, is_flag=True,  required=False, help='Does not modify the MAPQ of supplementary alignments, see BWA manual.')
@click.option('minimum_score', '--minimum_score', '-ms', default=30, type=int, required=False, help='Minimum score to ouput, see BWA manual. (Default 30).')
@click.option('alternative_score', '--alternative_score', '-h', default='5,200', type=str, required=False, help='Defines the number of reads that will be tagged as XA if they have 80 percent of the maximum score, see BWA manual. (Default 5,200).')
@click.option('all_alignments', '--all_alignments', '-aa', default=False, is_flag=True,  required=False, help='Outputs all SE or unpaired PE reads, see BWA manual.')    
@click.option('fasta_comment', '--fasta_comment', '-fc', default=None,  required=False, help='Adds FASTA/FASTQ comment to output, see BWA manual.')   
@click.option('fasta_header',  '--fasta_header', '-fh', default=None,  required=False, help='Outputs the reference header in the XR tag, see BWA manual.')   
@click.option('clip_supplementary',  '--clip_supplementary', '-cs', default=False, is_flag=True,  required=False, help='Enables soft clipping on supplementary alignments, see BWA manual.')   
@click.option('mark_splitted',  '--mark_splitted', '-ms', default=False, is_flag=True,  required=False, help='Labels short split reads as secondary, see BWA manual.')                                                               
@click.option('reads_distribution', '--reads_distribution', '-rd', default=None, required=False, help='Specifies the mean, standard deviation (10 percent of the mean if absent), max (4 sigma from the mean if absent) and min of the insert size distribution for FR orientation, see BWA manual. Must be provided as FLOAT,FLOAT,INT,INT (Default reads distribution metrics are inferred from the data).')
@click.option('add_underscores', '--add_underscores', '-au', default=False, is_flag=True, required=False, help='Indicates outputting a new reference fasta file with added underscores in the sequence names that is afterwards used for calling. (Default False).')
@click.option('use_underscores', '--use_underscores', '-uu', default=False, is_flag=True, required=False, help='Uses as a reference the fasta file with added underscores in the sequence names that is afterwards used for calling. (Default False).')
@click.option('temp_dir', '--temp_dir', '-td', required=False, help='Provides a custom directory to write temporary files. (Default the chosen directory for the output).')
@click.option('compress', '--compress', '-z', is_flag=True, default=False, required=False, help='Should the reference FASTA be compressed after the run (Default False).')
@click.option('sort_chunck_size', '--sort_chunck_size', '-cz', default='768M', required=False, help='WARNING: sorting large files >= 50 GB might require to increaase the maximum memory per thread parameter; recognised sufixes are  K/M/G (Default 768M).')
def align(fq1, fq2, reference, bwa_path, samtools_path, directory, method, output_format, minimum_mapping_quality, keep_unmapped, N_threads, minimum_seed_length, band_width, dropoff,
                 internal_seeds, reseeding_occurence, N_skip_seeds, drop_chains, discard_chains, N_mate_rescues, skip_mate_rescue, skip_pairing, match_score, mismatch_penalty,
                 gap_open_penalty, gap_extension_penalty, end_clipping_penalty, unpaired_penalty, read_type, single_end, smart_pairing, read_group, header_string, include_alt, split_alignment, supplementary_mapq, minimum_score, alternative_score, all_alignments, fasta_comment, fasta_header, clip_supplementary, mark_splitted, reads_distribution, add_underscores, temp_dir, use_underscores, compress, sort_chunck_size):
    """Align raw reads in fastq format to a reference genome. bwa is required to align TAPS reads, and bwa-meth fif you plan to process BS-seq data."""
    run_alignment(fq1, fq2, reference, bwa_path, samtools_path, directory, method, output_format, minimum_mapping_quality, keep_unmapped, N_threads, minimum_seed_length, band_width, dropoff,
                 internal_seeds, reseeding_occurence, N_skip_seeds, drop_chains, discard_chains, N_mate_rescues, skip_mate_rescue, skip_pairing, match_score, mismatch_penalty,
                 gap_open_penalty, gap_extension_penalty, end_clipping_penalty, unpaired_penalty, read_type, single_end, smart_pairing, read_group, header_string, include_alt, split_alignment, supplementary_mapq, minimum_score, alternative_score, all_alignments, fasta_comment, fasta_header, clip_supplementary, mark_splitted, reads_distribution, add_underscores,temp_dir, use_underscores, compress, sort_chunck_size)


warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

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


def check_index(use_bwa, reference, method, output_format, add_underscores, use_underscores):
    """Checks whether there are spaces in the reference names in the fasta file. In case such spaces exist, they will be replaced with underscores before building the index. Otherwise, checks if the provided reference is indexed, and creates an index if one is not found."""
    fasta_splitting_by_sequence(reference, None, 'last', add_underscores, 'all')
    reference_base_name = os.path.splitext(os.path.basename(reference))[0]
    reference_extension = os.path.splitext(os.path.basename(reference))[1]
    reference_dir = os.path.dirname(reference)
    if list(reference_dir)[-1]!="/":
        reference_dir = reference_dir + "/"
    if use_underscores and ((os.path.isfile(reference_dir + os.path.splitext(os.path.basename(reference_base_name))[0] + '_no_spaces.fa.gz') and reference_extension == '.gz') or (os.path.isfile(reference_dir + reference_base_name + '_no_spaces.fa.gz') and reference_extension != '.gz')):
        if reference_extension == '.gz':
            reference = reference_dir + os.path.splitext(os.path.basename(reference_base_name))[0]  + '_no_spaces.fa.gz'
            reference_base_name = os.path.splitext(os.path.basename(reference_base_name))[0] + '_no_spaces.fa'
        else:
            reference = reference_dir + reference_base_name + '_no_spaces.fa.gz'
            reference_base_name = reference_base_name + '_no_spaces.fa'
    try:
        if reference_extension == '.gz':
            bgzip_ = subprocess.Popen('bgzip -r {}'.format(reference), shell=True)
            exit_code = bgzip_.wait()
            if exit_code == 0:
                compressed_ = 'bgzip'
            else:
                compressed_ = 'gzip'
    except Exception:
        logs.info('No bgzip was found to create gzi index. Bgzip indexing may improve performance.', exc_info=True)
        compressed_ = 'gzip'
    if (os.path.isfile(reference + '.bwt') == False and method == 'mCtoT') \
            or (os.path.isfile(reference + '.bwameth.c2t') == False and method == 'CtoT'):
        if output_format == 'CRAM' and reference_extension == '.gz' and compressed_ == 'gzip':
            subprocess.Popen('gunzip {}'.format(reference), shell=True)
            reference = reference_dir + reference_base_name
        build_command = '{} index {}'.format(use_bwa, reference)
        index_fasta = subprocess.Popen(build_command, shell=True)
        index_fasta.wait()
    if isinstance(reference, tuple):
        reference = reference[0] + reference[1]
    return reference


def run_alignment(fq1, fq2, reference, bwa_path, samtools_path, directory, method, output_format, minimum_mapping_quality, keep_unmapped, N_threads, minimum_seed_length, band_width, dropoff,
                 internal_seeds, reseeding_occurence, N_skip_seeds, drop_chains, discard_chains, N_mate_rescues, skip_mate_rescue, skip_pairing, match_score, mismatch_penalty,
                 gap_open_penalty, gap_extension_penalty, end_clipping_penalty, unpaired_penalty, read_type, single_end, smart_pairing, read_group, header_string, include_alt, split_alignment, supplementary_mapq, minimum_score, alternative_score, all_alignments, fasta_comment, fasta_header, clip_supplementary, mark_splitted, reads_distribution, add_underscores, temp_dir, use_underscores, compress, sort_chunck_size):
    """Aligns the provided pair-end reads to the reference according to the method specified.
    Outputs a sorted and indexed file."""
    time_s = datetime.now()
    logs.info("asTair genome aligner started running. {} seconds".format((time_s - time_b).total_seconds()))
    name = os.path.splitext(os.path.basename(fq1))[0]
    if  os.path.splitext(os.path.basename(reference))[1] == '.gz':
        initial_compression = True
    else:
        initial_compression = False
    if single_end == False  or fq2 is not None:
        name = re.sub('(_R1|_1)', '', name)
    if  os.path.splitext(os.path.basename(fq1))[1] == '.gz':
        name = "_".join(name.split('.')[:-1])
    directory = os.path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    if os.path.exists(directory) == False:
        raise Exception("The output directory does not exist.")
        sys.exit(1)
    if keep_unmapped:
        aligned_string = ''
        minimum_mapping_quality = 0
    else:
        aligned_string = '-F 4 '
    if use_underscores:
        add_underscores = True
    use_bwa, use_samtools = which_path(bwa_path, samtools_path, method)
    reference = check_index(use_bwa, reference, method, output_format, add_underscores, use_underscores)
    if skip_mate_rescue:
        skip_mate_rescue = '-S'
    else:
        skip_mate_rescue = ''
    if temp_dir:
        if list(temp_dir)[-1]!="/":
            temp_dir = temp_dir + "/"
        temp_name = temp_dir + 'temp'
    else:
        temp_name = directory + 'temp'
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
    if single_end == True or fq2 is None:
        fq2 = ''
    if smart_pairing:
        smart_pairing = ' -p '
    else:
        smart_pairing = ''
    if read_group is None:
        read_group = ''
    else:
        read_group = ' -R ' + read_group
    if header_string is None:
        header_string = ''
    else:
        header_string = ' -H ' + header_string
    if include_alt:
        include_alt = ' -j '
    else:
        include_alt = ''
    if split_alignment:
        split_alignment = ' -5 '
    else:
        split_alignment = ''
    if supplementary_mapq:
        supplementary_mapq = ' -q ' 
    else:
        supplementary_mapq = ''
    if all_alignments:
        all_alignments = ' -a '
    else:
        all_alignments = ''
    if fasta_comment is None:
        fasta_comment = ''
    else:
        fasta_comment = ' -C ' + fasta_comment
    if fasta_header is None:
        fasta_header = ''
    else:
        fasta_header = ' -V ' + fasta_header
    if clip_supplementary:
        clip_supplementary = ' -Y '
    else:
        clip_supplementary = ''
    if mark_splitted:
        mark_splitted = ' -M '
    else:
        mark_splitted = ''
    if reads_distribution is None:
        reads_distribution = ''
    else:
        reads_distribution = ' -I '+ reads_distribution
    if output_format == 'CRAM':
        output_f = '-hC'
    else:
        output_f = '-hb'
    if method == 'mCtoT':
        alignment_command = '{} {} -t {} -k {} -w {} -d {} -r {} -y {} -c {} -D {} -W {} -m {} {} ' \
                            ' {} -A {} -B {} -O {} -E {} -L {} -U {} {} {} {} {} {} {} {} -T {} -h {} {} {} {} {} {} {} {} {} {} | {} view {} -T {} -q {} {} -O {} ' \
                            '| {} sort -m {} -@ {} -O {} > {}'.format(use_bwa, 'mem', N_threads, minimum_seed_length, band_width, dropoff,
                 internal_seeds, reseeding_occurence, N_skip_seeds, drop_chains, discard_chains, N_mate_rescues, skip_mate_rescue, skip_pairing, match_score, mismatch_penalty,
                 gap_open_penalty, gap_extension_penalty, end_clipping_penalty, unpaired_penalty, read_type, smart_pairing, read_group, header_string, include_alt, split_alignment, supplementary_mapq, minimum_score, alternative_score, all_alignments, fasta_comment, fasta_header, clip_supplementary, mark_splitted, reads_distribution, reference, fq1, fq2, use_samtools, output_f, reference, minimum_mapping_quality,
                   aligned_string, output_format, use_samtools, sort_chunck_size, N_threads, output_format, os.path.join(directory + name + '_' + method + "." + output_format.lower()))
    else:
        alignment_command = 'python {} -t {} --reference {} {} {} | {} view {} -T {} -q {} {} -O {} | {} sort -T {} -@ {} -O {} > {}'.\
            format(use_bwa, N_threads, reference, fq1, fq2, use_samtools, output_f, reference, minimum_mapping_quality,
                   aligned_string, output_format, use_samtools, temp_name, N_threads, output_format, os.path.join(directory + name + '_' + method + "." + output_format.lower()))
    try:
        if os.path.isfile(os.path.join(directory + name + '_' + method + "." + output_format.lower())):
            logs.error('The output files will not be overwritten. Please rename the input or the existing output files before rerunning if the input is different.')
        else:
            align = subprocess.Popen(alignment_command, shell=True)
            exit_code = align.wait()
            if exit_code == 0:
                indexing_command = '{} index {}'.format(use_samtools, os.path.join(directory + name + '_' + method + "." + output_format.lower()))
                index = subprocess.Popen(indexing_command, shell=True)
                index.wait()
        if os.path.splitext(os.path.basename(reference))[1] != '.gz' and (compress == True or initial_compression==True):
            try:
                gzip_ = subprocess.Popen('bgzip {}'.format(reference), shell=True)
            except Exception:
                logs.error('The output fasta will be GZIP compressed.')
                gzip_ = subprocess.Popen('gzip {}'.format(reference), shell=True)
            gzip_.wait()
        time_e = datetime.now()
        logs.info("asTair genome aligner finished running. {} seconds".format((time_e - time_b).total_seconds()))
    except IOError:
        logs.error('asTair cannot write to alignment file and its index.', exc_info=True)


if __name__ == '__main__':
    align()

