#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import re
import os
import sys
import pdb
import csv
import gzip
import click
import pysam
import shutil
import logging
import warnings
import subprocess
from os import path
from datetime import datetime
from collections import defaultdict

if sys.version[0] == '3':
    from itertools import zip_longest
elif sys.version[0] == '2':
    from itertools import izip_longest as zip_longest
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))

from astair.safe_division import non_zero_division
from astair.bam_file_parser import bam_file_opener
from astair.context_search import context_sequence_search
from astair.context_search import sequence_context_set_creation
from astair.simple_fasta_parser import fasta_splitting_by_sequence

@click.command()
@click.option('input_file', '--input_file', '-i', required=True, help='BAM|CRAM format file containing sequencing reads.')
# @click.option('control_file', '--control_file', '-c', required=False, help='BAM format file containing sequencing reads of a matched control.')
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format used for aligning of the sequencing reads and for pileup.')
@click.option('zero_coverage', '--zero_coverage', '-z', default=False, is_flag=True, help='When set to True, outputs positions not covered in the bam file. Uncovering zero coverage positions takes longer time than using the default option.')
@click.option('context', '--context', '-co', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis. (Default all).')
@click.option('user_defined_context', '--user_defined_context', '-uc', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
# @click.option('library', '--library', '-l', required=False, default = 'directional',  type=click.Choice(['directional', 'rr', 'ptat']), help='Provides information for the library preparation protocol (rr is reduced representation and ptat is post-method adapter-tagging).')
@click.option('method', '--method', '-m', required=False, default = 'mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT).')
@click.option('skip_clip_overlap', '--skip_clip_overlap', '-sc', required=False, default=False, type=bool, help='Skipping the random removal of overlapping bases between pair-end reads. Not recommended for pair-end libraries, unless the overlaps are removed prior to calling. (Default False)')
@click.option('single_end', '--se', '-se', default=False, is_flag=True, required=False, help='Indicates single-end sequencing reads (Default False).')
@click.option('minimum_base_quality', '--minimum_base_quality', '-bq', required=False, type=int, default=20, help='Set the minimum base quality for a read base to be used in the pileup (Default 20).')
@click.option('minimum_mapping_quality', '--minimum_mapping_quality', '-mq', required=False, type=int, default=0, help='Set the minimum mapping quality for a read to be used in the pileup (Default 0).')
@click.option('adjust_acapq_threshold', '--adjust_capq_threshold', '-amq', required=False, type=int, default=0, help='Used to adjust the mapping quality with default 0 for no adjustment and a recommended value for adjustment 50. (Default 0).')
@click.option('mark_matches', '--mark_matches', '-mm', required=False, default=True, type=bool, help='Output bases matching the reference per strand (Default True).')
@click.option('mark_ends', '--mark_ends', '-me', required=False, default=True, type=bool, help='Marks head and tail bases in the read (Default True).')
@click.option('add_indels', '--add_indels', '-ai', required=False, default=True, type=bool, help='Adds inserted bases and Ns for base skipped from the reference (Default True).')
@click.option('redo_baq', '--redo_baq', '-rbq', required=False, default=False, type=bool, help='Re-calculates per-Base Alignment Qualities ignoring existing base qualities (Default False).')
@click.option('compute_baq', '--compute_baq', '-cbq', required=False, default=True, type=bool, help='Performs re-alignment computing of per-Base Alignment Qualities (Default True).')
@click.option('ignore_orphans', '--ignore_orphans', '-io', required=False, default=True, type=bool, help='Ignore reads not in proper pairs (Default True).')
@click.option('max_depth', '--max_depth', '-md', required=False, type=int, default=250, help='Set the maximum read depth for the pileup, maximum value 8000 (Default 250).')
@click.option('per_chromosome', '--per_chromosome', '-chr', default=None, type=str, help='When used, it calculates the modification rates only per the chromosome given. (Default None).')
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('compress', '--gz', '-z', default=False, is_flag=True, required=False, help='Indicates whether the mods file output will be compressed with gzip (Default False).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
def call(input_file, reference, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, method, minimum_mapping_quality, adjust_acapq_threshold,mark_matches, mark_ends, add_indels, redo_baq, compute_baq, ignore_orphans, max_depth,per_chromosome, N_threads, directory, compress, single_end):
        """Call modified cytosines from a bam or cram file. The output consists of two files, one containing modification counts per nucleotide, the other providing genome-wide statistics per context."""
        cytosine_modification_finder(input_file, reference, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, method, minimum_mapping_quality, adjust_acapq_threshold,mark_matches, mark_ends, add_indels, redo_baq, compute_baq, ignore_orphans, max_depth, per_chromosome, N_threads, directory, compress, single_end)


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

#logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def modification_calls_writer(data_mods, compress, data_line, header=False):
    """Outputs the modification calls per position in a tab-delimited format."""
    try:
        if compress == False:
            if header:
                data_line.writerow(["CHROM", "START", "END", "MOD_LEVEL", "MOD", "UNMOD", "REF", "ALT", "SPECIFIC_CONTEXT", "CONTEXT", 'SNV', 'TOTAL_DEPTH'])
                data_line.writerow(data_mods)
            else:
                data_line.writerow(data_mods)
        else:
            if header:
                data_line.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('CHROM', 'START', 'END', 'MOD_LEVEL', 'MOD', 'UNMOD', 'REF', 'ALT', 'SPECIFIC_CONTEXT', 'CONTEXT', 'SNV', 'TOTAL_DEPTH'))
                data_line.write(
                    '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(data_mods[0], data_mods[1], data_mods[2], data_mods[3],
                                                                 data_mods[4], data_mods[5], data_mods[6], data_mods[7],
                                                                 data_mods[8], data_mods[9], data_mods[10], data_mods[11]))
            else:
                data_line.write(
                    '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(data_mods[0], data_mods[1], data_mods[2], data_mods[3],
                                                                 data_mods[4], data_mods[5], data_mods[6], data_mods[7],
                                                                 data_mods[8], data_mods[9], data_mods[10], data_mods[11]))
    except IOError:
        logs.error('asTair cannot write to modification calls file.', exc_info=True)


def statistics_calculator(mean_mod, mean_unmod, data_mod, user_defined_context, context_sample_counts):
    """Calculates the summary statistics of the cytosine modificaton levels."""
    context_sample_counts[data_mod[9]] += 1
    context_sample_counts[data_mod[8]] += 1
    for context in list(('CpG', 'CHH', 'CHG')):
        if re.match(context, data_mod[9]) and data_mod[10] == 'No':
            mean_mod[context] += data_mod[4]
            mean_unmod[context] += data_mod[5]
    for context in list(('CAG', 'CCG', 'CTG', 'CTT', 'CCT', 'CAT', 'CTA', 'CTC', 'CAC', 'CAA', 'CCA', 'CCC', 'CGA', 'CGT', 'CGC', 'CGG')):
        if re.match(context, data_mod[8]) and data_mod[10] == 'No':
            mean_mod[context] += data_mod[4]
            mean_unmod[context] += data_mod[5]
    if re.match(r"CN", data_mod[9]) and data_mod[10] == 'No':
        mean_mod['CNN'] += data_mod[4]
        mean_unmod['CNN'] += data_mod[5]
    if user_defined_context and re.match('user defined context', data_mod[9]) and data_mod[10] == 'No':
        mean_mod['user defined context'] += data_mod[4]
        mean_unmod['user defined context'] += data_mod[5]


def context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context, total_contexts, sub_contexts, header):
    """Writes the summary statistics of the cytosine modificaton levels."""
    with open(file_name, 'a') as statistics_output:
            write_file = csv.writer(statistics_output, delimiter='\t', lineterminator='\n')
            if header == True:
                write_file.writerow(["CONTEXT", "SPECIFIC_CONTEXT", "MEAN_MODIFICATION_RATE_PERCENT", "TOTAL_POSITIONS", "COVERED_POSITIONS", 'MODIFIED', 'UNMODIFIED'])
            write_file.writerow([context, "*", round(non_zero_division(mean_mod[context], mean_mod[context] + mean_unmod[context]) * 100, 3),
                        context_total_counts[total_contexts]+context_total_counts[total_contexts + 'b'], context_sample_counts[context], mean_mod[context], mean_unmod[context]])
            if len(sub_contexts) >= 1:
                for subcontext in sub_contexts:
                    write_file.writerow(["*", subcontext, round(non_zero_division(mean_mod[subcontext], mean_mod[subcontext] + mean_unmod[subcontext]) * 100, 3), context_total_counts[subcontext], context_sample_counts[subcontext], mean_mod[subcontext], mean_unmod[subcontext]])
            if user_defined_context:
                wr.writerow([user_defined_context, "*", round(non_zero_division(mean_mod['user defined context'], mean_mod['user defined context'] + mean_unmod['user defined context']) * 100, 3), context_total_counts['user defined context'], context_sample_counts['user defined context'], mean_mod['user defined context'], mean_unmod['user defined context']])
            

def final_statistics_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context):
    """Writes the summary statistics of the cytosine modificaton levels.
    Cytosine modification rate given as the percentage total modified cytosines
    divided by the total number of cytosines covered."""
    try:
        if context == 'all':
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, 'CpG', 'CG', list(('CGA','CGC', 'CGG', 'CGT')), True)
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, 'CHG', 'CHG', list(('CAG','CCG', 'CTG')), False)
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, 'CHH', 'CHH', list(('CTT', 'CAT', 'CCT', 'CTA', 'CAA', 'CCA', 'CTC', 'CAC', 'CCC')), False)
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, 'CNN', 'CN', list(), False)
        elif context == 'CpG':
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context, 'CG', list(('CGA','CGC', 'CGG', 'CGT')), True)
        elif context == 'CHG':
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context, 'CHG', list(('CAG','CCG', 'CTG')), True)
        elif context == 'CHH':
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context, 'CHH', list(('CTT', 'CAT', 'CCT', 'CTA', 'CAA', 'CCA', 'CTC', 'CAC', 'CCC')), True)
    except IOError:
        logs.error('asTair cannot write to statistics summary file.', exc_info=True)


def flags_expectation(modification_information_per_position, position, ignore_orphans, single_end):
    """Gives the expected flag-base couples, the reference and the modified base."""
    if single_end == True:
        if modification_information_per_position[position][3] == 'C':
            desired_tuples = [(0, 'C'), (0, 'T')]
            undesired_tuples = [(16, 'C'), (16, 'T')]
            modification = 'T'
            reference = 'C'
        elif modification_information_per_position[position][3] == 'G':
            desired_tuples = [(16, 'G'), (16, 'A')]
            undesired_tuples = [(0, 'G'), (0, 'A')]
            modification = 'A'
            reference = 'G'
    else:
        if modification_information_per_position[position][3] == 'C':
            desired_tuples = [(147, 'C'), (99, 'C'), (147, 'T'), (99, 'T')]
            undesired_tuples = [(163, 'C'), (83, 'C'), (163, 'T'), (83, 'T')]
            modification = 'T'
            reference = 'C'
            if ignore_orphans == False:
                desired_tuples.extend([(145, 'C'), (97, 'C'), (145, 'T'), (97, 'T')])
                undesired_tuples.extend([(161, 'C'), (81, 'C'), (161, 'T'), (81, 'T')])
        elif modification_information_per_position[position][3] == 'G':
            desired_tuples = [(163, 'G'), (83, 'G'), (163, 'A'), (83, 'A')]
            undesired_tuples = [(147, 'G'), (99, 'G'), (147, 'A'), (99, 'A')]
            modification = 'A'
            reference = 'G'
            if ignore_orphans == False:
                desired_tuples.extend([(161, 'G'), (81, 'G'), (161, 'A'), (81, 'A')])
                undesired_tuples.extend([(145, 'G'), (97, 'G'), (145, 'A'), (97, 'A')])
    return desired_tuples, undesired_tuples, modification, reference

        
def pileup_summary(modification_information_per_position, position, read_counts, mean_mod, mean_unmod, user_defined_context,
                   header, desired_tuples, undesired_tuples, modification, reference, depth, method, context_sample_counts, ignore_orphans, single_end, compress, data_line):
    """Creates the modification output per position in the format:
    [chrom, start, end, mod_level, mod, unmod, ref, alt, specific_context, context, snv, total_depth] 
    given the strand information and whether the library is pair-end or single-end. The key structure is read_counts
    that contains as dictionary items (read flag, base) tuples.
    Assigns heuristic snv categories of homozygous and not a snv by using the base ratios of the opposite strand."""
    if single_end == True:
        if non_zero_division(read_counts[undesired_tuples[1]], (read_counts[undesired_tuples[0]] + read_counts[undesired_tuples[1]])) < 0.8:
            snp = 'No'
        else:
            snp = 'homozygous'
        if method == 'mCtoT':
            all_data = list((position[0], position[1], position[1] + 1, round(non_zero_division(read_counts[desired_tuples[1]], (read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]])), 3),
                             read_counts[desired_tuples[1]], read_counts[desired_tuples[0]], reference, modification, modification_information_per_position[position][0],
                             modification_information_per_position[position][1], snp, depth))
        elif method =='CtoT':
            all_data = list((position[0], position[1], position[1] + 1, round(non_zero_division(read_counts[desired_tuples[0]], (read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]])), 3),
                             read_counts[desired_tuples[0]], read_counts[desired_tuples[1]], reference, modification, modification_information_per_position[position][0], modification_information_per_position[position][1], snp, depth))
    else:
        if ignore_orphans:
            if non_zero_division(read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]],
                                 (read_counts[undesired_tuples[0]] + read_counts[undesired_tuples[1]]
                                      + read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]])) < 0.8:
                snp = 'No'
            else:
                snp = 'homozygous'
            if method == 'mCtoT':
                all_data = list((position[0], position[1], position[1] + 1, round(
                non_zero_division(read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]], (
                    read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]] + read_counts[desired_tuples[0]] + read_counts[
                        desired_tuples[1]])), 3), read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]],
                             read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]], reference, modification,
                             modification_information_per_position[position][0],
                             modification_information_per_position[position][1], snp, depth))
            elif method == 'CtoT':
                all_data = list((position[0], position[1], position[1] + 1, round(
                    non_zero_division(read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]], (
                        read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]] + read_counts[desired_tuples[0]] + read_counts[
                            desired_tuples[1]])), 3), read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]],
                                 read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]], reference, modification,
                                 modification_information_per_position[position][0],
                                 modification_information_per_position[position][1], snp, depth))
        else:
            if non_zero_division(read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]] + read_counts[undesired_tuples[6]] + read_counts[undesired_tuples[7]],
                                 (read_counts[undesired_tuples[0]] + read_counts[undesired_tuples[1]]
                                      + read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]] + read_counts[undesired_tuples[4]] + read_counts[undesired_tuples[5]]
                                      + read_counts[undesired_tuples[6]] + read_counts[undesired_tuples[7]])) < 0.8:
                snp = 'No'
            else:
                snp = 'homozyguous'
            if method == 'mCtoT':
                all_data = list((position[0], position[1], position[1] + 1, round(
                non_zero_division(read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]] + read_counts[desired_tuples[6]] + read_counts[desired_tuples[7]], (
                    read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]] + read_counts[desired_tuples[0]] + read_counts[
                        desired_tuples[1]] + read_counts[desired_tuples[4]] + read_counts[desired_tuples[5]] + read_counts[desired_tuples[6]] + read_counts[desired_tuples[7]])), 3),
                                 read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]] + read_counts[desired_tuples[6]] + read_counts[desired_tuples[7]],
                             read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]] + read_counts[desired_tuples[4]] + read_counts[desired_tuples[5]], reference, modification,
                             modification_information_per_position[position][0],
                             modification_information_per_position[position][1], snp, depth))
            elif method == 'CtoT':
                all_data = list((position[0], position[1], position[1] + 1, round(
                    non_zero_division(read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]] + read_counts[desired_tuples[4]] + read_counts[desired_tuples[5]], (
                        read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]] + read_counts[desired_tuples[0]] + read_counts[
                            desired_tuples[1]]) + read_counts[desired_tuples[4]] + read_counts[desired_tuples[5]] + read_counts[desired_tuples[6]] +
                                      read_counts[desired_tuples[7]]), 3), read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]] + read_counts[desired_tuples[4]]
                                 + read_counts[desired_tuples[5]], read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]] + read_counts[desired_tuples[6]] + read_counts[desired_tuples[7]],
                                 reference, modification, modification_information_per_position[position][0],
                                 modification_information_per_position[position][1], snp, depth))
    statistics_calculator(mean_mod, mean_unmod, all_data, user_defined_context, context_sample_counts)
    modification_calls_writer(all_data, compress, data_line, header=header)


def clean_pileup(pileups, cycles, modification_information_per_position, mean_mod, mean_unmod, user_defined_context, file_name, method,
                 mark_matches, mark_ends, add_indels, context_sample_counts, ignore_orphans, single_end, compress, data_line):
    """Takes reads from the piled-up region and calculates modification levels."""
    for reads in pileups:
        if cycles == 0:
            header = True
        else:
            header = False
        if (reads.reference_name, reads.pos, reads.pos + 1) in modification_information_per_position:
            position = (reads.reference_name, reads.pos, reads.pos + 1)
            desired_tuples, undesired_tuples, modification, reference = flags_expectation(modification_information_per_position, position, ignore_orphans, single_end)
            read_counts = defaultdict(int)
            try:
                sequences = reads.get_query_sequences(mark_matches=mark_matches, mark_ends=mark_ends, add_indels=add_indels)
            except AssertionError:
                logs.exception("Failed getting query sequences (AssertionError, pysam). Please decrease the max_depth parameter.")
                continue
            for pileup, seq in zip_longest(reads.pileups, sequences, fillvalue='BLANK'):
                if pileup != 'BLANK':
                    read_counts[(pileup.alignment.flag, seq.upper())] += 1
            pileup_summary(modification_information_per_position, position, read_counts, mean_mod, mean_unmod, user_defined_context, header, desired_tuples, undesired_tuples, modification, reference, reads.get_num_aligned(), method, context_sample_counts, ignore_orphans, single_end, compress, data_line)
            modification_information_per_position.pop(position)
            cycles += 1

def cytosine_modification_finder(input_file, reference, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, method,
                                 minimum_mapping_quality, adjust_acapq_threshold, mark_matches, mark_ends, add_indels, redo_baq, compute_baq, ignore_orphans,
                                 max_depth, per_chromosome, N_threads, directory, compress, single_end):
    """Searches for cytosine modification positions in the desired contexts and calculates the modificaton levels."""
    time_s = datetime.now()
    logs.info("asTair modification finder started running. {} seconds".format((time_s - time_b).total_seconds()))
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if os.path.exists(directory) == False:
        raise Exception("The output directory does not exist.")
        sys.exit(1)
    if per_chromosome == None:
        file_name = path.join(directory, name + "_" + method + "_" + context + ".mods")
    else:
        file_name = path.join(directory, name + "_" + method + "_" + per_chromosome + "_" + context + ".mods")
    if not os.path.isfile(file_name) and not os.path.isfile(file_name + '.gz'):
        if user_defined_context:
            mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'CNN': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0,
                        'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'user defined context': 0,
                        'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
            mean_unmod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'CNN': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0,
                          'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'user defined context': 0,
                          'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
        else:
            mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'CNN': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0,
                        'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
            mean_unmod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'CNN': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0,
                          'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
        try:
            inbam = bam_file_opener(input_file, None, N_threads)
        except Exception:
            sys.exit(1)
        try:
            keys, fastas = fasta_splitting_by_sequence(reference, per_chromosome, None)
        except Exception:
            sys.exit(1)
        contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
        cycles = 0
        context_total_counts, context_sample_counts = defaultdict(int), defaultdict(int)
        if compress == False:
            calls_output = open(file_name, 'a')
            data_line = csv.writer(calls_output, delimiter='\t', lineterminator='\n')
        else:
            logs.info("Compressing output modification calls file.")
            data_line = gzip.open(file_name + '.gz', 'wt', compresslevel=9, encoding='utf8', newline='\n')
        if per_chromosome == None:
            for i in range(0, len(keys)):
                time_m = datetime.now()
                logs.info("Starting modification calling on {} chromosome (sequence). {} seconds".format(keys[i], (time_m - time_b).total_seconds()))
                modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys[i], user_defined_context, context_total_counts, None, None)
                pileups = inbam.pileup(keys[i], ignore_overlaps=skip_clip_overlap, min_base_quality=minimum_base_quality, stepper='samtools',
                                       max_depth=max_depth, redo_baq=redo_baq, ignore_orphans=ignore_orphans, compute_baq=compute_baq,
                                       min_mapping_quality=minimum_mapping_quality, adjust_acapq_threshold=adjust_acapq_threshold)
                clean_pileup(pileups, i, modification_information_per_position, mean_mod, mean_unmod, user_defined_context, file_name, method,
                             mark_matches, mark_ends, add_indels, context_sample_counts, ignore_orphans, single_end, compress, data_line)
        else:
            time_m = datetime.now()
            logs.info("Starting modification calling on {} chromosome (sequence). {} seconds".format(keys, (time_m - time_b).total_seconds()))
            modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys, user_defined_context, context_total_counts, None, None)
            pileups = inbam.pileup(keys, ignore_overlaps=skip_clip_overlap, min_base_quality=minimum_base_quality, stepper='samtools',
                                   max_depth=max_depth, redo_baq=redo_baq, ignore_orphans=ignore_orphans, compute_baq=compute_baq,
                                   min_mapping_quality=minimum_mapping_quality, adjust_acapq_threshold=adjust_acapq_threshold)
            clean_pileup(pileups, cycles, modification_information_per_position, mean_mod, mean_unmod, user_defined_context, file_name, method,
                         mark_matches, mark_ends, add_indels, context_sample_counts, ignore_orphans, single_end, compress, data_line)
        if zero_coverage:
            for position in modification_information_per_position.keys():
                if modification_information_per_position[position][3] == 'C':
                    all_data = list((position[0], position[1], position[1] + 1, 'NA', 0, 0, 'C', 'T',
                    modification_information_per_position[position][0], modification_information_per_position[position][1], 'No'))
                    modification_calls_writer(all_data, compress, data_line, header=False)
                elif modification_information_per_position[position][3] == 'G':
                    all_data = list((position[0], position[1], position[1] + 1, 'NA', 0, 0, 'G', 'A',
                    modification_information_per_position[position][0], modification_information_per_position[position][1], 'No'))
                    modification_calls_writer(all_data, compress, data_line, header=False)
        inbam.close()
        if per_chromosome == None:
            file_name = path.join(directory, name + "_" + method + "_" + context + ".stats")
        else:
            file_name = path.join(directory, name + "_" + method + "_" + per_chromosome + "_" + context + ".stats")
        final_statistics_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context)
        if compress == False:
             calls_output.close()
        else:
            data_line.close()
        time_e = datetime.now()
        logs.info("asTair modification finder finished running. {} seconds".format((time_e - time_b).total_seconds()))
    else:
        logs.error('Mods file with this name exists. Please rename before rerunning.')
        sys.exit(1)

if __name__ == '__main__':
    call()


