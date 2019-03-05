#!/usr/bin/env python

from __future__ import print_function

import re
import os
import sys
import pdb
import csv
import click
import logging
import warnings
from os import path
from datetime import datetime
from collections import defaultdict

if sys.version[0] == '3':
    from itertools import zip_longest
elif sys.version[0] == '2':
    from itertools import izip_longest as zip_longest
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))

from safe_division import non_zero_division
from bam_file_parser import bam_file_opener
from simple_fasta_parser import fasta_splitting_by_sequence
from reference_context_search_triad import sequence_context_set_creation
from reference_context_search_triad import context_sequence_search


@click.command()
@click.option('input_file', '--input_file', '-i', required=True, help='BAM format file containing sequencing reads.')
# @click.option('control_file', '--control_file', '-c', required=False, help='BAM format file containing sequencing reads of a matched control.')
@click.option('fasta_file', '--fasta_file', '-f', required=True, help='DNA sequence in fasta format used for aligning the sequencing reads and mpileup.')
@click.option('zero_coverage', '--zero_coverage', '-z', default=False, is_flag=True, help='When set to True, outputs positions not covered in the bam file. Uncovering zero coverage positions takes longer time than using the default option.')
@click.option('context', '--context', '-co', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis.')
@click.option('user_defined_context', '--user_defined_context', '-uc', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
# @click.option('library', '--library', '-l', required=False, default = 'directional',  type=click.Choice(['directional', 'rr', 'ptat']), help='Provides information for the library preparation protocol (rr is reduced representation and ptat is post-method adapter-tagging).')
@click.option('method', '--method', '-m', required=False, default = 'mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like).')
@click.option('skip_clip_overlap', '--skip_clip_overlap', '-sc', required=False, default=False, type=bool, help='Skipping the random removal of overlapping bases between paired-end reads. Not recommended for paired-end libraries, unless the overlaps are removed prior to calling. (Default False)')
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
@click.option('per_chromosome', '--per_chromosome', '-chr', default=None, type=str, help='When used, it calculates the modification rates only per the chromosome given. (Default None')
@click.option('N_threads', '--N_threads', '-t', default=10, required=True, help='The number of threads to spawn (the default value is 10).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
def modification_finder_exec(input_file, fasta_file, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, method, minimum_mapping_quality, adjust_acapq_threshold,
                             mark_matches, mark_ends, add_indels, redo_baq, compute_baq, ignore_orphans, max_depth,per_chromosome, N_threads, directory):
        cytosine_modification_finder(input_file, fasta_file, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, method, minimum_mapping_quality, adjust_acapq_threshold,
                                     mark_matches, mark_ends, add_indels, redo_baq, compute_baq, ignore_orphans, max_depth, per_chromosome, N_threads, directory)


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def modification_calls_writer(data_mods, file_name, header=False):
    """Outputs the modification calls per position in a tab-delimited format."""
    try:
        with open(file_name, 'a') as calls_output:
            data_line = csv.writer(calls_output, delimiter='\t', lineterminator='\n')
            if header:
                data_line.writerow(["CHROM", "START", "END", "MOD_LEVEL", "MOD", "UNMOD", "REF", "ALT", "CONTEXT",
                             "SPECIFIC_CONTEXT", 'SNV', 'TOTAL_DEPTH'])
                data_line.writerow(data_mods)
            else:
                data_line.writerow(data_mods)
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
        mean_mod['Unknown'] += data_mod[4]
        mean_unmod['Unknown'] += data_mod[5]
    if user_defined_context and re.match(user_defined_context, data_mod[9]) and data_mod[10] == 'No':
        mean_mod['user defined context'] += data_mod[4]
        mean_unmod['user defined context'] += data_mod[5]


def final_statistics_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts):
    """Writes the summary statistics of the cytosine modificaton levels.
    Cytosine modification rate given as the percentage total modified cytosines
    divided by the total number of cytosines covered."""
    try:
        with open(file_name, 'w') as statistics_output:
            wr = csv.writer(statistics_output, delimiter='\t', lineterminator='\n')
            wr.writerow(["CONTEXT", "SPECIFIC_CONTEXT", "MEAN_MODIFICATION_RATE_PERCENT", "TOTAL_POSITIONS", "COVERED_POSITIONS"])
            wr.writerow(["CpG", "", round(non_zero_division(mean_mod['CpG'], mean_mod['CpG'] + mean_unmod['CpG']) * 100, 3),
            context_total_counts['CG']+context_total_counts['CGb'], context_sample_counts['CpG']])
            wr.writerow(["", "CGA", round(non_zero_division(mean_mod['CGA'], mean_mod['CGA'] + mean_unmod['CGA']) * 100, 3),
                         context_total_counts['CGA'], context_sample_counts['CGA']])
            wr.writerow(["", "CGC", round(non_zero_division(mean_mod['CGC'], mean_mod['CGC'] + mean_unmod['CGC']) * 100, 3),
                         context_total_counts['CGC'], context_sample_counts['CGC']])
            wr.writerow(["", "CGG", round(non_zero_division(mean_mod['CGG'], mean_mod['CGG'] + mean_unmod['CGG']) * 100, 3),
                         context_total_counts['CGG'], context_sample_counts['CGG']])
            wr.writerow(["", "CGT", round(non_zero_division(mean_mod['CGT'], mean_mod['CGT'] + mean_unmod['CGT']) * 100, 3),
                         context_total_counts['CGT'], context_sample_counts['CGT']])
            wr.writerow(["CHG", "", round(non_zero_division(mean_mod['CHG'], mean_mod['CHG'] + mean_unmod['CHG']) * 100, 3),
                         context_total_counts['CHG']+context_total_counts['CHGb'], context_sample_counts['CHG']])
            wr.writerow(["", "CAG", round(non_zero_division(mean_mod['CAG'], mean_mod['CAG'] + mean_unmod['CAG']) * 100, 3),
                         context_total_counts['CAG'], context_sample_counts['CAG']])
            wr.writerow(["", "CCG", round(non_zero_division(mean_mod['CCG'], mean_mod['CCG'] + mean_unmod['CCG']) * 100, 3),
                         context_total_counts['CCG'], context_sample_counts['CCG']])
            wr.writerow(["", "CTG", round(non_zero_division(mean_mod['CTG'], mean_mod['CTG'] + mean_unmod['CTG']) * 100, 3),
                         context_total_counts['CTG'], context_sample_counts['CTG']])
            wr.writerow(["CHH", "", round(non_zero_division(mean_mod['CHH'], mean_mod['CHH'] + mean_unmod['CHH']) * 100, 3),
                         context_total_counts['CHH']+context_total_counts['CHHb'], context_sample_counts['CHH']])
            wr.writerow(["", "CTT", round(non_zero_division(mean_mod['CTT'], mean_mod['CTT'] + mean_unmod['CTT']) * 100, 3),
                         context_total_counts['CTT'], context_sample_counts['CTT']])
            wr.writerow(["", "CAT", round(non_zero_division(mean_mod['CAT'], mean_mod['CAT'] + mean_unmod['CAT']) * 100, 3),
                         context_total_counts['CAT'], context_sample_counts['CAT']])
            wr.writerow(["", "CCT", round(non_zero_division(mean_mod['CCT'], mean_mod['CCT'] + mean_unmod['CCT']) * 100, 3),
                         context_total_counts['CCT'], context_sample_counts['CCT']])
            wr.writerow(["", "CTA", round(non_zero_division(mean_mod['CTA'], mean_mod['CTA'] + mean_unmod['CTA']) * 100, 3),
                         context_total_counts['CTA'], context_sample_counts['CTA']])
            wr.writerow(["", "CAA", round(non_zero_division(mean_mod['CAA'], mean_mod['CAA'] + mean_unmod['CAA']) * 100, 3),
                         context_total_counts['CAA'], context_sample_counts['CAA']])
            wr.writerow(["", "CCA", round(non_zero_division(mean_mod['CCA'], mean_mod['CCA'] + mean_unmod['CCA']) * 100, 3),
                         context_total_counts['CCA'], context_sample_counts['CCA']])
            wr.writerow(["", "CTC", round(non_zero_division(mean_mod['CTC'], mean_mod['CTC'] + mean_unmod['CTC']) * 100, 3),
                         context_total_counts['CTC'], context_sample_counts['CTC']])
            wr.writerow(["", "CAC", round(non_zero_division(mean_mod['CAC'], mean_mod['CAC'] + mean_unmod['CAC']) * 100, 3),
                         context_total_counts['CAC'], context_sample_counts['CAC']])
            wr.writerow(["", "CCC", round(non_zero_division(mean_mod['CCC'], mean_mod['CCC'] + mean_unmod['CCC']) * 100, 3),
                         context_total_counts['CCC'], context_sample_counts['CCC']])

            wr.writerow(["CNN", "", round(non_zero_division(mean_mod['Unknown'], mean_mod['Unknown'] + mean_unmod['Unknown']) * 100, 3),
                         context_total_counts['CN']+context_total_counts['CNb'], context_sample_counts['CN']])
            if user_defined_context:
                wr.writerow([user_defined_context, "",round(non_zero_division(mean_mod['user defined context'], mean_mod['user defined context'] +
                        mean_unmod['user defined context']) * 100, 3), context_total_counts['user defined context'], context_sample_counts['user defined context']])
    except IOError:
        logs.error('asTair cannot write to statistics summary file.', exc_info=True)


def flags_expectation(modification_information_per_position, position):
    """Gives the expected flag-base couples, the reference and the modified base."""
    if modification_information_per_position[position][3] == 'C':
        desired_tuples = [(147, 'C'), (99, 'C'), (147, 'T'), (99, 'T')]
        undesired_tuples = [(163, 'C'), (83, 'C'), (163, 'T'), (83, 'T')]
        modification = 'T'
        reference = 'C'
    elif modification_information_per_position[position][3] == 'G':
        desired_tuples = [(163, 'G'), (83, 'G'), (163, 'A'), (83, 'A')]
        undesired_tuples = [(147, 'G'), (99, 'G'), (147, 'A'), (99, 'A')]
        modification = 'A'
        reference = 'G'
    return desired_tuples, undesired_tuples, modification, reference

        
def pillup_summary(modification_information_per_position, position, read_counts, mean_mod, mean_unmod, user_defined_context,
                   header, file_name, desired_tuples, undesired_tuples, modification, reference, depth, method, context_sample_counts):
    """Gives the modication call rows given strand information."""
    if non_zero_division(read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]],
                         (read_counts[undesired_tuples[0]] + read_counts[undesired_tuples[1]]
                              + read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]])) < 0.8:
        snp = 'No'
    elif non_zero_division(read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]],
                             (read_counts[undesired_tuples[0]] + read_counts[undesired_tuples[1]]
                                  + read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]])) >= 0.8:
        snp = 'homozyguous'
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
    statistics_calculator(mean_mod, mean_unmod, all_data, user_defined_context, context_sample_counts)
    modification_calls_writer(all_data, file_name, header=header)


def clean_pileup(pileups, cycles, modification_information_per_position, mean_mod, mean_unmod, user_defined_context, file_name, method,
                 mark_matches, mark_ends, add_indels, context_sample_counts):
    """Takes reads from the piled-up region and calculates modification levels."""
    for reads in pileups:
        if cycles == 0:
            header = True
        else:
            header = False
        if (reads.reference_name, reads.pos, reads.pos + 1) in modification_information_per_position:
            position = (reads.reference_name, reads.pos, reads.pos + 1)
            desired_tuples, undesired_tuples, modification, reference = flags_expectation(modification_information_per_position, position)
            read_counts = defaultdict(int)
            try:
                sequences = reads.get_query_sequences(mark_matches=mark_matches, mark_ends=mark_ends, add_indels=add_indels)
            except AssertionError:
                logs.exception("Failed getting query sequences (AssertionError, pysam). Please decrease the max_depth parameter.")
                continue
            for pileup, seq in zip_longest(reads.pileups, sequences, fillvalue='BLANK'):
                if pileup != 'BLANK':
                    read_counts[(pileup.alignment.flag, seq.upper())] += 1
            pillup_summary(modification_information_per_position, position, read_counts, mean_mod, mean_unmod, user_defined_context, header,
                           file_name, desired_tuples, undesired_tuples, modification, reference, reads.get_num_aligned(), method, context_sample_counts)
            cycles+=1
            modification_information_per_position.pop(position)

def cytosine_modification_finder(input_file, fasta_file, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, method,
                                 minimum_mapping_quality, adjust_acapq_threshold, mark_matches, mark_ends, add_indels, redo_baq, compute_baq, ignore_orphans,
                                 max_depth, per_chromosome, N_threads, directory):
    """Searches for cytosine modification positions in the desired contexts and calculates the modificaton levels."""
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if per_chromosome == None:
        file_name = path.join(directory, name + "_" + method + "_" + context + ".mods")
    else:
        file_name = path.join(directory, name + "_" + method + "_" + per_chromosome + "_" + context + ".mods")
    if not os.path.isfile(file_name):
        if user_defined_context:
            mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0,
                        'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'user defined context': 0,
                        'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
            mean_unmod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0,
                          'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'user defined context': 0,
                          'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
        else:
            mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0,
                        'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
            mean_unmod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0,
                          'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
        try:
            inbam = bam_file_opener(input_file, None, N_threads)
        except Exception:
            sys.exit(1)
        try:
            keys, fastas = fasta_splitting_by_sequence(fasta_file, per_chromosome)
        except Exception:
            sys.exit(1)
        contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
        cycles = 0
        context_total_counts, context_sample_counts = defaultdict(int), defaultdict(int)
        if per_chromosome == None:
            for i in range(0, len(keys)):
                time_m = datetime.now()
                logs.info("Starting modification calling on {} chromosome (sequence). {} seconds".format(keys[i], (time_m - time_b).total_seconds()))
                modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys[i], user_defined_context, context_total_counts, None)
                pileups = inbam.pileup(keys[i], ignore_overlaps=skip_clip_overlap, min_base_quality=minimum_base_quality, fasta_file=fasta_file, stepper='samtools',
                                       max_depth=max_depth, redo_baq=redo_baq, ignore_orphans=ignore_orphans, compute_baq=compute_baq,
                                       min_mapping_quality=minimum_mapping_quality, adjust_acapq_threshold=adjust_acapq_threshold)
                clean_pileup(pileups, cycles, modification_information_per_position, mean_mod, mean_unmod, user_defined_context, file_name, method,
                             mark_matches, mark_ends, add_indels, context_sample_counts)
        else:
            time_m = datetime.now()
            logs.info("Starting modification calling on {} chromosome (sequence). {} seconds".format(keys, (time_m - time_b).total_seconds()))
            modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys, user_defined_context, context_total_counts, None)
            pileups = inbam.pileup(keys, ignore_overlaps=skip_clip_overlap, min_base_quality=minimum_base_quality, fasta_file=fasta_file, stepper='samtools',
                                   max_depth=max_depth, redo_baq=redo_baq, ignore_orphans=ignore_orphans, compute_baq=compute_baq,
                                   min_mapping_quality=minimum_mapping_quality, adjust_acapq_threshold=adjust_acapq_threshold)
            clean_pileup(pileups, cycles, modification_information_per_position, mean_mod, mean_unmod, user_defined_context, file_name, method,
                         mark_matches, mark_ends, add_indels, context_sample_counts)
        if zero_coverage:
            for position in modification_information_per_position.keys():
                if modification_information_per_position[position][3] == 'C':
                    all_data = list((position[0], position[1], position[1] + 1, 'NA', 0, 0, 'C', 'T',
                    modification_information_per_position[position][0], modification_information_per_position[position][1], 'No'))
                    modification_calls_writer(all_data, file_name, header=False)
                elif modification_information_per_position[position][3] == 'G':
                    all_data = list((position[0], position[1], position[1] + 1, 'NA', 0, 0, 'G', 'A',
                    modification_information_per_position[position][0], modification_information_per_position[position][1], 'No'))
                    modification_calls_writer(all_data, file_name, header=False)
        if per_chromosome == None:
            file_name = path.join(directory, name + "_" + method + "_" + context + ".stats")
        else:
            file_name = path.join(directory, name + "_" + method + "_" + per_chromosome + "_" + context + ".stats")
        final_statistics_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts)
        time_e = datetime.now()
        logs.info("asTair modification finder finished running. {} seconds".format((time_e - time_b).total_seconds()))
    else:
        logs.error('Mods file with this name exists. Please rename before rerunning.', exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    modification_finder_exec()


