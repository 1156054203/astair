#!/usr/bin/env python3


import itertools
import csv
import re
import click
import warnings
import logging
import pdb
from os import path
from datetime import datetime
from collections import defaultdict


from safe_division import non_zero_division
from bam_file_parser import bam_file_opener
from simple_fasta_parser import fasta_splitting_by_sequence
from reference_context_search_triad import sequence_context_set_creation
from reference_context_search_triad import context_sequence_search

@click.command()
@click.option('input_file', '--input_file', required=True, help='BAM format file containing sequencing reads.')
# @click.option('control_file', '--control_file', required=False, help='BAM format file containing sequencing reads of a matched control.')
@click.option('fasta_file', '--fasta_file', required=True, help='DNA sequence in fasta format used for aligning the sequencing reads and mpileup.')
@click.option('zero_coverage', '--zero_coverage', default=False, is_flag=True, help='When set to True, outputs positions not covered in the bam file. Uncovering zero coverage positions takes longer time than using the default option.')
@click.option('context', '--context', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis.')
@click.option('user_defined_context', '--user_defined_context', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
# @click.option('library', '--library', required=False, default = 'directional',  type=click.Choice(['directional', 'RR', 'PMAT']), help='Provides information for the library preparation protocol (RR is reduced representation and PMAT is post-method adapter-tagging).')
# @click.option('method', '--method', required=False, default = 'CmtoT', type=click.Choice(['CtoT', 'CmtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and CmtoT (modified cytosines are converted to thymines, TAPS-like).')
@click.option('skip_clip_overlap', '--skip_clip_overlap', required=False, is_flag=True, help='Skipping the random removal of overlapping bases between paired-end reads. Not recommended for paired-end libraries, unless the overlaps are removed prior to calling.')
@click.option('minimum_base_quality', '--minimum_base_quality', required=False, type=int, default=20, help='Set the minimum base quality for a read base to be used in the pileup (Default 20).')
# @click.option('N_threads', '--N_threads', default = 1, required=True, help='The number of threads to spawn (the default value is 1).')
@click.option('directory', '--directory', required=True, type=str, help='Output directory to save files.')
def modification_finder_exec(input_file, fasta_file, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, directory):
        cytosine_modification_finder(input_file, fasta_file, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, directory)


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def modification_calls_writer(name, directory, data_mods, header=False):
    """Outputs the modification calls per position in a tab-delimited format."""
    with open(path.join(directory, name+".mods"), 'a', newline='') as calls_output:
        data_line = csv.writer(calls_output, delimiter='\t', lineterminator='\n')
        if header:
            data_line.writerow(["CHROM", "START", "END", "MOD_LEVEL", "MOD", "unmod", "ALT", "REF", "DEPTH", "CONTEXT",
                         "SPECIFIC_CONTEXT", 'SNV'])
            data_line.writerow(data_mods)
        else:
            data_line.writerow(data_mods)

def statistics_calculator(mean_mod, mean_unmod, data_mod, user_defined_context):
    """Calculates the summary statistics of the cytosine modificaton levels."""
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


def final_statistics_output(mean_mod, mean_unmod, directory, name, user_defined_context):
    """Writes the summary statistics of the cytosine modificaton levels."""
    with open(path.join(directory, name+".stats"), 'a', newline='') as statistics_output:
        wr = csv.writer(statistics_output, delimiter='\t', lineterminator='\n')
        wr.writerow([
                        "_____________________________________________________________________________________________________________________________"])
        wr.writerow([
                        "Cytosine modification rate given as the percentage total modified cytosines divided by the total number of cytosines covered."])
        wr.writerow([
                        "_____________________________________________________________________________________________________________________________"])
        wr.writerow(["mean CpG modification rate: ",
                     str(round(non_zero_division(mean_mod['CpG'], mean_mod['CpG'] + mean_unmod['CpG']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CGA modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CGA'], mean_mod['CGA'] + mean_unmod['CGA']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CGC modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CGC'], mean_mod['CGC'] + mean_unmod['CGC']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CGG modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CGG'], mean_mod['CGG'] + mean_unmod['CGG']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CGT modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CGT'], mean_mod['CGT'] + mean_unmod['CGT']) * 100, 3)) + ' %'])
        wr.writerow(["mean CHG modification rate: ",
                     str(round(non_zero_division(mean_mod['CHG'], mean_mod['CHG'] + mean_unmod['CHG']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CAG modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CAG'], mean_mod['CAG'] + mean_unmod['CAG']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CCG modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CCG'], mean_mod['CCG'] + mean_unmod['CCG']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CTG modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CTG'], mean_mod['CTG'] + mean_unmod['CTG']) * 100, 3)) + ' %'])
        wr.writerow(["mean CHH modification rate: ",
                     str(round(non_zero_division(mean_mod['CHH'], mean_mod['CHH'] + mean_unmod['CHH']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CTT modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CTT'], mean_mod['CTT'] + mean_unmod['CTT']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CAT modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CAT'], mean_mod['CAT'] + mean_unmod['CAT']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CCT modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CCT'], mean_mod['CCT'] + mean_unmod['CCT']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CTA modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CTA'], mean_mod['CTA'] + mean_unmod['CTA']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CAA modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CAA'], mean_mod['CAA'] + mean_unmod['CAA']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CCA modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CCA'], mean_mod['CCA'] + mean_unmod['CCA']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CTC modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CTC'], mean_mod['CTC'] + mean_unmod['CTC']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CAC modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CAC'], mean_mod['CAC'] + mean_unmod['CAC']) * 100, 3)) + ' %'])
        wr.writerow([" | " + "mean CCC modification rate:" + " | " + str(
            round(non_zero_division(mean_mod['CCC'], mean_mod['CCC'] + mean_unmod['CCC']) * 100, 3)) + ' %'])
        wr.writerow(["mean CNN modification rate: ", str(
            round(non_zero_division(mean_mod['Unknown'], mean_mod['Unknown'] + mean_unmod['Unknown']) * 100, 3)) + ' %'])
        if user_defined_context:
            wr.writerow(["mean " + user_defined_context + " modification rate: ", str(
                round(non_zero_division(mean_mod['user defined context'], mean_mod['user defined context'] + mean_unmod['user defined context']) * 100, 3)) + ' %'])
        wr.writerow([
                        "_____________________________________________________________________________________________________________________________"])
        
        
        
def pillup_summary(modification_information_per_position, position, read_counts, mean_mod, mean_unmod, name, directory, user_defined_context, header):
    """Gives the modication call rows given strand information."""
    if modification_information_per_position[position][3] == 'C':
        desired_tuples = [(147, 'C'), (99, 'C'), (147, 'T'), (99, 'T')]
        undesired_tuples = [(163, 'G'), (83, 'G'), (163, 'A'), (83, 'A')]
        modification = 'T'
        reference = 'C'
    elif modification_information_per_position[position][3] == 'G':
        desired_tuples = [(163, 'G'), (83, 'G'), (163, 'A'), (83, 'A')]
        undesired_tuples = [(147, 'C'), (99, 'C'), (147, 'T'), (99, 'T')]
        modification = 'A'
        reference = 'G'
    if non_zero_division(read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]],
                         (read_counts[undesired_tuples[0]] + read_counts[undesired_tuples[1]]
                              + read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]])) < 0.8:
        snp = 'No'
    elif non_zero_division(read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]],
                             (read_counts[undesired_tuples[0]] + read_counts[undesired_tuples[1]]
                                  + read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]])) >= 0.8:
        snp = 'homozyguous'
    all_data = list((position[0], position[1], position[1] + 1, round(
        non_zero_division(read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]], (
            read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]] + read_counts[desired_tuples[0]] + read_counts[
                desired_tuples[1]])), 3), read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]],
                     read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]], modification, reference,
                     modification_information_per_position[position][0],
                     modification_information_per_position[position][1], snp))
    statistics_calculator(mean_mod, mean_unmod, all_data, user_defined_context)
    modification_calls_writer(name, directory, all_data, header=header)
    

def cytosine_modification_finder(input_file, fasta_file, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, directory):
    """Searches for cytosine modification positions in the desired contexts and calculates the modificaton levels."""
    #header = True
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if user_defined_context:
        mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'user defined context': 0,  'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
        mean_unmod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'user defined context': 0, 'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
    else:
        mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
        mean_unmod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
    inbam = bam_file_opener(input_file, None)
    keys, fastas = fasta_splitting_by_sequence(fasta_file)
    contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
    for i in range(0, len(keys)):
        time_m = datetime.now()
        logs.info("Starting modification calling on {} chromosome (sequence). {} seconds".format(keys[i], (time_m - time_b).total_seconds()))
        modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys[i], user_defined_context)
        covered_positions = list()
        if skip_clip_overlap:
            clip_overlap = False
        else:
            clip_overlap = True
        pileups = inbam.pileup(keys[i], ignore_overlaps=clip_overlap, min_base_quality=minimum_base_quality, fasta_file=fasta_file, stepper='samtools', max_depth=100000)
        cycles=0
        for reads in pileups:
            if cycles==0:
                header = True
            else:
                header = False
            if (reads.reference_name, reads.pos, reads.pos + 1) in modification_information_per_position:
                position = (reads.reference_name, reads.pos, reads.pos + 1)
                covered_positions.append(position)
                read_counts = defaultdict(int)
                try:
                    sequences = reads.get_query_sequences()
                except AssertionError:
                    logs.exception("Failed getting query sequences (AssertionError, pysam)")
                    continue
                for pileup, seq in itertools.zip_longest(reads.pileups, sequences, fillvalue='BLANK'):
                    read_counts[(pileup.alignment.flag, seq)] += 1
                pillup_summary(modification_information_per_position, position, read_counts, mean_mod, mean_unmod, name, directory, user_defined_context, header)
                cycles+=1
    if zero_coverage:
        positions_not_covered = [pos for pos in modification_information_per_position if pos not in covered_positions]
        positions_not_covered.sort()
        for positions in positions_not_covered:
            if modification_information_per_position[positions][3] == 'C':
                all_data = list((positions[0], positions[1], positions[1] + 1, 0, 0, 0, 'T', 'C', modification_information_per_position[positions][0], modification_information_per_position[positions][1], 'No'))
                modification_calls_writer(name, directory, all_data, header=header)
                header = False
            elif modification_information_per_position[positions][3] == 'G':
                all_data = list((positions[0], positions[1], positions[1] + 1, 0, 0, 0, 'A', 'G', modification_information_per_position[positions][0], modification_information_per_position[positions][1], 'No'))
                modification_calls_writer(name, directory, all_data, header=header)
                header = False
    final_statistics_output(mean_mod, mean_unmod, directory, name, user_defined_context)
    time_e = datetime.now()
    logs.info("asTair modification finder finished running. {} seconds".format((time_e - time_b).total_seconds()))


if __name__ == '__main__':
    modification_finder_exec()


