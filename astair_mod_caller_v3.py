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
@click.option('input_file', '--input_file', '-i', required=True, help='BAM format file containing sequencing reads.')
# @click.option('control_file', '--control_file', '-c', required=False, help='BAM format file containing sequencing reads of a matched control.')
@click.option('fasta_file', '--fasta_file', '-f', required=True, help='DNA sequence in fasta format used for aligning the sequencing reads and mpileup.')
@click.option('zero_coverage', '--zero_coverage', '-z', default=False, is_flag=True, help='When set to True, outputs positions not covered in the bam file. Uncovering zero coverage positions takes longer time than using the default option.')
@click.option('context', '--context', '-co', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis.')
@click.option('user_defined_context', '--user_defined_context', '-uc', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
# @click.option('library', '--library', '-l', required=False, default = 'directional',  type=click.Choice(['directional', 'RR', 'PMAT']), help='Provides information for the library preparation protocol (RR is reduced representation and PMAT is post-method adapter-tagging).')
# @click.option('method', '--method', '-m', required=False, default = 'CmtoT', type=click.Choice(['CtoT', 'CmtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and CmtoT (modified cytosines are converted to thymines, TAPS-like).')
@click.option('skip_clip_overlap', '--skip_clip_overlap', '-sc', required=False, is_flag=True, help='Skipping the random removal of overlapping bases between paired-end reads. Not recommended for paired-end libraries, unless the overlaps are removed prior to calling.')
@click.option('minimum_base_quality', '--minimum_base_quality', '-bq', required=False, type=int, default=20, help='Set the minimum base quality for a read base to be used in the pileup (Default 20).')
@click.option('per_chromosome', '--per_chromosome', '-chr', default=None, type=str, help='When used, it calculates the modification rates only per the chromosome given. (Default None')
# @click.option('N_threads', '--N_threads', '-t', default = 1, required=True, help='The number of threads to spawn (the default value is 1).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
def modification_finder_exec(input_file, fasta_file, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, per_chromosome, directory):
        cytosine_modification_finder(input_file, fasta_file, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, per_chromosome, directory)


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def modification_calls_writer(data_mods, file_name, header=False):
    """Outputs the modification calls per position in a tab-delimited format."""
    with open(file_name, 'a', newline='') as calls_output:
        data_line = csv.writer(calls_output, delimiter='\t', lineterminator='\n')
        if header:
            data_line.writerow(["CHROM", "START", "END", "MOD_LEVEL", "MOD", "UNMOD", "ALT", "REF",  "CONTEXT",
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


def final_statistics_output(mean_mod, mean_unmod, user_defined_context, file_name):
    """Writes the summary statistics of the cytosine modificaton levels."""
    with open(file_name, 'a', newline='') as statistics_output:
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
        
        
        
def pillup_summary(modification_information_per_position, position, read_counts, mean_mod, mean_unmod, user_defined_context, header, file_name):
    """Gives the modication call rows given strand information."""
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
    if non_zero_division(read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]],
                         (read_counts[undesired_tuples[0]] + read_counts[undesired_tuples[1]]
                              + read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]])) < 0.8:
        snp = 'No'
    elif non_zero_division(read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]],
                             (read_counts[undesired_tuples[0]] + read_counts[undesired_tuples[1]]
                                  + read_counts[undesired_tuples[2]] + read_counts[undesired_tuples[3]])) >= 0.8:
        snp = 'homozyguous'
    all_data = list((position[0], position[1], position[1] + 1, round(
        non_zero_division(read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]], (
            read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]] + read_counts[desired_tuples[0]] + read_counts[
                desired_tuples[1]])), 3), read_counts[desired_tuples[2]] + read_counts[desired_tuples[3]],
                     read_counts[desired_tuples[0]] + read_counts[desired_tuples[1]], modification, reference,
                     modification_information_per_position[position][0],
                     modification_information_per_position[position][1], snp))
    statistics_calculator(mean_mod, mean_unmod, all_data, user_defined_context)
    modification_calls_writer(all_data, file_name, header=header)


def clean_pileup(pileups, cycles, modification_information_per_position, mean_mod, mean_unmod, user_defined_context, file_name):
    """Takes reads from the piled-up region and calculates modification levels."""
    for reads in pileups:
        if cycles == 0:
            header = True
        else:
            header = False
        if (reads.reference_name, reads.pos, reads.pos + 1) in modification_information_per_position:
            position = (reads.reference_name, reads.pos, reads.pos + 1)
            read_counts = defaultdict(int)
            try:
                sequences = reads.get_query_sequences()
            except AssertionError:
                logs.exception("Failed getting query sequences (AssertionError, pysam)")
                continue
            for pileup, seq in itertools.zip_longest(reads.pileups, sequences, fillvalue='BLANK'):
                if pileup.indel == 0 and pileup.is_del == 0 and pileup.is_refskip == 0:
                    read_counts[(pileup.alignment.flag, seq.upper())] += 1
            pillup_summary(modification_information_per_position, position, read_counts, mean_mod, mean_unmod, user_defined_context, header, file_name)
            cycles+=1
            modification_information_per_position.pop(position)

def cytosine_modification_finder(input_file, fasta_file, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, per_chromosome, directory):
    """Searches for cytosine modification positions in the desired contexts and calculates the modificaton levels."""
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if per_chromosome == None:
        file_name = path.join(directory, name + "_" + context + ".mods")
    else:
        file_name = path.join(directory, name + "_" + per_chromosome + "_" + context + ".mods")
    if user_defined_context:
        mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'user defined context': 0,  'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
        mean_unmod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'user defined context': 0, 'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
    else:
        mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
        mean_unmod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'CGA':0, 'CGT':0, 'CGC':0, 'CGG':0}
    inbam = bam_file_opener(input_file, None)
    keys, fastas = fasta_splitting_by_sequence(fasta_file, per_chromosome)
    contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
    if skip_clip_overlap:
        clip_overlap = False
    else:
        clip_overlap = True
    cycles = 0
    if per_chromosome == None:
        for i in range(0, len(keys)):
            time_m = datetime.now()
            logs.info("Starting modification calling on {} chromosome (sequence). {} seconds".format(keys[i], (time_m - time_b).total_seconds()))
            modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys[i], user_defined_context)
            pileups = inbam.pileup(keys[i], ignore_overlaps=clip_overlap, min_base_quality=minimum_base_quality, fasta_file=fasta_file, stepper='samtools', max_depth=250)
            clean_pileup(pileups, cycles, modification_information_per_position, mean_mod, mean_unmod, user_defined_context, file_name)
    else:
        time_m = datetime.now()
        logs.info("Starting modification calling on {} chromosome (sequence). {} seconds".format(keys, (time_m - time_b).total_seconds()))
        modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys, user_defined_context)
        pileups = inbam.pileup(keys, ignore_overlaps=clip_overlap, min_base_quality=minimum_base_quality, fasta_file=fasta_file, stepper='samtools', max_depth=250)
        clean_pileup(pileups, cycles, modification_information_per_position, mean_mod, mean_unmod, user_defined_context, file_name)
    if zero_coverage:
        for position in modification_information_per_position.keys():
            if modification_information_per_position[position][3] == 'C':
                all_data = list((position[0], position[1], position[1] + 1, 'NA', 0, 0, 'T', 'C',
                modification_information_per_position[position][0], modification_information_per_position[position][1], 'No'))
                modification_calls_writer(all_data, file_name, header=False)
            elif modification_information_per_position[position][3] == 'G':
                all_data = list((position[0], position[1], position[1] + 1, 'NA', 0, 0, 'A', 'G',
                modification_information_per_position[position][0], modification_information_per_position[position][1], 'No'))
                modification_calls_writer(all_data, file_name, header=False)
    if per_chromosome == None:
        file_name = path.join(directory, name + "_" + context + ".stats")
    else:
        file_name = path.join(directory, name + "_" + per_chromosome + "_" + context + ".stats")
    final_statistics_output(mean_mod, mean_unmod, user_defined_context, file_name)
    time_e = datetime.now()
    logs.info("asTair modification finder finished running. {} seconds".format((time_e - time_b).total_seconds()))


if __name__ == '__main__':
    modification_finder_exec()


