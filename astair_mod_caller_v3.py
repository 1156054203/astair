#!/usr/bin/env python3

import pysam
import ahocorasick
import itertools
import csv
import re
import sys
import click
import warnings
import logging
from os import path
from datetime import datetime
from collections import defaultdict


from DNA_sequences_operations import complementary
from safe_division import non_zero_division

@click.command()
@click.option('input_file', '--input_file', required=True, help='BAM format file containing sequencing reads.')
# @click.option('control_file', '--control_file', required=False, help='BAM format file containing sequencing reads of a matched control.')
@click.option('fasta_file', '--fasta_file', required=True, help='DNA sequence in fasta format used for aligning the sequencing reads and mpileup.')
@click.option('zero_coverage', '--zero_coverage', default=False, is_flag=True, help='When set to True, outputs positions not covered in the bam file. Uncovering zero coverage positions takes longer time than using the default option.')
@click.option('context', '--context', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis.')
@click.option('user_defined_context', '--user_defined_context', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
# @click.option('library', '--library', required=False, default = 'directional',  type=click.Choice(['directional', 'RR', 'PMAT']), help='Provides information for the library preparation protocol (RR is reduced representation and PMAT is post-method adapter-tagging).')
# @click.option('method', '--method', required=False, default = 'CmtoT', type=click.Choice(['CtoT', 'CmtoT']), help='Specify sequencing method, possible options are CtoT (uunmodified cytosines are converted to thymines, bisulfite sequencing-like) and CmtoT (modified cytosines are converted to thymines, TAPS-like).')
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


def fasta_splitting_by_sequence(fasta_file):
    """Parses fasta files with multiple genomes."""
    fastas = {}
    keys, sequences = list(), list()
    with open(fasta_file, 'r', newline='') as fasta_handle:
        fasta_sequence = fasta_handle.read()
        if re.match(r".*(?=\r\n)", fasta_sequence):
            keys = re.findall(r"(?<=>).*(?=\r\n)", fasta_sequence)
            sequences = re.findall(r"(?<=\r\n)(?!>).*(?=\r\n)", fasta_sequence)
            if len(keys) < len(sequences):
                new_strings = [string for string in sequences if len(string) <= 75]
                sequences = [string for string in sequences if string not in new_strings]
                joined_sequence = "".join(new_strings)
                sequences.insert(0, joined_sequence)
        elif re.match(r".*(?=\n)", fasta_sequence):
            keys = re.findall(r"(?<=>).*(?=\n)", fasta_sequence)
            sequences = re.findall(r"(?<=\n)(?!>).*(?=\n)", fasta_sequence)
    for i in range(0, len(keys)):
        fastas[keys[i]] = sequences[i]
    return keys, fastas


def sequence_context_set_creation(desired_sequence, user_defined_context):
    """Prepares sets of possible cytosine contexts."""
    letters_top = ['A', 'C', 'T', 'a', 'c', 't']
    if user_defined_context:
        user = list(map(''.join, itertools.product(*zip(user_defined_context.upper(), user_defined_context.lower()))))
    if desired_sequence == 'all':
        CHG = [y + x + z for x in letters_top for y in ['C', 'c'] for z in ['G', 'g']]
        CHGb = [y + x + z for x in letters_top for z in ['C', 'c'] for y in ['G', 'g']]
        CHH = [y + x + z for x in letters_top for y in ['C', 'c'] for z in letters_top]
        CHHb = [y + x + z for x in letters_top for z in ['C', 'c'] for y in letters_top]
        CG = [y + z for y in ['C', 'c'] for z in ['G', 'g']]
        CN = [y + x + z for x in ['N', 'n'] for y in ['C', 'c'] for z in ['N', 'n']]
        CNb = [y + x + z for x in ['N', 'n'] for z in ['C', 'c'] for y in ['N', 'n']]
        if user_defined_context:
            contexts = {'CHG': list(CHG), 'CHGb': list(CHGb), 'CHH': list(CHH), 'CHHb': list(CHHb), 'CG': list(CG),'CN': list(CN), 'CNb': list(CNb), 'user': list(user)}
            all_keys = list(('CHG', 'CHGb', 'CHH', 'CHHb', 'CG', 'CN', 'CNb', 'user'))
        else:
            contexts = {'CHG': list(CHG), 'CHGb': list(CHGb), 'CHH': list(CHH), 'CHHb': list(CHHb), 'CG': list(CG),'CN': list(CN), 'CNb': list(CNb)}
            all_keys = list(('CHG', 'CHGb', 'CHH', 'CHHb', 'CG', 'CN', 'CNb'))
    elif desired_sequence == 'CpG':
        CG = [y + z for y in ['C', 'c'] for z in ['G', 'g']]
        if user_defined_context:
            contexts = {'CG': list(CG), 'user': list(user)}
            all_keys = list(( 'CG', 'user'))
        else:
            contexts = {'CG': list(CG)}
            all_keys = list(('CG'))
    elif desired_sequence == 'CHG':
        CHG = [y + x + z for x in letters_top for y in ['C', 'c'] for z in ['G', 'g']]
        CHGb = [y + x + z for x in letters_top for z in ['C', 'c'] for y in ['G', 'g']]
        if user_defined_context:
            contexts = {'CHG': list(CHG), 'CHGb': list(CHGb), 'user': list(user)}
            all_keys = list(( 'CHG', 'CHGb', 'user'))
        else:
            contexts = {'CHG': list(CHG), 'CHGb': list(CHGb)}
            all_keys = list(('CHG', 'CHGb'))
    elif desired_sequence == 'CHH':
        CHH = [y + x + z for x in letters_top for y in ['C', 'c'] for z in letters_top]
        CHHb = [y + x + z for x in letters_top for z in ['C', 'c'] for y in letters_top]
        if user_defined_context:
            contexts = {'CHH': list(CHH), 'CHHb': list(CHHb), 'user': list(user)}
            all_keys = list(( 'CHH', 'CHHb', 'user'))
        else:
            contexts = {'CHH': list(CHH), 'CHHb': list(CHHb)}
            all_keys = list(('CHH', 'CHHb'))
    return contexts, all_keys


def ahocorasick_search(objects, context, string, string_name, user_defined_context, data_context):
    """Looks for cytosine contexts in the reference fasta file."""
    auto = ahocorasick.Automaton()
    for pattern in context[objects]:
        auto.add_word(pattern, pattern)
    auto.make_automaton()
    if objects[-1] == 'b':
        for end_ind, found in auto.iter(complementary(string)):
            reversed = list(found.upper())
            reversed.reverse()
            data_context[(string_name, end_ind, end_ind + 1)] = tuple(
                ("".join(reversed), objects[0:-1], 'A', 'G'))
    elif objects == 'CG':
        for end_ind, found in auto.iter(string):
            data_context[(string_name, end_ind - 1, end_ind)] = tuple((found.upper(), 'CpG', 'T', 'C'))
            data_context[(string_name, end_ind, end_ind + 1)] = tuple((found.upper(), 'CpG', 'A', 'G'))
    elif objects == 'CHG' or objects == 'CHH':
        for end_ind, found in auto.iter(string):
            data_context[(string_name, end_ind - 2, end_ind - 1)] = tuple((found.upper(), objects, 'T', 'C'))
    elif objects == 'CN':
        for end_ind, found in auto.iter(string):
            data_context[(string_name, end_ind - 1, end_ind)] = tuple((found.upper(), 'CN', 'T', 'C'))
    elif objects == 'user':
        index_c = [user_defined_context.find('C') if user_defined_context.find('C') >= 0 else user_defined_context.find('c') if user_defined_context.find('c') >= 0 else print('The user defined context does not contain cytosines.')][0]
        for end_ind, found in auto.iter(string):
            data_context[(string_name, end_ind - index_c - 1, end_ind - index_c)] = tuple(
                (found.upper(), 'user defined context', 'T', 'C'))
        for end_ind, found in auto.iter(complementary(string)):
            reversed = list(found.upper())
            reversed.reverse()
            data_context[(string_name, end_ind - index_c, end_ind - index_c + 1)] = tuple(
                ("".join(reversed), 'user defined context', 'A', 'G'))


def context_sequence_search(context, key, fastas, string_name, user_defined_context):
    """Starts the search for cytosine contexts in the reference fasta file."""
    data_context = {}
    string = fastas[string_name]
    if key.count('C') == 0:
        for objects in key:
            ahocorasick_search(objects, context, string, string_name, user_defined_context, data_context)
    else:
        objects = "".join(key)
        ahocorasick_search(objects, context, string, string_name, user_defined_context, data_context)
    return data_context


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
    for context in list(('CAG', 'CCG', 'CTG', 'CTT', 'CCT', 'CAT', 'CTA', 'CTC', 'CAC', 'CAA', 'CCA', 'CCC')):
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

def cytosine_modification_finder(input_file, fasta_file, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, directory):
    """Searches for cytosine modification positions in the desired contexts and calculates the modificaton levels."""
    header = True
    try:
        open(input_file, 'rb')
    except (SystemExit, KeyboardInterrupt, IOError, FileNotFoundError):
        logs.error('The input file does not exist.', exc_info=True)
        sys.exit(1)
    try:
        open(fasta_file, 'r')
    except (SystemExit, KeyboardInterrupt, IOError, FileNotFoundError):
        logs.error('The genome reference fasta file does not exist.', exc_info=True)
        sys.exit(1)
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)

    if user_defined_context:
        mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'user defined context': 0}
        mean_unmod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0, 'user defined context': 0}
    else:
        mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0}
        mean_unmod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0, 'CAG': 0, 'CCG': 0, 'CTG': 0, 'CTT': 0, 'CCT': 0, 'CAT': 0, 'CTA': 0, 'CTC': 0, 'CAC': 0, 'CAA': 0, 'CCA': 0, 'CCC': 0}
    inbam = pysam.AlignmentFile(input_file, "rb")
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
        pileups = inbam.pileup(keys[i], ignore_overlaps=clip_overlap, min_base_quality=minimum_base_quality, fasta_file=fasta_file, stepper='samtools')
        for reads in pileups:
            if (reads.reference_name, reads.pos, reads.pos + 1) in modification_information_per_position:
                position = (reads.reference_name, reads.pos, reads.pos + 1)
                covered_positions.append(position)
                read_counts = defaultdict(int)
                sequences = reads.get_query_sequences()
                for pileup, seq in itertools.zip_longest(reads.pileups, sequences, fillvalue='BLANK'):
                    read_counts[(pileup.alignment.flag, seq)] += 1
                if modification_information_per_position[position][3] == 'C' and non_zero_division(read_counts[(163, 'A')] + read_counts[(83, 'A')], (read_counts[(163, 'G')] + read_counts[(83, 'G')] + read_counts[(163, 'A')] + read_counts[(83, 'A')])) < 0.8:
                    all_data = list((position[0], position[1], position[1] + 1, round(non_zero_division(read_counts[(99, 'T')] + read_counts[(147, 'T')], (read_counts[(147, 'C')] + read_counts[(99, 'C')] + read_counts[(147, 'T')] + read_counts[(99, 'T')])), 3), read_counts[(99, 'T')] + read_counts[(147, 'T')], read_counts[(147, 'C')] + read_counts[(99, 'C')], 'T', 'C', modification_information_per_position[position][0],modification_information_per_position[position][1], 'No'))
                    statistics_calculator(mean_mod, mean_unmod, all_data, user_defined_context)
                    modification_calls_writer(name, directory, all_data, header=header)
                    header = False
                elif modification_information_per_position[position][3] == 'C' and non_zero_division(read_counts[(163, 'A')] + read_counts[(83, 'A')],(read_counts[(163, 'G')] + read_counts[(83, 'G')] + read_counts[(163, 'A')] + read_counts[(83, 'A')])) >= 0.8:
                    all_data = list((position[0], position[1], position[1] + 1, round(non_zero_division(read_counts[(99, 'T')] + read_counts[(147, 'T')], (read_counts[(147, 'C')] + read_counts[(99, 'C')] + read_counts[(147, 'T')] + read_counts[(99, 'T')])), 3), read_counts[(99, 'T')] + read_counts[(147, 'T')],read_counts[(147, 'C')] + read_counts[(99, 'C')], 'T', 'C', modification_information_per_position[position][0], modification_information_per_position[position][1], 'homozyguous'))
                    statistics_calculator(mean_mod, mean_unmod, all_data, user_defined_context)
                    modification_calls_writer(name, directory, all_data, header=header)
                    header = False
                elif modification_information_per_position[position][3] == 'G' and non_zero_division(read_counts[(99, 'T')] + read_counts[(147, 'T')],(read_counts[(147, 'C')] + read_counts[(99, 'C')] + read_counts[(147, 'T')] + read_counts[(99, 'T')])) < 0.8:
                    all_data = list((position[0], position[1], position[1] + 1, round(non_zero_division(read_counts[(163, 'A')] + read_counts[(83, 'A')], (read_counts[(163, 'G')] + read_counts[(83, 'G')] + read_counts[(163, 'A')] + read_counts[(83, 'A')])), 3), read_counts[(163, 'A')] + read_counts[(83, 'A')],read_counts[(163, 'G')] + read_counts[(83, 'G')], 'A', 'G', modification_information_per_position[position][0],modification_information_per_position[position][1], 'No'))
                    statistics_calculator(mean_mod, mean_unmod, all_data, user_defined_context)
                    modification_calls_writer(name, directory, all_data, header=header)
                    header = False
                elif modification_information_per_position[position][3] == 'G' and non_zero_division(read_counts[(99, 'T')] + read_counts[(147, 'T')],(read_counts[(147, 'C')] + read_counts[(99, 'C')] + read_counts[(147, 'T')] + read_counts[(99, 'T')])) >= 0.8:
                    all_data = list((position[0], position[1], position[1] + 1, round(non_zero_division(read_counts[(163, 'A')] + read_counts[(83, 'A')], (read_counts[(163, 'G')] + read_counts[(83, 'G')] + read_counts[(163, 'A')] + read_counts[(83, 'A')])), 3), read_counts[(163, 'A')] + read_counts[(83, 'A')],read_counts[(163, 'G')] + read_counts[(83, 'G')], 'A', 'G', modification_information_per_position[position][0], modification_information_per_position[position][1], 'homozyguous'))
                    statistics_calculator(mean_mod, mean_unmod, all_data, user_defined_context)
                    modification_calls_writer(name, directory, all_data, header=header)
                    header = False
    if zero_coverage:
        positions_not_covered = [pos for pos in modification_information_per_position if pos not in covered_positions]
        positions_not_covered.sort()
        for positions in positions_not_covered:
            if modification_information_per_position[positions][3] == 'C':
                all_data = list((positions[0], positions[1], positions[1] + 1, 0, 0, 0, 'T', 'C', modification_information_per_position[positions][0], modification_information_per_position[positions][1], 'No'))
                modification_calls_writer(name, directory, all_data, header=header)
                header = False
            elif modification_information_per_position[positions][3] == 'G':
                all_data = list((positions[0], positions[1], positions[1] + 1, 0, 0, 0, 'A', 'G', modification_information_per_position[position][0], modification_information_per_position[positions][1], 'No'))
                modification_calls_writer(name, directory, all_data, header=header)
                header = False
    final_statistics_output(mean_mod, mean_unmod, directory, name, user_defined_context)
    time_e = datetime.now()
    logs.info("asTair modification finder finished running. {} seconds".format((time_e - time_b).total_seconds()))


if __name__ == '__main__':
    modification_finder_exec()


