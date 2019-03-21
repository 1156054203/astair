#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import division

import re
import os
import sys
import csv
import pdb
import click
import pysam
import random
import logging
import warnings
from os import path
from datetime import datetime
from collections import defaultdict

from astair.bam_file_parser import bam_file_opener
from astair.simple_fasta_parser import fasta_splitting_by_sequence
from astair.context_search import context_sequence_search
from astair.context_search import sequence_context_set_creation


@click.command()
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format used for generation and modification of the sequencing reads at desired contexts.')
@click.option('read_length', '--read_length', '-l', type=int, required=True, help='Desired length of pair-end sequencing reads.')
@click.option('input_file', '--input_file', '-i', required=True, help='Sequencing reads as a BAM|CRAM file or fasta sequence to generate reads.')
@click.option('simulation_input', '--simulation_input', '-si', type=click.Choice(['bam']), default='bam', required=False, help='Input file format according to the desired outcome. BAM|CRAM files can be generated with other WGS simulators allowing for sequencing errors and read distributions or can be real-life sequencing data.')
@click.option('method', '--method', '-m', required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT)')
@click.option('modification_level', '--modification_level', '-ml',  type=int, required=False, help='Desired modification level; can take any value between 0 and 100.')
@click.option('library', '--library', '-lb',  type=click.Choice(['directional']), default='directional', required=False, help='Provide the correct library construction method. NB: Non-directional methods under development.')
@click.option('modified_positions', '--modified_positions', '-mp', required=False, default=None, help='Provide a tab-delimited list of positions to be modified. By default the simulator randomly modifies certain positions. Please use seed for replication if no list is given.')
@click.option('context', '--context', '-co', required=False, default='all', type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be modified in the output file. Default behaviour is all, which modifies positions in CpG, CHG, CHH contexts. (Default all)')
@click.option('user_defined_context', '--user_defined_context', '-uc', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
@click.option('coverage', '--coverage', '-cv', required=False, type=int, help='Desired depth of sequencing coverage.')
@click.option('region', '--region', '-r', nargs=3, type=click.Tuple([str, int, int]), default=(None, None, None), required=False, help='The one-based genomic coordinates of the specific region of interest given in the form chromosome, start position, end position, e.g. chr1 100 2000.')
@click.option('user_defined_context', '--user_defined_context', '-uc', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
@click.option('overwrite', '--overwrite', '-ov', required=False, default=False, type=bool, help='Indicates whether existing output files with matching names will be overwritten. (Default False)')
#@click.option('per_chromosome', '--per_chromosome', '-chr', default=None, type=str, help='When used, it modifies the chromosome given only. (Default None')
@click.option('GC_bias', '--GC_bias', '-gc', default=0.3, required=True, type=float, help='The value of total GC levels in the read above which lower coverage will be observed in Ns and fasta modes. (Default 0.5)')
@click.option('sequence_bias', '--sequence_bias', '-sb', default=0.1, required=True, type=float, help='The proportion of lower-case letters in the read string for the Ns and fasta modes that will decrease the chance of the read being output. (Default 0.1)')
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
@click.option('seed', '--seed', '-s', type=int, required=False, help='An integer number to be used as a seed for the random generators to ensure replication.')


def simulate(reference, read_length, input_file, method, library, simulation_input, modification_level,
                   modified_positions, coverage, context, region, directory, seed, user_defined_context, N_threads, GC_bias, sequence_bias, overwrite):
    """Simulate TAPS/BS conversion on top of an existing bam/cram file."""
    modification_simulator(reference, read_length, input_file, method, library, simulation_input, modification_level,
              modified_positions, coverage, context, region, directory, seed, user_defined_context, N_threads, GC_bias, sequence_bias, overwrite)



warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

#logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()


def cytosine_modification_lookup(reference, context, user_defined_context, modified_positions, region):
    """Finds all required cytosine contexts or takes positions from a tab-delimited file containing
     the list of positions to be modified."""
    keys, fastas = fasta_splitting_by_sequence(reference, None)
    context_total_counts = defaultdict(int)
    if modified_positions is None:
        contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
        if region is None:
            for i in range(0, len(keys)):
                modification_information = context_sequence_search(contexts, all_keys, fastas, keys[i], user_defined_context, context_total_counts, region)
        else:
            if region[0] in keys:
                 modification_information = context_sequence_search(contexts, all_keys, fastas, region[0], user_defined_context, context_total_counts, region)
        try:
            return modification_information
        except UnboundLocalError:
            logs.error('There is no reference sequence of this name in the provided fasta file.', exc_info=True)
            sys.exit(1)
    else:
        tupler = list()
        try:
            with open(modified_positions) as csvfile:
                position_reader = csv.reader(csvfile, delimiter='\t', lineterminator='\n')
                for row in position_reader:
                    tupler.append(tuple((str(row[0]), int(row[1]), int(row[2]))))
            return tupler
        except Exception:
            logs.error('The cytosine positions file does not exist.', exc_info=True)
            sys.exit(1)


def modification_level_transformation(modification_level, modified_positions):
    """Transforms the user-provided modification level to a float, string or bool."""
    if modification_level != 0 and modified_positions is None:
         modification_level = modification_level / 100
    elif modified_positions is not None:
        modification_level = 'user_provided_list'
    else:
        modification_level = None
    return modification_level


def random_position_modification(modification_information, modification_level, modified_positions, library, seed, context):
    """Creates lists of positions per context that will be modified according to the method."""
    modification_level = modification_level_transformation(modification_level, modified_positions)
    if modification_level == None:
        modification_level = 0
    if context == 'all' and modified_positions == None:
        modification_list_by_context = set()
        all_keys = list(('CHG','CHH','CpG'))
        for context_string in all_keys:
            modification_list_by_context = modification_list_by_context.union(set((keys) for keys, vals in modification_information.items() if vals[1] == context_string))
        required = round((len(modification_list_by_context)) * modification_level)
    elif context != 'all' and modified_positions == None:
        modification_list_by_context = set((keys) for keys, vals in modification_information.items() if vals[1] == context)
        required = round((len(modification_list_by_context)) * modification_level)
    else:
        random_sample = set(modification_information)
        modification_level = 'custom'
    if seed is not None and modified_positions == None:
        random.seed(seed)
        random_sample = set(random.sample(modification_list_by_context, int(required)))
    else:
        if modified_positions == None:
            random_sample = set(random.sample(modification_list_by_context, int(required)))
    return modification_level, random_sample


def general_read_information_output(name, directory, read, modification_level, header, region, method, context):
    """Writes to fastq file."""
    if read.is_read1 == True:
        orientation = '/1'
    else:
        orientation = '/2'
    try:
        if not isinstance(modification_level, str):
            modification_level = int(modification_level*100)
        with open(path.join(directory, name + '_' + method + '_' + str(modification_level) + '_' + context + '_read_information.txt'), 'a') as reads_info_output:
            line = csv.writer(reads_info_output, delimiter='\t', lineterminator='\n')
            if header == True:
                line.writerow(['Read ID', 'reference', 'start', 'end'])
                line.writerow([read.qname + orientation, read.reference_name, read.reference_start,
                                    read.reference_start + read.query_length])
            else:
                line.writerow([read.qname + orientation, read.reference_name, read.reference_start,
                                    read.reference_start + read.query_length])
    except IOError:
        logs.error('asTair cannot write to read information file.', exc_info=True)

def modification_by_strand(read, library):
    """Outputs read positions that may be modified."""
    if read.flag == 99 or read.flag == 147 and library == 'directional':
        positC = [val.start() + read.reference_start for val in re.finditer('C', read.query_sequence)]
        positions = set((read.reference_name, pos, pos + 1) for pos in positC)
        base = 'T'
        return positions, base
    elif read.flag == 83 or read.flag == 163 and library == 'directional':
        positG = [val.start() + read.reference_start for val in re.finditer('G', read.query_sequence)]
        positions = set((read.reference_name, pos, pos + 1) for pos in positG)
        base = 'A'
    return positions, base


def absolute_modification_information(modified_positions_data, modification_information, modified_positions, name, directory, modification_level, context, method):
    """Gives a statistics summary file about the modified positions."""
    modified_positions_data = set(modified_positions_data)
    modified_positions_data = list(modified_positions_data)
    modified_positions_data.sort()
    if modified_positions is None and modification_level is not None:
        if context == 'all':
            context_list_length = len(modification_information)
        else:
            context_list_length = len(set(keys for keys, vals in modification_information.items() if vals[1] == context))
        mod_level = round((len(modified_positions_data) / context_list_length) * 100, 3)
        try:
            with open(path.join(directory, name + '_' + method + '_' + str(modification_level) + '_' + context + '_modified_positions_information.txt'), 'w') as reads_info_output:
                line = csv.writer(reads_info_output, delimiter='\t', lineterminator='\n')
                line.writerow(['__________________________________________________________________________________________________'])
                line.writerow(['Absolute modified positions: ' + str(len(modified_positions_data)) + '   |   ' +
                               'Percentage to all positions of the desired context: ' + str(mod_level) + ' %'])
                line.writerow(['__________________________________________________________________________________________________'])
                for row in modified_positions_data:
                    line.writerow(row)
        except IOError:
            logs.error('asTair cannot write to modified positions summary file.', exc_info=True)
    else:
        pass


def bam_input_simulation(directory, name, modification_level, context, input_file, reference, user_defined_context,
    modified_positions, library, seed, region, modified_positions_data, method, N_threads, header, overwrite, extension):
    """Inserts modification information acording to method and context to a bam or cram file."""
    if not os.path.isfile(path.join(directory, name + '_' + method + '_' + str(modification_level) + '_' + context + extension)) or overwrite is True:
        if pysam.AlignmentFile(input_file).is_cram:
            outbam = pysam.AlignmentFile(path.join(directory, name + '_' + method + '_' + str(modification_level) + '_' + context + extension),
                "wc", reference_filename=reference, template=bam_file_opener(input_file, None, N_threads), header=header)
        else:
            outbam = pysam.AlignmentFile(path.join(directory, name + '_' + method + '_' + str(modification_level) + '_' + context + extension),
                                     "wb", template=bam_file_opener(input_file, None, N_threads), header=header)
        modification_information = cytosine_modification_lookup(reference, 'all', user_defined_context, modified_positions, region)
        modification_level, random_sample = random_position_modification(modification_information, modification_level,
                                                                         modified_positions, library, seed, context)

        for read in bam_file_opener(input_file, 'fetch', N_threads):
            general_read_information_output(name, directory, read, modification_level, header, region, method, context)
            quals = read.query_qualities
            if read.flag in [99, 147, 83, 163]:
                positions, base = modification_by_strand(read, library)
                modified_positions_data.extend(list(random_sample.intersection(positions)))
                if method == 'mCtoT':
                    if random_sample.intersection(positions):
                        indices = [position[1] - read.reference_start for position in random_sample.intersection(positions)]
                        strand = list(read.query_sequence)
                        replace = list(base * len(indices))
                        for (index, replacement) in zip(indices, replace):
                            strand[index] = replacement
                        read.query_sequence = "".join(strand)
                        read.query_qualities = quals
                        outbam.write(read)
                    else:
                        outbam.write(read)

                else:
                    indices = [position[1] - read.reference_start for position in positions if
                               position not in random_sample.intersection(positions)]
                    strand = list(read.query_sequence)
                    replace = list(base * len(indices))
                    for (index, replacement) in zip(indices, replace):
                        strand[index] = replacement
                    read.query_sequence = "".join(strand)
                    read.query_qualities = quals
                    outbam.write(read)
                header = False
        return modification_information


def modification_simulator(reference, read_length, input_file, method, library, simulation_input, modification_level,
                           modified_positions, coverage, context, region, directory, seed, user_defined_context, N_threads, GC_bias, sequence_bias, overwrite):
    "Assembles the whole modification simulator and runs per mode, method, library and context."
    header = True
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    if region.count(None)!=0:
        region = None
    modified_positions_data = list()
    if simulation_input == 'bam':
        if pysam.AlignmentFile(input_file).is_cram:
            extension = '.cram'
        else:
            extension = '.bam'
        try:
            modification_information = bam_input_simulation(directory, name, modification_level, context, input_file,
            reference, user_defined_context, modified_positions, library, seed, region, modified_positions_data, method, N_threads, header, overwrite, extension)
            absolute_modification_information(modified_positions_data, modification_information, modified_positions,name,directory, modification_level, context, method)
            pysam.index(path.join(directory, name + '_' + method + '_' + str(modification_level) + '_' + context + extension))
        except AttributeError:
            logs.error(
                'The output files will not be overwritten. Please rename the input or the existing output files before rerunning if the input is different.',
                exc_info=True)
    time_m = datetime.now()
    logs.info("asTair's cytosine modification simulator has finished running. {} seconds".format((
    time_m - time_b).total_seconds()))


if __name__ == '__main__':
    simulate()
