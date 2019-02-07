#!/usr/bin/env python3


import random
import sys
import logging
import warnings
import re
import csv
from os import path
import click
import pysam
from datetime import datetime

from bam_file_parser import bam_file_opener
from simple_fasta_parser import fasta_splitting_by_sequence
from reference_context_search_triad import sequence_context_set_creation
from reference_context_search_triad import context_sequence_search


@click.command()
@click.option('fasta_file', '--fasta_file', required=True, help='DNA sequence in fasta format used for aligning the sequencing reads and mpileup.')
@click.option('read_length', '--read_length', type=int, required=True, help='Desired length of pair-end sequencing reads.')
@click.option('input_file', '--input_file', required=True, help='Sequencing reads as a Bam file or fasta sequence to generate reads.')
@click.option('simulation_input', '--simulation_input', type=click.Choice(['bam', 'fasta', 'Ns']), default='bam', required=False, help='Input file format according to the desired outcome. Bam files can be generated with other WGS simulators allowing for sequencing errors and read distributions or can be real-life sequencing data; fasta can be used to generate simple simulated data where only read positions and modification positions are known; Ns is of use to generate statistically possible sequences with changes in N nucleotides at certain positions.')
@click.option('method', '--method', required=False, default = 'CmtoT', type=click.Choice(['CtoT', 'CmtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and CmtoT (modified cytosines are converted to thymines, TAPS-like).')
@click.option('modification_level', '--modification_level', type=int, required=False, help='Desired modification level; can take any value between 0 and 100.')
@click.option('library', '--library', type=click.Choice(['directional', 'ptat']), default='directional', required=False, help='Provide the correct library construction method.')
@click.option('modified_positions', '--modified_positions', required=False, default=None, help='Provide a tab-delimited list of positions to be modified. By default the simulator randomly modifies certain positions. Please use seed for replication if no list is given.')
@click.option('context', '--context', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be modified in the output file. Default behaviour is all, which modifies positions in CpG, CHG, CHH contexts.')
@click.option('user_defined_context', '--user_defined_context', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
@click.option('coverage', '--coverage', required=False, help='Desired depth of sequencing coverage.')
@click.option('region', '--region', nargs=3, type=click.Tuple([str, int, int]), default=(None, None, None), required=False, help='The one-based genomic coordinates of the specific region of interest given in the form chromosome, start position, end position, e.g. chr1 100 2000.')
@click.option('directory', '--directory', required=True, type=str, help='Output directory to save files.')
@click.option('seed', '--seed', type=int, required=False, help='An integer number to be used as a seed for the random generators to ensure replication.')


def simulator_exec(fasta_file, read_length, input_file, method, library, simulation_input, modification_level,
                   modified_positions, coverage, context, region, directory, seed, user_defined_context):
    modification_simulator(fasta_file, read_length, input_file, method, library, simulation_input, modification_level,
              modified_positions, coverage, context, region, directory, seed, user_defined_context)



warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()


def cytosine_modification_lookup(fasta_file, context, user_defined_context, modified_positions):
    keys, fastas = fasta_splitting_by_sequence(fasta_file)
    if modified_positions is None:
        contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
        for i in range(0, len(keys)):
            modification_information = context_sequence_search(contexts, all_keys, fastas, keys[i], user_defined_context)
        return modification_information
    else:
        tupler, list_to_use = list(), set()
        try:
            with open(modified_positions, newline='') as csvfile:
                position_reader = csv.reader(csvfile, delimiter='\t', lineterminator='\n')
                for row in position_reader:
                    tupler.append(row)
            list_to_use.union(set((x[0], int(x[1]), int(x[2])) for x in tupler))
            return list_to_use
        except (SystemExit, KeyboardInterrupt, IOError, FileNotFoundError):
            logs.error('The cytosine positions file does not exist.', exc_info=True)
            sys.exit(1)


def modification_level_transformation(modification_level, modified_positions):
    if modification_level != 0 and modified_positions is None:
         modification_level = modification_level / 100
    elif modified_positions is not None:
        modification_level = 'user_provided_list'
    else:
        modification_level = None
    return modification_level


def random_position_modification(modification_information, modification_level, modified_positions, library, seed, context):
    modification_level = modification_level_transformation(modification_level, modified_positions)
    if context == 'all':
        modification_list_by_context = set()
        all_keys = list(('CHG','CHH','CpG'))
        for context_string in all_keys:
            modification_list_by_context = modification_list_by_context.union(set((keys) for keys, vals in modification_information.items() if vals[1] == context_string))
    else:
        modification_list_by_context = set((keys) for keys, vals in modification_information.items() if vals[1] == context)
    required = round((len(modification_list_by_context)) * modification_level)
    if seed is not None:
        random.seed(seed)
        random_sample = set(random.sample(modification_list_by_context, required))
    else:
        random_sample = set(random.sample(modification_list_by_context, required))
    return modification_level, random_sample


def general_read_information_output(name, directory, read, modification_level, header, region):
    if read.is_read1 == True:
        orientation = '/1'
    else:
        orientation = '/2'
    with open(directory + name + '_read_information_' + str(modification_level) + '.txt', 'a',
              newline='') as reads_info_output:
        line = csv.writer(reads_info_output, delimiter='\t', lineterminator='\n')
        if header == True:
            line.writerow(['Read ID', 'reference', 'start', 'end'])
            if region is not None:
                line.writerow([read.qname + orientation, read.reference_name, read.reference_start + region,
                                read.reference_start + read.query_length + region])
            else:
                line.writerow([read.qname + orientation, read.reference_name, read.reference_start,
                                read.reference_start + read.query_length])
        else:
            if region is not None:
                line.writerow([read.qname + orientation, read.reference_name, read.reference_start + region,
                                read.reference_start + read.query_length + region])
            else:
                line.writerow([read.qname + orientation, read.reference_name, read.reference_start,
                                read.reference_start + read.query_length])


def modification_by_strand(read, library):
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


def absolute_modification_information(modified_positions_data, modification_information, modified_positions, name, directory, modification_level, context, header):
    modified_positions_data = set(modified_positions_data)
    modified_positions_data = list(modified_positions_data)
    modified_positions_data.sort()
    if modified_positions is None and modification_level is not None:
        if context == 'all':
            context_list_length = len(modification_information)
        else:
            context_list_length = len(set(keys for keys, vals in modification_information.items() if vals[1] == context))
        mod_level = round((len(modified_positions_data) / context_list_length) * 100, 3)
        with open(path.join(directory, name + '_modified_' + 'read_information_' + str(modification_level) + "_" + context + '.txt'), 'w',
                  newline='') as reads_info_output:
            line = csv.writer(reads_info_output, delimiter='\t', lineterminator='\n')
            if header == True:
                line.writerow(['__________________________________________________________________________________________________'])
                line.writerow(['Absolute modified positions: ' + str(len(modified_positions_data)) + '   |   ' + 'Percentage to all positions of the desired context: ' + str(mod_level) + ' %'])
                line.writerow(['__________________________________________________________________________________________________'])
                for row in modified_positions_data:
                    line.writerow(row)
    else:
        pass


def modification_simulator(fasta_file, read_length, input_file, method, library, simulation_input, modification_level,
                           modified_positions, coverage, context, region, directory, seed, user_defined_context):
    header = True
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    if region.count(None)==0:
        region = region[2]-region[1]
    else:
        region = None
    modified_positions_data = list()
    outbam = pysam.AlignmentFile(path.join(directory, name + '_' + str(modification_level) + ".bam"), "wb", template=bam_file_opener(input_file, None), header=True)
    modification_information = cytosine_modification_lookup(fasta_file, 'all', user_defined_context, modified_positions)
    modification_level, random_sample = random_position_modification(modification_information, modification_level, modified_positions, library, seed, context)
    for read in bam_file_opener(input_file, 'fetch'):
        general_read_information_output(name, directory, read, modification_level, header, region)
        quals = read.query_qualities
        if read.flag in [99, 147, 83, 163]:
            positions, base = modification_by_strand(read, library)
            modified_positions_data.extend(list(random_sample.intersection(positions)))
            if method == 'CmtoT':
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
                indices = [position[1] - read.reference_start for position in positions if position not in random_sample.intersection(positions)]
                strand = list(read.query_sequence)
                replace = list(base * len(indices))
                for (index, replacement) in zip(indices, replace):
                    strand[index] = replacement
                read.query_sequence = "".join(strand)
                read.query_qualities = quals
                outbam.write(read)
            header = False
        else:
            pass
    absolute_modification_information(modified_positions_data, modification_information, modified_positions, name,
                                      directory, modification_level, context, True)
    time_m = datetime.now()
    logs.info("asTair's cytosine modification simulator has finished running. {} seconds".format((
    time_m - time_b).total_seconds()))


if __name__ == '__main__':
    simulator_exec()
