#!/usr/bin/env python3

import re
import sys
import csv
import pdb
import gzip
import click
import pysam
import random
import logging
import warnings
from os import path
from datetime import datetime

from bam_file_parser import bam_file_opener
from DNA_sequences_operations import reverse_complementary
from simple_fasta_parser import fasta_splitting_by_sequence
from reference_context_search_triad import context_sequence_search
from reference_context_search_triad import sequence_context_set_creation


@click.command()
@click.option('fasta_file', '--fasta_file', '-f', required=True, help='DNA sequence in fasta format used for aligning the sequencing reads and mpileup.')
@click.option('read_length', '--read_length', '-l', type=int, required=True, help='Desired length of pair-end sequencing reads.')
@click.option('input_file', '--input_file', '-i', required=True, help='Sequencing reads as a Bam file or fasta sequence to generate reads.')
@click.option('simulation_input', '--simulation_input', '-si', type=click.Choice(['bam', 'fasta', 'Ns']), default='bam', required=False, help='Input file format according to the desired outcome. Bam files can be generated with other WGS simulators allowing for sequencing errors and read distributions or can be real-life sequencing data; fasta can be used to generate simple simulated data where only read positions and modification positions are known; Ns is of use to generate statistically possible sequences with changes in N nucleotides at certain positions.')
@click.option('method', '--method', '-m', required=False, default='CmtoT', type=click.Choice(['CtoT', 'CmtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and CmtoT (modified cytosines are converted to thymines, TAPS-like).')
@click.option('modification_level', '--modification_level', '-ml',  type=int, required=False, help='Desired modification level; can take any value between 0 and 100.')
@click.option('library', '--library', '-lb',  type=click.Choice(['directional', 'ptat', 'rr']), default='directional', required=False, help='Provide the correct library construction method.')
@click.option('modified_positions', '--modified_positions', '-mp', required=False, default=None, help='Provide a tab-delimited list of positions to be modified. By default the simulator randomly modifies certain positions. Please use seed for replication if no list is given.')
@click.option('context', '--context', '-co', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be modified in the output file. Default behaviour is all, which modifies positions in CpG, CHG, CHH contexts.')
@click.option('user_defined_context', '--user_defined_context', '-uc', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
@click.option('coverage', '--coverage', '-cv', required=False, help='Desired depth of sequencing coverage.')
@click.option('region', '--region', '-r', nargs=3, type=click.Tuple([str, int, int]), default=(None, None, None), required=False, help='The one-based genomic coordinates of the specific region of interest given in the form chromosome, start position, end position, e.g. chr1 100 2000.')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
@click.option('seed', '--seed', '-s', type=int, required=False, help='An integer number to be used as a seed for the random generators to ensure replication.')


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
    keys, fastas = fasta_splitting_by_sequence(fasta_file, None)
    if modified_positions is None:
        contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
        for i in range(0, len(keys)):
            modification_information = context_sequence_search(contexts, all_keys, fastas, keys[i], user_defined_context)
        return modification_information
    else:
        tupler = list()
        try:
            with open(modified_positions, newline='') as csvfile:
                position_reader = csv.reader(csvfile, delimiter='\t', lineterminator='\n')
                for row in position_reader:
                    tupler.append(tuple((str(row[0]), int(row[1]), int(row[2]))))
            return tupler
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
        random_sample = set(random.sample(modification_list_by_context, required))
    else:
        if modified_positions == None:
            random_sample = set(random.sample(modification_list_by_context, required))
    return modification_level, random_sample

##################################### write ptat and rrbs simulators


##################################################### make this probabilistic

def sequence_biaser(reads_, gc_bias, sequence_bias, read_length):
    simple_reads = list()
    for single_read in reads_:
        gc_content = (len(re.findall('C', single_read)) + len(re.findall('G', single_read)) + len(
            re.findall('c', single_read)) + len(re.findall('g', single_read))) / read_length
        low_quality = sum([1 for x in single_read if x.islower()]) / read_length
        if gc_content < gc_bias and low_quality < sequence_bias:
            simple_reads.append(single_read)
        elif gc_content <= gc_bias and low_quality < sequence_bias * 5:
            simple_reads.append(single_read)
        else:
            if random.choice([1, 4, 5]) == 4:
                simple_reads.append(single_read)
    return simple_reads


def read_sequences_generator(key, string, read_length, depth, unknown, gc_bias, sequence_bias, seed):
    sequence = list(string)
    reads = set()
    desired_depth = round((len(string) / (8 * read_length)) * depth)
    if unknown != None:
        unknown_bases = [pos for pos, char in enumerate(string) if char == 'N']
        for i in range(len(unknown_bases)):
            for k in range(0, len(string) - read_length + 1):
                element = unknown_bases[i]
                random.seed(seed)
                sequence[element] = random.choice(['A', 'G', 'T', 'C'])
                reads.add("".join(sequence)[k:k + read_length])
    else:
        for k in range(0, len(string) - read_length + 1):
            reads.add(string[k:k + read_length])
    read = [x for x in reads]
    read_ = [reverse_complementary(x) for x in reads]
    read.extend(read_)
    read = sequence_biaser(read, gc_bias, sequence_bias, read_length)
    random.seed(seed)
    read = [random.choice(read) for x in range(desired_depth)]
    return read


############################# make a function to output meaningful phred scores

def fastq_output(directory, name, key, string, read_length, depth, orientation, unknown, gc_bias, sequence_bias, seed):
    if orientation == 'F':
        name_info = "/1"
        read_info = "_R1"
        read = read_sequences_generator(key, string, read_length, depth, unknown, gc_bias, sequence_bias, seed)
    else:
        name_info = "/2"
        read_info = "_R2"
        read = [reverse_complementary(x) for x in
                read_sequences_generator(key, string, read_length, depth, unknown, gc_bias, sequence_bias, seed)]
    index = list(range(0, len(read) + 1))
    with gzip.open(directory + name + read_info + '.fq.gz', 'at') as myfile:
        for row in read:
            myfile.write('{}\n'.format("@{}_".format(key) + str(index.pop(-1)) + name_info))
            myfile.write('{}\n'.format(row))
            myfile.write('{}\n'.format('+'))
            myfile.write('{}\n'.format('I' * read_length))

############################# make a function to output a sam file


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
                line.writerow(['Absolute modified positions: ' + str(len(modified_positions_data)) + '   |   ' +
                               'Percentage to all positions of the desired context: ' + str(mod_level) + ' %'])
                line.writerow(['__________________________________________________________________________________________________'])
                for row in modified_positions_data:
                    line.writerow(row)
    else:
        pass


def bam_input_simulation(directory, name, modification_level, context, input_file, fasta_file, user_defined_context,
    modified_positions, library, seed, region, modified_positions_data, method):
    outbam = pysam.AlignmentFile(path.join(directory, name + '_' + str(modification_level) + '_' + context + ".bam"),
                                 "wb", template=bam_file_opener(input_file, None), header=True)
    modification_information = cytosine_modification_lookup(fasta_file, 'all', user_defined_context, modified_positions)
    modification_level, random_sample = random_position_modification(modification_information, modification_level,
                                                                     modified_positions, library, seed, context)
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
    if simulation_input == 'bam':
        modification_information = bam_input_simulation(directory, name, modification_level, context, input_file,
        fasta_file, user_defined_context, modified_positions, library, seed, region, modified_positions_data, method)
        absolute_modification_information(modified_positions_data, modification_information, modified_positions, name,
                                      directory, modification_level, context, True)
    else:
        keys, fastas = fasta_splitting_by_sequence(fasta_file, None)
        if simulation_input == 'Ns':
            for i in range(0, len(keys)):
                fastq_output(directory, name, keys[i], fastas[keys[i]], read_length, coverage, 'F', None, 0.5, 0.1, seed)
                fastq_output(directory, name, keys[i], fastas[keys[i]], read_length, coverage, 'R', None, 0.5, 0.1, seed)
        elif simulation_input == 'fasta':
            for i in range(0, len(keys)):
                fastq_output(directory, name, keys[i], fastas[keys[i]], read_length, coverage, 'F', 'Unknown', 0.5, 0.1, seed)
                fastq_output(directory, name, keys[i], fastas[keys[i]], read_length, coverage, 'R','Unknown', 0.5, 0.1, seed)
    time_m = datetime.now()
    logs.info("asTair's cytosine modification simulator has finished running. {} seconds".format((
    time_m - time_b).total_seconds()))


if __name__ == '__main__':
    simulator_exec()
