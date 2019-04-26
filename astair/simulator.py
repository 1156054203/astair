#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import division

import re
import os
import sys
import csv
import pdb
import mmap
import click
import pysam
import numpy
import random
import logging
import warnings
import itertools
from os import path
from datetime import datetime
from collections import defaultdict
from threading import Thread as Thread

if sys.version[0] == '3':
    from queue import Queue as Queue
elif sys.version[0] == '2':
    from Queue import Queue as Queue
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))

from astair.bam_file_parser import bam_file_opener
from astair.simple_fasta_parser import fasta_splitting_by_sequence
from astair.context_search import context_sequence_search
from astair.context_search import sequence_context_set_creation


@click.command()
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format used for generation and modification of the sequencing reads at desired contexts.')
@click.option('read_length', '--read_length', '-l', type=int, required=True, help='Desired length of pair-end sequencing reads.')
@click.option('input_file', '--input_file', '-i', required=True, help='Sequencing reads as a BAM|CRAMfile or fasta sequence to generate reads.')
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
@click.option('overwrite', '--overwrite', '-ov', required=False, default=False, is_flag=True, help='Indicates whether existing output files with matching names will be overwritten. (Default False)')
@click.option('per_chromosome', '--per_chromosome', '-chr', default=None, type=str, help='When used, it modifies the chromosome given only. (Default None')
@click.option('GC_bias', '--GC_bias', '-gc', default=0.3, required=True, type=float, help='The value of total GC levels in the read above which lower coverage will be observed in Ns and fasta modes. (Default 0.5)')
@click.option('sequence_bias', '--sequence_bias', '-sb', default=0.1, required=True, type=float, help='The proportion of lower-case letters in the read string for the Ns and fasta modes that will decrease the chance of the read being output. (Default 0.1)')
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('reverse_modification', '--rev', '-rv', default=False, is_flag=True, required=False, help='Returns possible or known modified position to their unmodified expected state. NB: Works only on files with MD tags (Default False).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
@click.option('seed', '--seed', '-s', type=int, required=False, help='An integer number to be used as a seed for the random generators to ensure replication.')


def simulate(reference, read_length, input_file, method, library, simulation_input, modification_level,
                   modified_positions, coverage, context, region, directory, seed, user_defined_context, N_threads, per_chromosome, GC_bias, sequence_bias, overwrite, reverse_modification):
    """Simulate TAPS/BS conversion on top of an existing bam/cram file."""
    modification_simulator(reference, read_length, input_file, method, library, simulation_input, modification_level,
              modified_positions, coverage, context, region, directory, seed, user_defined_context, N_threads, per_chromosome, GC_bias, sequence_bias, overwrite, reverse_modification)



warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

#logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def csv_line_skipper(csvfile, start, end, count, key):
    """Random access for tab-delimited file reading."""
    if count == 0:
        list_of_rows = list()
        csvfile.seek(start)
        for read in csvfile.readlines():
            if read.split()[3].lower().islower() == False and float(read.split()[3]) != 0 and read.split()[0]==key:
                list_of_rows.append(tuple(read.split()))
    else:
        csvfile.seek(start)
        list_of_rows = list(csv.reader(csvfile, delimiter='\t', lineterminator='\n'))
    return list_of_rows


def simple_reader(tupler, line, input_file, keys):
    """Simplified function used for tab-delimited file reading."""
    inbam = bam_file_opener(input_file, None, 1)
    for row in line:
        if row[0] != 'CHROM' and row[0] == keys and float(row[3]) != 0:
            tupler[tuple((str(row[0]), int(row[1]), int(row[2])))] = numpy.array(object=[0, float(row[3]) * len([i.flag for i in inbam.fetch(contig=str(row[0]), start=int(row[1]), stop=int(row[2])) if (i.flag in [99, 147] and row[7] in ['C', 'T']) or (i.flag in [163, 83] and row[7] in ['A', 'G'])])], dtype=numpy.int8, copy=False)


def chunck_sizer(input_list, threadN, input_file, keys, tupler):
    """Threading for reading of larger files."""
    chunk_list = list()
    queue = Queue(maxsize=0)
    chunk = round(len(input_list)/threadN)
    for i in range(0, threadN):
        chunk_list.append(input_list[i*chunk:i*chunk+chunk])
    threads = [Thread(target=simple_reader, args=(tupler, line, input_file, keys), ) for line in chunk_list]
    for thread in threads:
        thread.daemon = True
        thread.start()
        queue.put(thread, block=False, timeout=10)
        thread.join()
        queue.task_done()
    return queue


def cytosine_modification_lookup(context, user_defined_context, modified_positions, region, fastas, keys, context_total_counts, input_file, N_threads, csvfile):
    """Finds all required cytosine contexts or takes positions from a tab-delimited file containing
     the list of positions to be modified."""
    if modified_positions is None:
        contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
        if region == None:
            modification_information = context_sequence_search(contexts, all_keys, fastas, keys, user_defined_context, context_total_counts, region)
        else:
            if region[0] in keys:
                 modification_information = context_sequence_search(contexts, all_keys, fastas, region[0], user_defined_context, context_total_counts, region)
        try:
            return modification_information
        except UnboundLocalError:
            logs.error('There is no reference sequence of this name in the provided fasta file.', exc_info=True)
            sys.exit(1)
    else:
        tupler = dict()
        try:
            memory_map = mmap.mmap(csvfile.fileno(), 0)
            start = memory_map.find(keys.encode('utf8'))
            bam_to_view = bam_file_opener(input_file, None, 1)
            count = bam_to_view.header.references.index(keys)
            end = memory_map.rfind(keys.encode('utf8'))
            position_reader = csv_line_skipper(csvfile, start, end, count, keys)
            queue = chunck_sizer(position_reader, N_threads, input_file, keys, tupler)
            queue.join()
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


def general_read_information_output(read, header, line):
    """Writes to fastq file."""
    if read.is_read1 == True:
        orientation = '/1'
    else:
        orientation = '/2'
    try:
        if header == True:
            line.writerow(['Read ID', 'reference', 'start', 'end'])
            line.writerow([read.qname + orientation, read.reference_name, read.reference_start,
                                read.reference_start + read.query_length])
        else:
            line.writerow([read.qname + orientation, read.reference_name, read.reference_start,
                                read.reference_start + read.query_length])
    except IOError:
        logs.error('asTair cannot write to read information file.', exc_info=True)


def cigar_search(read_data):
    """Looks whether there are indels, soft clipping or pads the CIGAR string"""
    changes = [int(s) for s in re.findall(r'\d+', read_data)]
    non_overlap = [x + 1 if x == 0 else x for x in changes]
    names = list(re.findall(r'[^\W\d_]+', read_data))
    positions = [x for x in [sum(non_overlap[0:i]) for i in range(1, len(non_overlap)+1)]]
    return names, positions, changes


def position_correction_cigar(read, method, random_sample, positions, reverse_modification):
    """Uses the CIGAR string information to correct the expected cytosine positions."""
    names, positions_cigar, changes = cigar_search(read.cigarstring)
    index = 0
    for change in names:
        if len(positions) != 0:
            if change == 'D':
                if isinstance(list(positions)[0], tuple):
                    subsample = random_sample.intersection(positions)
                    if method == 'CtoT' and reverse_modification == False:
                        corrected_positions = [x[1] - abs(read.qstart - read.reference_start) if (x[1] - abs(read.qstart - read.reference_start)) < positions_cigar[index] else x[1] - abs(read.qstart - read.reference_start) - changes[index] for x in positions if x not in subsample]
                    else:
                        corrected_positions = [x[1] - abs(read.qstart - read.reference_start) if (x[1] - abs(
                            read.qstart - read.reference_start)) < positions_cigar[index] else x[1] - abs(
                            read.qstart - read.reference_start) - changes[index] for x in subsample]
                else:
                    corrected_positions = [x if x < positions_cigar[index] else x - changes[index] for x in positions]
                index += 1
                positions = corrected_positions
            elif change == 'I':
                if isinstance(list(positions)[0], tuple):
                    subsample = random_sample.intersection(positions)
                    if method == 'CtoT' and reverse_modification == False:
                        corrected_positions = [x[1] - abs(read.qstart-read.reference_start) if (x[1] - abs(read.qstart-read.reference_start)) < positions_cigar[index] else x[1] - abs(read.qstart-read.reference_start) + changes[index] for x in positions if x not in subsample]
                    else:
                        corrected_positions = [x[1] - abs(read.qstart - read.reference_start) if (x[1] - abs(
                            read.qstart - read.reference_start)) < positions_cigar[index] else x[1] - abs(
                            read.qstart - read.reference_start) + changes[index] for x in subsample]
                else:
                    corrected_positions = [x if x < positions_cigar[index] else x + changes[index] for x in positions]
                index += 1
                positions = corrected_positions
            elif change == 'S':
                if isinstance(list(positions)[0], tuple):
                    subsample = random_sample.intersection(positions)
                    if method == 'CtoT' and reverse_modification == False:
                        corrected_positions = [x[1] - abs(read.qstart - read.reference_start) for x in positions if
                                               x not in subsample]
                    else:
                        corrected_positions = [x[1] - abs(read.qstart-read.reference_start) for x in subsample]
                else:
                    if index == 0:
                        corrected_positions = [x for x in positions if x > positions_cigar[index]]
                    else:
                        corrected_positions = [x for x in positions if x < positions_cigar[index]]
                index += 1
                positions = corrected_positions
            else:
                index += 1
    return positions


def modification_by_strand(read, library, reverse_modification, fastas):
    """Outputs read positions that may be modified."""
    if library == 'directional':
        if reverse_modification == False:
            if read.flag in [99,147]:
                base, ref = 'T', 'C'
            elif read.flag in [83, 163]:
                base, ref = 'A', 'G'
        elif reverse_modification == True:
            if read.flag in [99,147]:
                base, ref = 'C', 'C'
            elif read.flag in [83, 163]:
                base, ref = 'G', 'G'
        posit = [val.start() + read.reference_start for val in re.finditer(ref, fastas[read.reference_name][read.reference_start:read.reference_start+read.qlen].upper())]
        positions = set((read.reference_name, pos, pos + 1) for pos in posit)
        return positions, base



def absolute_modification_information(modified_positions_data, modification_information, modified_positions, name, directory, modification_level, context, method, per_chromosome, line):
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
    else:
        mod_level = 'Custom'
    try:
        line.writerow(['__________________________________________________________________________________________________'])
        line.writerow(['Absolute modified positions: ' + str(len(modified_positions_data)) + '   |   ' +
                       'Percentage to all positions of the desired context: ' + str(mod_level) + ' %'])
        line.writerow(['__________________________________________________________________________________________________'])
        for row in modified_positions_data:
            line.writerow(row)
    except IOError:
        logs.error('asTair cannot write to modified positions summary file.', exc_info=True)



def modification_information_and_reads_fetching(context, user_defined_context, modified_positions, region, fastas, key, context_total_counts, modification_level, library, seed, input_file, N_threads, csvfile):
    """Prepares the required cytosine positions vectors."""
    modification_information = cytosine_modification_lookup(context, user_defined_context, modified_positions, region, fastas, key, context_total_counts, input_file, N_threads, csvfile)
    modification_level, random_sample = random_position_modification(modification_information, modification_level, modified_positions, library, seed, context)
    return modification_information, random_sample, modification_level



def read_modification(input_file, fetch, N_threads, name, directory, modification_level, header, region, method, context, modified_positions_data, random_sample, fastas, library, reverse_modification, outbam, line, modified_positions, modification_information):
    """Looks at the sequencing reads and modifies/removes modifications at the required cytosine positions."""
    for read in bam_file_opener(input_file, fetch, N_threads):
        general_read_information_output(read, header, line)
        quals = read.query_qualities
        if read.flag in [99, 147, 83, 163] and read.reference_length != 0:
            positions, base = modification_by_strand(read, library, reverse_modification, fastas)
            modified_positions_data.extend(list(random_sample.intersection(positions)))
            if method == 'mCtoT':
                if len(read.tags) != 0 and (
                    re.findall('I', read.cigarstring, re.IGNORECASE) or re.findall('D', read.cigarstring,
                                                                                   re.IGNORECASE)) or re.findall('S',
                                                                                                                 read.cigarstring,
                                                                                                                 re.IGNORECASE):
                    indices = position_correction_cigar(read, method, random_sample, positions, reverse_modification)
                else:
                    indices = [position[1] - read.reference_start for position in random_sample.intersection(positions)]
                if len(indices) > 0:
                    strand = list(read.query_sequence)
                    indices.sort()
                    if modified_positions:
                        not_changed = list(random_sample.intersection(positions))
                        not_changed.sort()
                        for index in range(0, len(not_changed)):
                           if modification_information[not_changed[index]][0] <= modification_information[not_changed[index]][1] and len(indices) > index:
                               if len(strand) > indices[index] and (
                                   ((strand[index] in ['C', 'T'] and read.flag in [99, 147])) or (
                                       (strand[index] in ['G', 'A'] and read.flag in [83, 163]))):
                                    strand[indices[index]] = base
                                    modification_information[not_changed[index]][0] += 1

                    else:
                        replace = list(base * len(indices))
                        for (index, replacement) in zip(indices, replace):
                            if len(strand) > index:
                                if ((strand[index] in ['C', 'T'] and read.flag in [99, 147])) or (
                                (strand[index] in ['G', 'A'] and read.flag in [83, 163])):
                                    strand[index] = replacement
                    read.query_sequence = "".join(strand)
                    read.query_qualities = quals
                    outbam.write(read)
                else:
                    outbam.write(read)
            else:
                if len(read.tags) != 0 and (
                    re.findall('I', read.cigarstring, re.IGNORECASE) or re.findall('D', read.cigarstring,
                                                                                   re.IGNORECASE)) or re.findall('S',
                                                                                                                 read.cigarstring,
                                                                                                                 re.IGNORECASE):
                    indices = position_correction_cigar(read, method, random_sample, positions, reverse_modification)
                else:
                    if reverse_modification == True:
                        indices = [position[1] - read.reference_start for position in random_sample.intersection(positions)]
                    else:
                        indices = [position[1] - read.reference_start for position in positions if
                                   position not in random_sample.intersection(positions)]
                if len(indices) > 0:
                    strand = list(read.query_sequence)
                    indices.sort()
                    if modified_positions:
                        not_changed = list(random_sample.intersection(positions))
                        not_changed.sort()
                        for index in range(0, len(not_changed)):
                            if modification_information[not_changed[index]][0] <= modification_information[not_changed[index]][1] and len(indices) > index:
                                if len(strand) > indices[index] and (
                                    ((strand[index] in ['C', 'T'] and read.flag in [99, 147])) or (
                                        (strand[index] in ['G', 'A'] and read.flag in [83, 163]))):
                                    strand[indices[index]] = base
                                    modification_information[not_changed[index]][0] += 1
                    else:
                        replace = list(base * len(indices))
                        for (index, replacement) in zip(indices, replace):
                            if len(strand) > index:
                                if ((strand[index] in ['C', 'T'] and read.flag in [99, 147])) or (
                                (strand[index] in ['G', 'A'] and read.flag in [83, 163])):
                                    strand[index] = replacement
                    read.query_sequence = "".join(strand)
                    read.query_qualities = quals
                    outbam.write(read)
                else:
                    outbam.write(read)
            header = False
        elif read.flag in [99, 147, 83, 163]:
            outbam.write(read)


def bam_input_simulation(directory, name, modification_level, context, input_file, reference, user_defined_context, per_chromosome,
    modified_positions, library, seed, region, modified_positions_data, method, N_threads, header, overwrite, extension, reverse_modification):
    """Inserts modification information acording to method and context to a bam or cram file."""
    if not os.path.isfile(path.join(directory, name + '_' + method + '_' + str(modification_level) + '_' + context + extension)) or overwrite is True:
        if pysam.AlignmentFile(input_file).is_cram:
            file_type = 'wc'
        else:
            file_type = 'wb'
        if not isinstance(modification_level, str):
            modification_level_= int(modification_level*100)
        else:
            modification_level_ = modification_level
        if modified_positions:
            csvfile = open(modified_positions, 'r+')
        else:
            csvfile = None
        if reverse_modification == False:
            if per_chromosome == None:
                outbam = pysam.AlignmentFile(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + extension),
                file_type, reference_filename=reference, template=bam_file_opener(input_file, None, N_threads), header=header)
            else:
                outbam = pysam.AlignmentFile(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context  + '_' + per_chromosome  + extension),
                file_type, reference_filename=reference, template=bam_file_opener(input_file, None, N_threads), header=header)
        else:
            if per_chromosome == None:
                outbam = pysam.AlignmentFile(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + '_reversed' + extension),
                file_type, reference_filename=reference, template=bam_file_opener(input_file, None, N_threads), header=header)
            else:
                outbam = pysam.AlignmentFile(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context  + '_reversed_' + per_chromosome  + extension),
                file_type, reference_filename=reference, template=bam_file_opener(input_file, None, N_threads), header=header)
        keys, fastas = fasta_splitting_by_sequence(reference, per_chromosome)
        if per_chromosome == None:
            name_to_use = path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + '_read_information.txt')
            name_to_use_absolute = path.join(directory,name + '_' + method + '_' + str(modification_level_) + '_' + context + '_modified_positions_information.txt')
        else:
            name_to_use = path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + per_chromosome + '_' + context + '_read_information.txt')
            name_to_use_absolute = path.join(directory,name + '_' + method + '_' + str(modification_level_) + '_' + context + '_' + per_chromosome + '_modified_positions_information.txt')
        reads_info_output = open(name_to_use, 'w')
        line = csv.writer(reads_info_output, delimiter='\t', lineterminator='\n')
        reads_info_absolute = open(name_to_use_absolute, 'w')
        line_ = csv.writer(reads_info_absolute, delimiter='\t', lineterminator='\n')
        context_total_counts = defaultdict(int)
        if region == None and per_chromosome == None:
            for i in range(0, len(keys)):
                modified_positions_data = list()
                modification_information, random_sample, modification_level = modification_information_and_reads_fetching(context, user_defined_context, modified_positions, region, fastas, keys[i], context_total_counts, modification_level, library, seed, input_file, N_threads, csvfile)
                fetch = tuple((keys[i], 0, pysam.AlignmentFile(input_file).get_reference_length(keys[i])))
                read_modification(input_file, fetch, N_threads, name, directory, modification_level, header, region, method, context, modified_positions_data, random_sample, fastas, library, reverse_modification, outbam, line, modified_positions, modification_information)
                absolute_modification_information(modified_positions_data, modification_information, modified_positions,name,directory, modification_level, context, method, per_chromosome, line_)
            reads_info_output.close()
            reads_info_absolute.close()
            if modified_positions:
                csvfile.close()
            return modification_information
        elif per_chromosome == None and region != None:
            modified_positions_data = list()
            modification_information, random_sample, modification_level = modification_information_and_reads_fetching(context, user_defined_context, modified_positions, region, fastas, region[0], context_total_counts, modification_level, library, seed, input_file, N_threads, csvfile)
            fetch = tuple((region[0], region[1], region[2]))
            read_modification(input_file, fetch, N_threads, name, directory, modification_level, header, region, method, context, modified_positions_data, random_sample, fastas, library, reverse_modification, outbam, line, modified_positions, modification_information)
            absolute_modification_information(modified_positions_data, modification_information, modified_positions,name,directory, modification_level, context, method, per_chromosome, line_)
            reads_info_output.close()
            reads_info_absolute.close()
            if modified_positions:
                csvfile.close()
            return modification_information
        else:
            modified_positions_data = list()
            modification_information, random_sample, modification_level = modification_information_and_reads_fetching(context, user_defined_context, modified_positions, region, fastas, per_chromosome, context_total_counts, modification_level, library, seed, input_file, N_threads, csvfile)
            fetch = tuple((per_chromosome, 0, pysam.AlignmentFile(input_file).get_reference_length(per_chromosome)))
            read_modification(input_file, fetch, N_threads, name, directory, modification_level, header, region, method, context, modified_positions_data, random_sample, fastas, library, reverse_modification, outbam, line, modified_positions, modification_information)
            absolute_modification_information(modified_positions_data, modification_information, modified_positions,name,directory, modification_level, context, method, per_chromosome, line_)
            reads_info_output.close()
            reads_info_absolute.close()
            if modified_positions:
                csvfile.close()
            return modification_information



def modification_simulator(reference, read_length, input_file, method, library, simulation_input, modification_level,
                           modified_positions, coverage, context, region, directory, seed, user_defined_context, N_threads, per_chromosome, GC_bias, sequence_bias, overwrite, reverse_modification):
    "Assembles the whole modification simulator and runs per mode, method, library and context."
    time_s = datetime.now()
    logs.info("asTair's cytosine modification simulator started running. {} seconds".format((time_s - time_b).total_seconds()))
    header = True
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    if path.exists(directory) == False:
        raise Exception("The output directory does not exist.")
        sys.exit(1)
    if region.count(None)!=0:
        region = None
    if simulation_input == 'bam':
        if pysam.AlignmentFile(input_file).is_cram:
            extension = '.cram'
        else:
            extension = '.bam'
        try:
            modification_level = modification_level_transformation(modification_level, modified_positions)
            if not isinstance(modification_level, str):
                modification_level_ = int(modification_level*100)
            else:
                modification_level_ = modification_level
            modification_information = bam_input_simulation(directory, name, modification_level, context, input_file, reference, user_defined_context, per_chromosome, modified_positions, library, seed, region, None, method, N_threads, header, overwrite, extension, reverse_modification)
            if reverse_modification == False:
                if per_chromosome == None:
                    pysam.index(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + extension))
                else:
                    pysam.index(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + '_' + per_chromosome + extension))
            else:
                if per_chromosome == None:
                    pysam.index(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + '_reversed' + extension))
                else:
                    pysam.index(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + '_reversed_' + per_chromosome + extension))
        except AttributeError:
            logs.error(
                'The output files will not be overwritten. Please rename the input or the existing output files before rerunning if the input is different.',
                exc_info=True)
    time_m = datetime.now()
    logs.info("asTair's cytosine modification simulator finished running. {} seconds".format((
    time_m - time_b).total_seconds()))


if __name__ == '__main__':
    simulate()

