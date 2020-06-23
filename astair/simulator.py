#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import re
import os
import sys
import csv
import pdb
import zlib
import gzip
import mmap
import click
import pysam
import numpy
import random
import logging
import warnings
from os import path
from datetime import datetime
from collections import defaultdict


from astair.vcf_reader import read_vcf
from astair.cigar_search import cigar_search
from astair.safe_division import safe_rounder
from astair.bam_file_parser import bam_file_opener
from astair.context_search import context_sequence_search
from astair.cigar_search import position_correction_cigar
from astair.context_search import sequence_context_set_creation
from astair.simple_fasta_parser import fasta_splitting_by_sequence


@click.command()
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format used for generation and modification of the sequencing reads at desired contexts.')
@click.option('control_file', '--control_file', '-c', required=False, help='A VCF file with SNP status returned after WGS genotyping or a publicly available SNP list (dbSNP, 1000 genomes, etc).')
@click.option('model', '--model', '-mo', default='none',  type=click.Choice(['none']), required=False, help='Decide on model for class estimation.')
@click.option('read_length', '--read_length', '-l', type=int, required=True, help='Desired length of pair-end sequencing reads.')
@click.option('input_file', '--input_file', '-i', required=True, help='Sequencing reads as a BAM|CRAMfile or fasta sequence to generate reads.')
@click.option('simulation_input', '--simulation_input', '-si', type=click.Choice(['bam']), default='bam', required=False, help='Input file format according to the desired outcome. BAM|CRAM files can be generated with other WGS simulators allowing for sequencing errors and read distributions or can be real-life sequencing data.')
@click.option('method', '--method', '-m', required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT).')
@click.option('modification_level', '--modification_level', '-ml',  type=int, required=False, help='Desired modification level; can take any value between 0 and 100.')
@click.option('library', '--library', '-lb',  type=click.Choice(['directional', 'reverse']), default='directional', required=False, help='Provide the correct library construction method. NB: Non-directional methods under development.')
@click.option('modified_positions', '--modified_positions', '-mp', required=False, default=None, help='Provide a tab-delimited list of positions to be modified. By default the simulator randomly modifies certain positions. Please use seed for replication if no list is given.')
@click.option('context', '--context', '-co', required=False, default='all', type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be modified in the output file. Default behaviour is all, which modifies positions in CpG, CHG, CHH contexts. (Default all).')
@click.option('user_defined_context', '--user_defined_context', '-uc', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
@click.option('coverage', '--coverage', '-cv', required=False, type=int, help='Desired depth of sequencing coverage.')
@click.option('region', '--region', '-r', nargs=3, type=click.Tuple([str, int, int]), default=(None, None, None), required=False, help='The one-based genomic coordinates of the specific region of interest given in the form chromosome, start position, end position, e.g. chr1 100 2000.')
@click.option('overwrite', '--overwrite', '-ov', required=False, default=False, is_flag=True, help='Indicates whether existing output files with matching names will be overwritten. (Default False).')
@click.option('per_chromosome', '--per_chromosome', '-chr', type=str, help='When used, it calculates the modification rates only per the chromosome given.')
@click.option('GC_bias', '--GC_bias', '-gc', default=0.3, required=True, type=float, help='The value of total GC levels in the read above which lower coverage will be observed in Ns and fasta modes. (Default 0.5).')
@click.option('sequence_bias', '--sequence_bias', '-sb', default=0.1, required=True, type=float, help='The proportion of lower-case letters in the read string for the Ns and fasta modes that will decrease the chance of the read being output. (Default 0.1).')
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('reverse_modification', '--rev', '-rv', default=False, is_flag=True, required=False, help='Returns possible or known modified position to their unmodified expected state. NB: Works only on files with MD tags (Default False).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
@click.option('seed', '--seed', '-s', type=int, required=False, help='An integer number to be used as a seed for the random generators to ensure replication.')
@click.option('skip_clip_overlap', '--skip_clip_overlap', '-sc', required=False, default=False, type=bool, help='Skipping the random removal of overlapping bases between pair-end reads. Not recommended for pair-end libraries, unless the overlaps are removed prior to calling. (Default False)')
@click.option('single_end', '--single_end', '-se', default=False, is_flag=True, required=False, help='Indicates single-end sequencing reads (Default False).')
@click.option('minimum_base_quality', '--minimum_base_quality', '-bq', required=False, type=int, default=20, help='Set the minimum base quality for a read base to be used in the pileup (Default 20).')
@click.option('minimum_mapping_quality', '--minimum_mapping_quality', '-mq', required=False, type=int, default=0, help='Set the minimum mapping quality for a read to be used in the pileup (Default 0).')
@click.option('adjust_acapq_threshold', '--adjust_acapq_threshold', '-amq', required=False, type=int, default=0, help='Used to adjust the mapping quality with default 0 for no adjustment and a recommended value for adjustment 50. (Default 0).')
@click.option('add_indels', '--add_indels', '-ai', required=False, default=True, type=bool, help='Adds inserted bases and Ns for base skipped from the reference (Default True).')
@click.option('redo_baq', '--redo_baq', '-rbq', required=False, default=False, type=bool, help='Re-calculates per-Base Alignment Qualities ignoring existing base qualities (Default False).')
@click.option('compute_baq', '--compute_baq', '-cbq', required=False, default=True, type=bool, help='Performs re-alignment computing of per-Base Alignment Qualities (Default True).')
@click.option('ignore_orphans', '--ignore_orphans', '-io', required=False, default=True, type=bool, help='Ignore reads not in proper pairs (Default True).')
@click.option('max_depth', '--max_depth', '-md', required=False, type=int, default=250, help='Set the maximum read depth for the pileup. Please increase the maximum value for spike-ins and other highly-covered sequences. (Default 250).')
def simulate(reference,  control_file, model, read_length, input_file, simulation_input, method, modification_level, library, modified_positions, context, user_defined_context, coverage, region, overwrite, per_chromosome, GC_bias, sequence_bias, N_threads, reverse_modification, directory, seed, skip_clip_overlap, single_end, minimum_base_quality, minimum_mapping_quality, adjust_acapq_threshold, add_indels, redo_baq, compute_baq, ignore_orphans, max_depth):
    """Simulate TAPS/BS conversion on top of an existing bam/cram file."""
    modification_simulator(reference,  control_file, model, read_length, input_file, simulation_input, method, modification_level, library, modified_positions, context, user_defined_context, coverage, region, overwrite, per_chromosome, GC_bias, sequence_bias, N_threads, reverse_modification, directory, seed, skip_clip_overlap, single_end, minimum_base_quality, minimum_mapping_quality,adjust_acapq_threshold, add_indels, redo_baq, compute_baq, ignore_orphans, max_depth)


warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

#logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def csv_line_skipper(csvfile, start, key, tupler, input_file):
    """Random access for tab-delimited file reading."""
    cycles = 0
    csvfile.seek(start)
    inbam = bam_file_opener(input_file, None, 1)
    if isinstance(csvfile, gzip.GzipFile):
        for read in csvfile.readlines():
            cycles += 1
            if read.decode('utf8').split()[0] == key:
                if read.decode('utf8').split()[3].lower().islower()==False and float(
                        read.decode('utf8').split()[3]) != 0:
                    tupler[tuple((str(read.decode('utf8').split()[0]), int(read.decode('utf8').split()[1]), int(read.decode('utf8').split()[2])))] = numpy.array(object=[0, float(read.decode('utf8').split()[3]) * len(
                        [i.flag for i in inbam.fetch(contig=str(read.decode('utf8').split()[0]), start=int(read.decode('utf8').split()[1]), stop=int(read.decode('utf8').split()[2])) if
                         (i.flag in [99, 147] and read.decode('utf8').split()[7] in ['C', 'T']) or (
                         i.flag in [163, 83] and read.decode('utf8').split()[7] in ['A', 'G'])])], dtype=numpy.int8, copy=False)
            else:
                break
    else:
        for read in csvfile.readlines():
            cycles += 1
            if read.split()[0] == key:
                if read.split()[3].lower().islower()==False and float(read.split()[3]) != 0 and read.split()[0] == key:
                    tupler[tuple((str(read.split()[0]), int(read.split()[1]),
                                  int(read.split()[2])))] = numpy.array(
                        object=[0, float(read.split()[3]) * len(
                            [i.flag for i in inbam.fetch(contig=str(read.split()[0]),
                                                         start=int(read.split()[1]),
                                                         stop=int(read.split()[2])) if
                             (i.flag in [99, 147] and read.split()[7] in ['C', 'T']) or (
                                 i.flag in [163, 83] and read.split()[7] in ['A', 'G'])])], dtype=numpy.int8,
                        copy=False)
            else:
                break
    time_r = datetime.now()
    logs.info("Random access of positions is completed. {} seconds".format((
    time_r - time_b).total_seconds()))


def cytosine_modification_lookup(context, user_defined_context, modified_positions, region, fastas, keys, context_total_counts, input_file, N_threads, csvfile):
    """Finds all required cytosine contexts or takes positions from a tab-delimited file containing
     the list of positions to be modified."""
    if modified_positions is None:
        contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
        if region is None:
            modification_information = context_sequence_search(contexts, all_keys, fastas, keys, user_defined_context, context_total_counts, region, None)
        else:
            if region[0] in keys:
                 modification_information = context_sequence_search(contexts, all_keys, fastas, region[0], user_defined_context, context_total_counts, None)
        try:
            return modification_information
        except UnboundLocalError:
            logs.error('There is no reference sequence of this name in the provided fasta file.', exc_info=True)
            sys.exit(1)
    else:
        tupler = dict()
        try:
            if isinstance(csvfile, gzip.GzipFile):
                decompressed = zlib.decompress(csvfile, 16 + zlib.MAX_WBITS)
                memory_map = mmap.mmap(decompressed.fileno(), 0)
            else:
                memory_map = mmap.mmap(csvfile.fileno(), 0)
            start = memory_map.find(keys.encode('utf8'))
            csv_line_skipper(csvfile, start, keys, tupler, input_file)
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
    if modification_level is None:
        modification_level = 0
    if context == 'all' and modified_positions is None:
        modification_list_by_context = set()
        all_keys = list(('CHG','CHH','CpG'))
        for context_string in all_keys:
            modification_list_by_context = modification_list_by_context.union(set((keys) for keys, vals in modification_information.items() if vals[1] == context_string))
        required = safe_rounder(len(modification_list_by_context) * modification_level, 1, False)
    elif context != 'all' and modified_positions is None:
        modification_list_by_context = set((keys) for keys, vals in modification_information.items() if vals[1] == context)
        required = safe_rounder(len(modification_list_by_context) * modification_level, 1, False)
    else:
        random_sample = set(modification_information)
        modification_level = 'custom'
    if seed is not None and modified_positions is None:
        random.seed(seed)
        random_sample = set(random.sample(modification_list_by_context, int(required)))
    else:
        if modified_positions is None:
            random_sample = set(random.sample(modification_list_by_context, int(required)))
    return modification_level, random_sample


def general_read_information_output(read, header, line):
    """Writes to fastq file."""
    if read.is_read1 == True:
        orientation = '/1'
    else:
        orientation = '/2'
    try:
        if header:
            line.write(
                '{}\t{}\t{}\t{}\n'.format('#Read ID', 'reference', 'start', 'end'))
            line.write(
                '{}\t{}\t{}\t{}\n'.format(read.qname + orientation, read.reference_name, read.reference_start,
                                read.reference_start + read.query_length))
        else:
            line.write(
                '{}\t{}\t{}\t{}\n'.format(read.qname + orientation, read.reference_name, read.reference_start,
                                          read.reference_start + read.query_length))
    except IOError:
        logs.error('asTair cannot write to read information file.', exc_info=True)


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
        mod_level = safe_rounder((len(modified_positions_data) / context_list_length), 3, True)
    else:
        mod_level = 'Custom'
    try:
        line.write('__________________________________________________________________________________________________\n')
        line.write('Absolute modified positions: ' + str(len(modified_positions_data)) + '   |   ' +
                       'Percentage to all positions of the desired context: ' + str(mod_level) + ' %\n')
        line.write('__________________________________________________________________________________________________\n')
        for row in modified_positions_data:
            line.write('{}\t{}\t{}\n'.format(row[0], row[1], row[2]))
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
                                                                                                                 re.IGNORECASE) or re.findall('H',
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
                                                                                                                 re.IGNORECASE) or re.findall('H',
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


def bam_input_simulation(directory, name, modification_level, context, input_file, reference, user_defined_context, per_chromosome, modified_positions, library, seed, region, modified_positions_data, method, N_threads, header, overwrite, extension, reverse_modification):
    """Inserts modification information acording to method and context to a bam or cram file."""
    if not os.path.isfile(path.join(directory, name + '_' + method + '_' + str(modification_level) + '_' + context + extension)) or overwrite  == True:
        if pysam.AlignmentFile(input_file).is_cram:
            file_type = 'wc'
        else:
            file_type = 'wb'
        if not isinstance(modification_level, str):
            modification_level_= int(modification_level*100)
        else:
            modification_level_ = modification_level
        if modified_positions:
            if modified_positions[-3:] == '.gz':
                if sys.version[0] == '3':
                    csvfile = gzip.open(modified_positions, 'rt+', encoding='utf8', compresslevel=9)
                else:
                    csvfile = gzip.open(modified_positions, 'rt+', compresslevel=9)
            else:
                csvfile = open(modified_positions, 'r+')
        else:
            csvfile = None
        if reverse_modification == False:
            if per_chromosome is None:
                outbam = pysam.AlignmentFile(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + extension),
                file_type, reference_filename=reference, template=bam_file_opener(input_file, None, N_threads), header=header)
            else:
                outbam = pysam.AlignmentFile(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context  + '_' + per_chromosome  + extension),
                file_type, reference_filename=reference, template=bam_file_opener(input_file, None, N_threads), header=header)
        else:
            if per_chromosome is None:
                outbam = pysam.AlignmentFile(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + '_reversed' + extension),
                file_type, reference_filename=reference, template=bam_file_opener(input_file, None, N_threads), header=header)
            else:
                outbam = pysam.AlignmentFile(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context  + '_reversed_' + per_chromosome  + extension),
                file_type, reference_filename=reference, template=bam_file_opener(input_file, None, N_threads), header=header)
        keys, fastas = fasta_splitting_by_sequence(reference, per_chromosome, None, False, 'all')
        if per_chromosome is None:
            name_to_use = path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + '_read_information.txt')
            name_to_use_absolute = path.join(directory,name + '_' + method + '_' + str(modification_level_) + '_' + context + '_modified_positions_information.txt')
        else:
            name_to_use = path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + per_chromosome + '_' + context + '_read_information.txt')
            name_to_use_absolute = path.join(directory,name + '_' + method + '_' + str(modification_level_) + '_' + context + '_' + per_chromosome + '_modified_positions_information.txt')
        line = gzip.open(name_to_use + '.gz', 'wt', compresslevel=9)
        line_ = gzip.open(name_to_use_absolute + '.gz', 'wt', compresslevel=9)
        context_total_counts = defaultdict(int)
        if region is None:
            for i in range(0, len(keys)):
                modified_positions_data = list()
                modification_information, random_sample, modification_level = modification_information_and_reads_fetching(context, user_defined_context, modified_positions, region, fastas, keys[i], context_total_counts, modification_level, library, seed, input_file, N_threads, csvfile)
                fetch = tuple((keys[i], 0, pysam.AlignmentFile(input_file).get_reference_length(keys[i])))
                read_modification(input_file, fetch, N_threads, name, directory, modification_level, header, region, method, context, modified_positions_data, random_sample, fastas, library, reverse_modification, outbam, line, modified_positions, modification_information)
                absolute_modification_information(modified_positions_data, modification_information, modified_positions,name,directory, modification_level, context, method, per_chromosome, line_)
            line.close()
            line_.close()
            if modified_positions:
                csvfile.close()
        elif region is not None:
            modified_positions_data = list()
            modification_information, random_sample, modification_level = modification_information_and_reads_fetching(context, user_defined_context, modified_positions, region, fastas, region[0], context_total_counts, modification_level, library, seed, input_file, N_threads, csvfile)
            fetch = tuple((region[0], region[1], region[2]))
            read_modification(input_file, fetch, N_threads, name, directory, modification_level, header, region, method, context, modified_positions_data, random_sample, fastas, library, reverse_modification, outbam, line, modified_positions, modification_information)
            absolute_modification_information(modified_positions_data, modification_information, modified_positions,name,directory, modification_level, context, method, per_chromosome, line_)
            line.close()
            line_.close()
            if modified_positions:
                csvfile.close()



def modification_simulator(reference,  control_file, model, read_length, input_file, simulation_input, method, modification_level, library, modified_positions, context, user_defined_context, coverage, region, overwrite, per_chromosome, GC_bias, sequence_bias, N_threads, reverse_modification, directory, seed,
                           skip_clip_overlap, single_end, minimum_base_quality, minimum_mapping_quality, adjust_acapq_threshold, add_indels, redo_baq, compute_baq, ignore_orphans, max_depth):
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
    if model == 'none':
        model, labels, model_name = None, None, None
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
            bam_input_simulation(directory, name, modification_level, context, input_file, reference, user_defined_context, per_chromosome, modified_positions, library, seed, region, None, method, N_threads, header, overwrite, extension, reverse_modification)
            if reverse_modification == False:
                if per_chromosome is None:
                    pysam.index(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + extension))
                else:
                    pysam.index(path.join(directory, name + '_' + method + '_' + str(modification_level_) + '_' + context + '_' + per_chromosome + extension))
            else:
                if per_chromosome is None:
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

