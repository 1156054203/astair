#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import division

import re
import sys
import csv
import pdb
import numpy
import click
import logging
import warnings
from os import path
from math import ceil
from datetime import datetime

if sys.version[0] == '3':
    from itertools import zip_longest
elif sys.version[0] == '2':
    from itertools import izip_longest as zip_longest
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))

try:
    import matplotlib as mplot
    mplot.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    plt.style.use('seaborn-whitegrid')
    plt.ioff()
except Exception:
    warnings.warn("Matplotlib was not found, visualisation output will not be supported.", ImportWarning)



from astair.mbias import positions_discovery
from astair.safe_division import safe_rounder
from astair.safe_division import non_zero_division
from astair.bam_file_parser import bam_file_opener
from astair.DNA_sequences_operations import complementary
from astair.simple_fasta_parser import fasta_splitting_by_sequence


@click.command()
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format.')
@click.option('input_file', '--input_file', '-i', required=True, help='BAM|CRAM format file containing sequencing reads.')
@click.option('directory', '--directory', '-d', required=True, help='Output directory to save files.')
@click.option('read_length', '--read_length', '-l', type=int, required=True, help='The read length is needed to calculate the IDbias.')
@click.option('method', '--method', '-m',  required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT).')
@click.option('per_chromosome', '--per_chromosome', '-chr', type=str, help='When used, it calculates the modification rates only per the chromosome given.')
@click.option('single_end', '--single_end', '-se', default=False, is_flag=True, required=False, help='Indicates single-end sequencing reads (Default False).')
@click.option('plot', '--plot', '-p', required=False, is_flag=True, help='Phred scores will be visualised and output as a pdf file. Requires installed matplotlib.')
@click.option('colors', '--colors', '-c', default=['teal', 'deepskyblue', 'mediumblue', 'orange', 'gold', 'sienna'], type=list, required=False, help="List of color values used for visualistion of CpG, CHG and CHH modification levels per read, which are given as color1,color2,color3. Accepts valid matplotlib color names, RGB and RGBA hex strings and  single letters denoting color {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}. (Default 'teal', 'deepskyblue', 'mediumblue', 'orange', 'gold', 'sienna').")
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('no_information', '--no_information', '-ni', default='0', type=click.Choice(['.', '0', '*', 'NA']), required=False, help='What symbol should be used for a value where no enough quantative information is used. (Default 0).')
def idbias(reference, input_file, directory, read_length, method, single_end, plot, colors, N_threads, per_chromosome, no_information):
    """Generate indel count per read length information (IDbias). This is a quality-control measure."""
    IDbias_plotting(reference, input_file, directory, read_length, method, single_end, plot, colors, N_threads, per_chromosome, no_information)

warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

#logging.basicConfig(level=logging.WARNING)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def initialise_data_counters(read_length):
    """Clean initialisation of empty dictionaries used for counters."""
    all_read_data = list(({}, {}, {}, {}, {}, {}))
    for read_data in all_read_data:
        for i in range(0, read_length):
            read_data[i] = 0
    return all_read_data[0], all_read_data[1], all_read_data[2], all_read_data[3], all_read_data[4], all_read_data[5]


def idbias_calculator(flag, ref_sequence, read_sequence, read_length, read_mods_CHH_insert, read_mods_CHG_insert, read_mods_CpG_insert, read_mods_CHH_deletion, read_mods_CHG_deletion, read_mods_CpG_deletion, read_umod_CHH_insert, read_umod_CHG_insert, read_umod_CpG_insert, read_umod_CHH_deletion, read_umod_CHG_deletion, read_umod_CpG_deletion, method, single_end, inserts, deletions, read_insertions, read_deletions, cytosines, modified, unmodified, full_sequence):
    """Calculates the indel count per read position, pair orientation and cytosine context."""
    for i in range(0, len(read_sequence)):
        if len(inserts) >= 1 and i in inserts:
            read_insertions[i] += 1
        elif len(deletions) >= 1 and i in deletions:
            read_deletions[i] += 1
    if single_end == False:
        OT, OB = [99, 147], [83, 163]
    else:
        OT, OB = [0], [16]
    if flag in OT:
        cpg_all = [pos[0] for pos in cytosines if full_sequence[pos[1]:pos[1]+2].upper()=="CG"]
        chg_all = [pos[0] for pos in cytosines if full_sequence[pos[1]:pos[1]+3].upper() in ['CAG', 'CCG', 'CTG']]
        chh_all = [pos[0] for pos in cytosines if full_sequence[pos[1]:pos[1]+3].upper() in ['CAA', 'CAC', 'CAT', 'CCA', 'CCC', 'CCT', 'CTA', 'CTC', 'CTT']]
    else:
        cpg_all = [pos[0] for pos in cytosines if full_sequence[pos[1]-1:pos[1]+1].upper()=="CG"]
        contexts = [full_sequence[pos[1]-2:pos[1]+1].upper() for pos in cytosines]
        chg_all = [cytosines[ind][0] for ind in range(len(contexts)) if contexts[ind] in ['CAG', 'CGG', 'CTG']]
        chh_all = [cytosines[ind][0] for ind in range(len(contexts)) if contexts[ind] in ['AAG', 'AGG', 'ATG', 'GAG', 'GGG', 'GTG', 'TAG', 'TGG', 'TTG']]
    cpg_mods, chg_mods, chh_mods, cpg_umods, chg_umods, chh_umods = [x for x in modified if x in cpg_all], [x for x in modified if x in chg_all], [x for x in modified if x in chh_all], [x for x in unmodified if x in cpg_all], [x for x in unmodified if x in chg_all], [x for x in unmodified if x in chh_all]
    for i in range(0, len(read_sequence)):
        if len(inserts) >= 1 and i in inserts and any([k for k in range(i-5, i+5)]) in chh_mods:
            read_mods_CHH_insert[i] += 1
        elif len(inserts) >= 1 and i  in inserts and any([k for k in range(i-5, i+5)]) in chg_mods:
            read_mods_CHG_insert[i] += 1
        elif len(inserts) >= 1 and i  in inserts and any([k for k in range(i-5, i+5)]) in cpg_mods:
            read_mods_CpG_insert[i] += 1
        elif len(inserts) >= 1 and i  in inserts and any([k for k in range(i-5, i+5)]) in chh_umods:
            read_umod_CHH_insert[i] += 1
        elif len(inserts) >= 1 and i  in inserts and any([k for k in range(i-5, i+5)]) in chg_umods:
            read_umod_CHG_insert[i] += 1
        elif len(inserts) >= 1 and i  in inserts and any([k for k in range(i-5, i+5)]) in cpg_umods:
            read_umod_CpG_insert[i] += 1
        elif len(deletions) >= 1 and i in deletions and any([k for k in range(i-5, i+5)]) in chh_mods:
            read_mods_CHH_deletion[i] += 1
        elif len(deletions) >= 1 and i  in deletions and any([k for k in range(i-5, i+5)]) in chg_mods:
            read_mods_CHG_deletion[i] += 1
        elif len(deletions) >= 1 and i in deletions and any([k for k in range(i-5, i+5)]) in cpg_mods:
            read_mods_CpG_deletion[i] += 1
        elif len(deletions) >= 1 and i  in deletions and any([k for k in range(i-5, i+5)]) in chh_umods:
            read_umod_CHH_deletion[i] += 1
        elif len(deletions) >= 1 and i  in deletions and any([k for k in range(i-5, i+5)]) in chg_umods:
            read_umod_CHG_deletion[i] += 1
        elif len(deletions) >= 1 and i  in deletions and any([k for k in range(i-5, i+5)]) in cpg_umods:
            read_umod_CpG_deletion[i] += 1
    return read_mods_CHH_insert, read_mods_CHG_insert, read_mods_CpG_insert, read_mods_CHH_deletion, read_mods_CHG_deletion, read_mods_CpG_deletion, read_umod_CHH_insert, read_umod_CHG_insert, read_umod_CpG_insert, read_umod_CHH_deletion, read_umod_CHG_deletion, read_umod_CpG_deletion, read_insertions, read_deletions

def per_flag_counter(r_sequence, ref_name, r_start, read_length, exp1, exp2, flag, inserts,  method, single_end, deletions, total_read_count_R1, total_read_count_R2, insert_read_count_R1, insert_read_count_R2, deletion_read_count_R1, deletion_read_count_R2, indel_read_count_R1, indel_read_count_R2, read1_insertions, read1_deletions, read2_insertions, read2_deletions, read1_coverage, read2_coverage, read1_mods_CHH_insert, read1_mods_CHG_insert, read1_mods_CpG_insert, read1_mods_CHH_deletion, read1_mods_CHG_deletion, read1_mods_CpG_deletion, read1_umod_CHH_insert, read1_umod_CHG_insert, read1_umod_CpG_insert, read1_umod_CHH_deletion, read1_umod_CHG_deletion, read1_umod_CpG_deletion, read2_mods_CHH_insert, read2_mods_CHG_insert, read2_mods_CpG_insert, read2_mods_CHH_deletion, read2_mods_CHG_deletion, read2_mods_CpG_deletion, read2_umod_CHH_insert, read2_umod_CHG_insert, read2_umod_CpG_insert, read2_umod_CHH_deletion, read2_umod_CHG_deletion, read2_umod_CpG_deletion, cytosines, modified, unmodified, fastas):
    """Does a correct indel and mod calculation per read flag."""
    if flag in exp1:
        total_read_count_R1 += 1
        if len(inserts) >= 1 and len(deletions) == 0:
            insert_read_count_R1 += 1
        if len(deletions)>= 1 and len(inserts) == 0:
            deletion_read_count_R1 += 1
        if len(inserts)>= 1 and len(deletions) >= 1:
            indel_read_count_R1 += 1
        for i in range(0, len(r_sequence)):
            read1_coverage[i] += 1
        idbias_calculator(flag, fastas[ref_name][r_start:r_start+len(r_sequence)], r_sequence, read_length, read1_mods_CHH_insert, read1_mods_CHG_insert, read1_mods_CpG_insert, read1_mods_CHH_deletion, read1_mods_CHG_deletion, read1_mods_CpG_deletion, read1_umod_CHH_insert, read1_umod_CHG_insert, read1_umod_CpG_insert, read1_umod_CHH_deletion, read1_umod_CHG_deletion, read1_umod_CpG_deletion, method, single_end, inserts, deletions, read1_insertions, read1_deletions, cytosines, modified, unmodified, fastas[ref_name])
    else:
        total_read_count_R2 += 1
        if len(inserts) >= 1 and len(deletions) == 0:
            insert_read_count_R2 += 1
        if len(deletions)>= 1 and len(inserts) == 0:
            deletion_read_count_R2 += 1
        if len(inserts)>= 1 and len(deletions) >= 1:
            indel_read_count_R2 += 1
        for i in range(0, len(r_sequence)):
            read2_coverage[i] += 1
        idbias_calculator(flag, fastas[ref_name][r_start:r_start+len(r_sequence)], r_sequence, read_length, read2_mods_CHH_insert, read2_mods_CHG_insert, read2_mods_CpG_insert, read2_mods_CHH_deletion, read2_mods_CHG_deletion, read2_mods_CpG_deletion, read2_umod_CHH_insert, read2_umod_CHG_insert, read2_umod_CpG_insert, read2_umod_CHH_deletion, read2_umod_CHG_deletion, read2_umod_CpG_deletion, method, single_end, inserts, deletions, read2_insertions, read2_deletions, cytosines, modified, unmodified, fastas[ref_name])
    return total_read_count_R1, total_read_count_R2, insert_read_count_R1, insert_read_count_R2, deletion_read_count_R1, deletion_read_count_R2, indel_read_count_R1, indel_read_count_R2, read1_insertions, read1_deletions, read2_insertions, read2_deletions, read1_coverage, read2_coverage, read1_mods_CHH_insert, read1_mods_CHG_insert, read1_mods_CpG_insert, read1_mods_CHH_deletion, read1_mods_CHG_deletion, read1_mods_CpG_deletion, read1_umod_CHH_insert, read1_umod_CHG_insert, read1_umod_CpG_insert, read1_umod_CHH_deletion, read1_umod_CHG_deletion, read1_umod_CpG_deletion, read2_mods_CHH_insert, read2_mods_CHG_insert, read2_mods_CpG_insert, read2_mods_CHH_deletion, read2_mods_CHG_deletion, read2_mods_CpG_deletion, read2_umod_CHH_insert, read2_umod_CHG_insert, read2_umod_CpG_insert, read2_umod_CHH_deletion, read2_umod_CHG_deletion, read2_umod_CpG_deletion
            


def idbias_evaluater(input_file, read_length, method, single_end, N_threads, fastas, per_chromosome):
    """Outputs the modification levels per read position, pair orientation and cytosine context."""
    read1_mods_CHH_insert, read1_mods_CHG_insert, read1_mods_CpG_insert, read1_mods_CHH_deletion, read1_mods_CHG_deletion, read1_mods_CpG_deletion = initialise_data_counters(read_length)
    read1_umod_CHH_insert, read1_umod_CHG_insert, read1_umod_CpG_insert, read1_umod_CHH_deletion, read1_umod_CHG_deletion, read1_umod_CpG_deletion = initialise_data_counters(read_length)
    read2_mods_CHH_insert, read2_mods_CHG_insert, read2_mods_CpG_insert, read2_mods_CHH_deletion, read2_mods_CHG_deletion, read2_mods_CpG_deletion = initialise_data_counters(read_length)
    read2_umod_CHH_insert, read2_umod_CHG_insert, read2_umod_CpG_insert, read2_umod_CHH_deletion, read2_umod_CHG_deletion, read2_umod_CpG_deletion = initialise_data_counters(read_length)
    total_read_count_R1, total_read_count_R2, insert_read_count_R1, insert_read_count_R2, deletion_read_count_R1, deletion_read_count_R2, indel_read_count_R1, indel_read_count_R2, read1_insertions, read1_deletions, read2_insertions, read2_deletions, read1_coverage, read2_coverage = 0, 0, 0, 0, 0, 0, 0, 0, {}, {}, {}, {}, {}, {}
    for i in range(0, read_length):
        read1_insertions[i], read1_deletions[i], read1_coverage[i], read2_insertions[i], read2_deletions[i], read2_coverage[i] = 0, 0, 0, 0, 0, 0
    for read in bam_file_opener(input_file, 'fetch', N_threads):
        if len(re.findall("(?:.*D.*)",str(read.cigarstring))) > 0 or len(re.findall("(?:.*I.*)",str(read.cigarstring))) > 0:
            inserts = [i[0] for i in read.get_aligned_pairs(with_seq=True)  if i[1] is None]
            deletions = [i[1]-read.reference_start for i in read.get_aligned_pairs(with_seq=True) if i[0] is None]
        else:
             inserts, deletions = [], []
        if single_end == False:
            exp1, exp2 = [99, 83], [147, 163]
        else:
            exp1, exp2 = [0], [16]
        if read.reference_length != 0 and read.flag in exp1 + exp2 and ((per_chromosome is None) or (per_chromosome is not None and read.reference_name == per_chromosome)):
            flag, ref_name, r_start, r_sequence = read.flag, read.reference_name, read.reference_start, read.query_sequence
            if single_end == False:
                if flag == 99 or flag == 147:
                    ref, alt = 'C', 'T'
                elif flag == 163 or flag == 83:
                    ref, alt, = 'G', 'A'
            else:
                if flag == 0:
                    ref, alt = 'C', 'T'
                elif flag == 16:
                    ref, alt, = 'G', 'A'
            cytosines, modified, unmodified = positions_discovery(read, fastas, ref_name, ref, alt, r_sequence, method)
            total_read_count_R1, total_read_count_R2, insert_read_count_R1, insert_read_count_R2, deletion_read_count_R1, deletion_read_count_R2, indel_read_count_R1, indel_read_count_R2, read1_insertions, read1_deletions, read2_insertions, read2_deletions, read1_coverage, read2_coverage, read1_mods_CHH_insert, read1_mods_CHG_insert, read1_mods_CpG_insert, read1_mods_CHH_deletion, read1_mods_CHG_deletion, read1_mods_CpG_deletion, read1_umod_CHH_insert, read1_umod_CHG_insert, read1_umod_CpG_insert, read1_umod_CHH_deletion, read1_umod_CHG_deletion, read1_umod_CpG_deletion, read2_mods_CHH_insert, read2_mods_CHG_insert, read2_mods_CpG_insert, read2_mods_CHH_deletion, read2_mods_CHG_deletion, read2_mods_CpG_deletion, read2_umod_CHH_insert, read2_umod_CHG_insert, read2_umod_CpG_insert, read2_umod_CHH_deletion, read2_umod_CHG_deletion, read2_umod_CpG_deletion = per_flag_counter(r_sequence, ref_name, r_start, read_length, exp1, exp2, flag, inserts,  method, single_end, deletions, total_read_count_R1, total_read_count_R2, insert_read_count_R1, insert_read_count_R2, deletion_read_count_R1, deletion_read_count_R2, indel_read_count_R1, indel_read_count_R2, read1_insertions, read1_deletions, read2_insertions, read2_deletions, read1_coverage, read2_coverage, read1_mods_CHH_insert, read1_mods_CHG_insert, read1_mods_CpG_insert, read1_mods_CHH_deletion, read1_mods_CHG_deletion, read1_mods_CpG_deletion, read1_umod_CHH_insert, read1_umod_CHG_insert, read1_umod_CpG_insert, read1_umod_CHH_deletion, read1_umod_CHG_deletion, read1_umod_CpG_deletion, read2_mods_CHH_insert, read2_mods_CHG_insert, read2_mods_CpG_insert, read2_mods_CHH_deletion, read2_mods_CHG_deletion, read2_mods_CpG_deletion, read2_umod_CHH_insert, read2_umod_CHG_insert, read2_umod_CpG_insert, read2_umod_CHH_deletion, read2_umod_CHG_deletion, read2_umod_CpG_deletion, cytosines, modified, unmodified, fastas)
    return read1_mods_CHH_insert, read1_mods_CHG_insert, read1_mods_CpG_insert, read1_mods_CHH_deletion, read1_mods_CHG_deletion, read1_mods_CpG_deletion, read1_umod_CHH_insert, read1_umod_CHG_insert, read1_umod_CpG_insert, read1_umod_CHH_deletion, read1_umod_CHG_deletion, read1_umod_CpG_deletion, read1_insertions, read1_deletions, read2_insertions, read2_deletions, read2_mods_CHH_insert, read2_mods_CHG_insert, read2_mods_CpG_insert, read2_mods_CHH_deletion, read2_mods_CHG_deletion, read2_mods_CpG_deletion, read2_umod_CHH_insert, read2_umod_CHG_insert, read2_umod_CpG_insert, read2_umod_CHH_deletion, read2_umod_CHG_deletion, read2_umod_CpG_deletion, total_read_count_R1, total_read_count_R2, insert_read_count_R1, insert_read_count_R2, deletion_read_count_R1, deletion_read_count_R2, indel_read_count_R1, indel_read_count_R2, read1_coverage, read2_coverage



def indel_calculator(i, read_mods, read_umods, read_values, no_information):
    """Calculates summary statistics per context and read orientation."""
    read_values[i] = non_zero_division(read_mods[i], read_mods[i] + read_umods[i], no_information)
    values = [(keys + 1, safe_rounder(vals[0], 12, True)) if isinstance(vals, list) else (keys + 1, safe_rounder(vals, 12, True)) for keys, vals in read_values.items()]
    mod_counts = sum(list(read_mods.values()))
    umod_counts = sum(list(read_umods.values())) + mod_counts
    return read_values, values, umod_counts, mod_counts


def idbias_statistics_calculator(fastas, input_file, name, directory, read_length, method, single_end, N_threads, per_chromosome, no_information):
    """Creates a summary statistics of the modification levels per read position, pair orientation and cytosine context,
    and then writes them as a text file that can be used for independent visualisation."""
    read1_mods_CHH_insert, read1_mods_CHG_insert, read1_mods_CpG_insert, read1_mods_CHH_deletion, read1_mods_CHG_deletion, read1_mods_CpG_deletion, read1_umod_CHH_insert, read1_umod_CHG_insert, read1_umod_CpG_insert, read1_umod_CHH_deletion, read1_umod_CHG_deletion, read1_umod_CpG_deletion, read1_insertions, read1_deletions, read2_insertions, read2_deletions, read2_mods_CHH_insert, read2_mods_CHG_insert, read2_mods_CpG_insert, read2_mods_CHH_deletion, read2_mods_CHG_deletion, read2_mods_CpG_deletion, read2_umod_CHH_insert, read2_umod_CHG_insert, read2_umod_CpG_insert, read2_umod_CHH_deletion, read2_umod_CHG_deletion, read2_umod_CpG_deletion, total_read_count_R1, total_read_count_R2, insert_read_count_R1, insert_read_count_R2, deletion_read_count_R1, deletion_read_count_R2, indel_read_count_R1, indel_read_count_R2, read1_coverage, read2_coverage = idbias_evaluater(input_file, read_length, method, single_end, N_threads, fastas, per_chromosome)
    read_values_1_CpG_insert, read_values_1_CHG_insert, read_values_1_CHH_insert, read_values_1_CpG_deletion, read_values_1_CHG_deletion, read_values_1_CHH_deletion = initialise_data_counters(read_length)
    read_values_2_CpG_insert, read_values_2_CHG_insert, read_values_2_CHH_insert, read_values_2_CpG_deletion, read_values_2_CHG_deletion, read_values_2_CHH_deletion  = initialise_data_counters(read_length)
    if single_end == True:
        total_read_count_R1 += total_read_count_R2
        insert_read_count_R1 += insert_read_count_R2
        deletion_read_count_R1 += deletion_read_count_R2
        indel_read_count_R1 += indel_read_count_R2
        for i in range(0, read_length):
            read1_mods_CpG_insert[i] += read2_mods_CpG_insert[i]
            read1_mods_CpG_deletion[i] += read2_mods_CpG_deletion[i]
            read1_mods_CHG_insert[i] += read2_mods_CHG_insert[i]
            read1_mods_CHG_deletion[i] += read2_mods_CHG_deletion[i]
            read1_mods_CHH_insert[i] += read2_mods_CHH_insert[i]
            read1_mods_CHH_deletion[i] += read2_mods_CHH_deletion[i]
            read1_umod_CpG_insert[i] += read2_umod_CpG_insert[i]
            read1_umod_CpG_deletion[i] += read2_umod_CpG_deletion[i]
            read1_umod_CHG_insert[i] += read2_umod_CHG_insert[i]
            read1_umod_CHG_deletion[i] += read2_umod_CHG_deletion[i]
            read1_umod_CHH_insert[i] += read2_umod_CHH_insert[i]
            read1_umod_CHH_deletion[i] += read2_umod_CHH_deletion[i]
            read1_insertions[i] += read2_insertions[i]
            read1_deletions[i] += read2_deletions[i] 
            read1_coverage[i] += read2_coverage[i]
    for i in range(0, read_length):
        read_values_1_CHH_insert, values_1_CHH_insert, umod_counts_1_CHH_insert, mod_counts_1_CHH_insert = indel_calculator(i, read1_mods_CHH_insert, read1_umod_CHH_insert, read_values_1_CHH_insert, no_information)
        read_values_1_CHH_deletion, values_1_CHH_deletion, umod_counts_1_CHH_deletion, mod_counts_1_CHH_deletion = indel_calculator(i, read1_mods_CHH_deletion, read1_umod_CHH_deletion, read_values_1_CHH_deletion, no_information)
        read_values_1_CHG_insert, values_1_CHG_insert, umod_counts_1_CHG_insert, mod_counts_1_CHG_insert = indel_calculator(i, read1_mods_CHG_insert, read1_umod_CHG_insert, read_values_1_CHG_insert, no_information)
        read_values_1_CHG_deletion, values_1_CHG_deletion, umod_counts_1_CHG_deletion, mod_counts_1_CHG_deletion = indel_calculator(i, read1_mods_CHG_deletion, read1_umod_CHG_deletion, read_values_1_CHG_deletion, no_information)
        read_values_1_CpG_insert, values_1_CpG_insert, umod_counts_1_CpG_insert, mod_counts_1_CpG_insert = indel_calculator(i, read1_mods_CpG_insert, read1_umod_CpG_insert, read_values_1_CpG_insert, no_information)
        read_values_1_CpG_deletion, values_1_CpG_deletion, umod_counts_1_CpG_deletion, mod_counts_1_CpG_deletion = indel_calculator(i, read1_mods_CpG_deletion, read1_umod_CpG_deletion, read_values_1_CpG_deletion, no_information)
        read_values_2_CHH_insert, values_2_CHH_insert, umod_counts_2_CHH_insert, mod_counts_2_CHH_insert = indel_calculator(i, read2_mods_CHH_insert, read2_umod_CHG_insert, read_values_2_CHH_insert, no_information)
        read_values_2_CHH_deletion, values_2_CHH_deletion, umod_counts_2_CHH_deletion, mod_counts_2_CHH_deletion = indel_calculator(i, read2_mods_CHH_deletion, read2_umod_CHH_deletion, read_values_2_CHH_deletion, no_information)
        read_values_2_CHG_insert, values_2_CHG_insert, umod_counts_2_CHG_insert, mod_counts_2_CHG_insert = indel_calculator(i, read2_mods_CHG_insert, read2_umod_CHG_insert, read_values_2_CHG_insert, no_information)
        read_values_2_CHG_deletion, values_2_CHG_deletion, umod_counts_2_CHG_deletion, mod_counts_2_CHG_deletion = indel_calculator(i, read2_mods_CHG_deletion, read2_umod_CHG_deletion, read_values_2_CHG_deletion, no_information)
        read_values_2_CpG_insert, values_2_CpG_insert, umod_counts_2_CpG_insert, mod_counts_2_CpG_insert = indel_calculator(i, read2_mods_CpG_insert, read2_umod_CpG_insert, read_values_2_CpG_insert, no_information)
        read_values_2_CpG_deletion, values_2_CpG_deletion, umod_counts_2_CpG_deletion, mod_counts_2_CpG_deletion = indel_calculator(i, read2_mods_CpG_deletion, read2_umod_CpG_deletion, read_values_2_CpG_deletion, no_information)
        if per_chromosome is not None:
            calls_output = open(directory + name + '_' + per_chromosome +  '_ID-bias.stats', 'w')
        else:
            calls_output = open(directory + name + '_ID-bias.stats', 'w')
        data_line = csv.writer(calls_output, delimiter='\t', lineterminator='\n')
        data_line.writerow(["#POS","R1_total_reads", "R1_total_insertions", "R1_total_deletions", "R1_total_both", "R1_CHH_insert_total", "R1_CHG_insert_total", "R1_CpG_insert_total", "R1_CHH_insert_modified", "R1_CHG_insert_modified", "R1_CpG_insert_modified", "R1_CHH_deletion_total", "R1_CHG_deletion_total", "R1_CpG_deletion_total", "R1_CHH_deletion_modified", "R1_CHG_deletion_modified", "R1_CpG_deletion_modified", "R2_total_reads", "R2_total_insertions", "R2_total_deletions", "R2_total_both", "R2_CHH_insert_total", "R2_CHG_insert_total", "R2_CpG_insert_total", "R2_CHH_insert_modified", "R2_CHG_insert_modified", "R2_CpG_insert_modified", "R2_CHH_deletion_total", "R2_CHG_deletion_total", "R2_CpG_deletion_total", "R2_CHH_deletion_modified", "R2_CHG_deletion_modified", "R2_CpG_deletion_modified"])
        data_line.writerow([no_information, total_read_count_R1, insert_read_count_R1, deletion_read_count_R1, indel_read_count_R1, umod_counts_1_CHH_insert, umod_counts_1_CHG_insert, umod_counts_1_CpG_insert, mod_counts_1_CHH_insert, mod_counts_1_CHG_insert, mod_counts_1_CpG_insert, umod_counts_1_CHH_deletion, umod_counts_1_CHG_deletion, umod_counts_1_CpG_deletion, mod_counts_1_CHH_deletion, mod_counts_1_CHG_deletion, mod_counts_1_CpG_deletion, total_read_count_R2, insert_read_count_R2, deletion_read_count_R2, indel_read_count_R2, umod_counts_2_CHH_insert, umod_counts_2_CHG_insert, umod_counts_2_CpG_insert, mod_counts_2_CHH_insert, mod_counts_2_CHG_insert, mod_counts_2_CpG_insert, umod_counts_2_CHH_deletion, umod_counts_2_CHG_deletion, umod_counts_2_CpG_deletion, mod_counts_2_CHH_deletion, mod_counts_2_CHG_deletion, mod_counts_2_CpG_deletion])
        for i in range(read_length):
            data_line.writerow([i+1, read1_coverage[i], read1_insertions[i], read1_deletions[i], no_information, list(read1_umod_CHH_insert.values())[i], list(read1_umod_CHG_insert.values())[i], list(read1_umod_CpG_insert.values())[i], list(read1_mods_CHH_insert.values())[i], list(read1_mods_CHG_insert.values())[i], list(read1_mods_CpG_insert.values())[i], list(read1_umod_CHH_deletion.values())[i], list(read1_umod_CHG_deletion.values())[i], list(read1_umod_CpG_deletion.values())[i], list(read1_mods_CHH_deletion.values())[i], list(read1_mods_CHG_deletion.values())[i], list(read1_mods_CpG_deletion.values())[i], read2_coverage[i], read2_insertions[i], read2_deletions[i], no_information, list(read2_umod_CHH_insert.values())[i], list(read2_umod_CHG_insert.values())[i], list(read2_umod_CpG_insert.values())[i], list(read2_mods_CHH_insert.values())[i], list(read2_mods_CHG_insert.values())[i], list(read2_mods_CpG_insert.values())[i], list(read2_umod_CHH_deletion.values())[i], list(read2_umod_CHG_deletion.values())[i], list(read2_umod_CpG_deletion.values())[i], list(read2_mods_CHH_deletion.values())[i], list(read2_mods_CHG_deletion.values())[i], list(read2_mods_CpG_deletion.values())[i]])
    return read_values_1_CHH_insert, values_1_CHH_insert, umod_counts_1_CHH_insert, mod_counts_1_CHH_insert, read_values_1_CHH_deletion, values_1_CHH_deletion, umod_counts_1_CHH_deletion, mod_counts_1_CHH_deletion, read_values_1_CHG_insert, values_1_CHG_insert, umod_counts_1_CHG_insert, mod_counts_1_CHG_insert, read_values_1_CHG_deletion, values_1_CHG_deletion, umod_counts_1_CHG_deletion, mod_counts_1_CHG_deletion, read_values_1_CpG_insert, values_1_CpG_insert, umod_counts_1_CpG_insert, mod_counts_1_CpG_insert, read_values_1_CpG_deletion, values_1_CpG_deletion, umod_counts_1_CpG_deletion, mod_counts_1_CpG_deletion, read_values_2_CHH_insert, values_2_CHH_insert, umod_counts_2_CHH_insert, mod_counts_2_CHH_insert, read_values_2_CHH_deletion, values_2_CHH_deletion, umod_counts_2_CHH_deletion, mod_counts_2_CHH_deletion, read_values_2_CHG_insert, values_2_CHG_insert, umod_counts_2_CHG_insert, mod_counts_2_CHG_insert, read_values_2_CHG_deletion, values_2_CHG_deletion, umod_counts_2_CHG_deletion, mod_counts_2_CHG_deletion, read_values_2_CpG_insert, values_2_CpG_insert, umod_counts_2_CpG_insert, mod_counts_2_CpG_insert, read_values_2_CpG_deletion, values_2_CpG_deletion, umod_counts_2_CpG_deletion, mod_counts_2_CpG_deletion, total_read_count_R1, total_read_count_R2, insert_read_count_R1, insert_read_count_R2, deletion_read_count_R1, deletion_read_count_R2, indel_read_count_R1, indel_read_count_R2, read1_coverage, read2_coverage, read1_insertions, read1_deletions, read2_insertions, read2_deletions 

def IDbias_plotting(reference, input_file, directory, read_length, method, single_end, plot, colors, N_threads, per_chromosome, no_information):
    """The general ID-bias calculation and statistics output function, which might be also visualised if the plotting module is enabled."""
    time_s = datetime.now()
    logs.info("asTair's ID-bias summary function started running. {} seconds".format((time_s - time_b).total_seconds()))
    name = path.splitext(path.basename(input_file))[0]
    if list(directory)[-1]!="/":
        directory = directory + "/"
    if path.exists(directory) == False:
        raise Exception("The output directory does not exist.")
        sys.exit(1)
    keys, fastas = fasta_splitting_by_sequence(reference, per_chromosome, None, False, 'all')
    if no_information == '0':
        no_information = int(no_information)
    read_values_1_CHH_insert, values_1_CHH_insert, umod_counts_1_CHH_insert, mod_counts_1_CHH_insert, read_values_1_CHH_deletion, values_1_CHH_deletion, umod_counts_1_CHH_deletion, mod_counts_1_CHH_deletion, read_values_1_CHG_insert, values_1_CHG_insert, umod_counts_1_CHG_insert, mod_counts_1_CHG_insert, read_values_1_CHG_deletion, values_1_CHG_deletion, umod_counts_1_CHG_deletion, mod_counts_1_CHG_deletion, read_values_1_CpG_insert, values_1_CpG_insert, umod_counts_1_CpG_insert, mod_counts_1_CpG_insert, read_values_1_CpG_deletion, values_1_CpG_deletion, umod_counts_1_CpG_deletion, mod_counts_1_CpG_deletion, read_values_2_CHH_insert, values_2_CHH_insert, umod_counts_2_CHH_insert, mod_counts_2_CHH_insert, read_values_2_CHH_deletion, values_2_CHH_deletion, umod_counts_2_CHH_deletion, mod_counts_2_CHH_deletion, read_values_2_CHG_insert, values_2_CHG_insert, umod_counts_2_CHG_insert, mod_counts_2_CHG_insert, read_values_2_CHG_deletion, values_2_CHG_deletion, umod_counts_2_CHG_deletion, mod_counts_2_CHG_deletion, read_values_2_CpG_insert, values_2_CpG_insert, umod_counts_2_CpG_insert, mod_counts_2_CpG_insert, read_values_2_CpG_deletion, values_2_CpG_deletion, umod_counts_2_CpG_deletion, mod_counts_2_CpG_deletion, total_read_count_R1, total_read_count_R2, insert_read_count_R1, insert_read_count_R2, deletion_read_count_R1, deletion_read_count_R2, indel_read_count_R1, indel_read_count_R2, read1_coverage, read2_coverage, read1_insertions, read1_deletions, read2_insertions, read2_deletions  = idbias_statistics_calculator(fastas, input_file, name, directory, read_length, method, single_end, N_threads, per_chromosome, no_information)
    try:
        read1_normalised_insert_values, read1_normalised_deletion_values, read2_normalised_insert_values, read2_normalised_deletion_values = {}, {}, {}, {}
        for i in range(read_length):
            read1_normalised_insert_values[i] = non_zero_division(read1_insertions[i], read1_coverage[i], no_information)*100
            read1_normalised_deletion_values[i] = non_zero_division(read1_deletions[i], read1_coverage[i], no_information)*100
            read2_normalised_insert_values[i] = non_zero_division(read2_insertions[i], read2_coverage[i], no_information)*100
            read2_normalised_deletion_values[i] = non_zero_division(read2_deletions[i], read2_coverage[i], no_information)*100
        read1_normalised_insert_values = list(read1_normalised_insert_values.values())
        read1_normalised_deletion_values = list(read1_normalised_deletion_values.values())
        read2_normalised_insert_values = list(read2_normalised_insert_values.values())
        read2_normalised_deletion_values = list(read2_normalised_deletion_values.values())
        if plot:
            if colors != ['teal', 'deepskyblue', 'mediumblue', 'orange', 'gold', 'sienna']:
                colors = "".join(colors).split(',') 
            plt.figure(figsize=(10, 20))
            if single_end == False:
                fig, fq = plt.subplots(1, 2, sharey=True)
            else:
                fig, fq = plt.subplots(1, 1)
            fig.suptitle('Sequencing ID-bias: relative indel abundance', fontsize=14)
            if single_end == False:
                plt.subplots_adjust(wspace=0.4)
                plt.subplots_adjust(right=1)
                fq[0].set_ylabel('Relative abundance, %', fontsize=12)
                fq[0].set_xlabel('First in pair indel type', fontsize=12)
                fq[0].bar([0,1,2], [non_zero_division(insert_read_count_R1, total_read_count_R1, 0)*100, non_zero_division(deletion_read_count_R1, total_read_count_R1, 0)*100, non_zero_division(indel_read_count_R1, total_read_count_R1, 0)*100], color=['lightgray', 'deepskyblue', 'mediumblue'])
                fq[0].set_xticklabels(['',  'insert', 'deletion', 'both'],  fontsize=12)
                fq[0].grid(color='lightgray', linestyle='solid', linewidth=1)
                fq[1].set_ylabel('Relative abundance, %', fontsize=12)
                fq[1].set_xlabel('Second in pair indel type', fontsize=12)
                fq[1].bar([0,1,2], [non_zero_division(insert_read_count_R2, total_read_count_R2, 0)*100, non_zero_division(deletion_read_count_R2, total_read_count_R2, 0)*100, non_zero_division(indel_read_count_R2, total_read_count_R2, 0)*100], color=['lightgray', 'deepskyblue', 'mediumblue'])
                fq[1].set_xticklabels(['', 'insert', 'deletion', 'both'],  fontsize=12)
                fq[1].grid(color='lightgray', linestyle='solid', linewidth=1)
            else:
                fq.set_ylabel('Relative abundance, %', fontsize=12)
                fq.set_xlabel('Indel type', fontsize=12)
                fq.bar([0,1,2], [non_zero_division(insert_read_count_R1, total_read_count_R1, 0)*100, non_zero_division(deletion_read_count_R1, total_read_count_R1, 0)*100, non_zero_division(indel_read_count_R1, total_read_count_R1, 0)*100], color=['lightgray', 'deepskyblue', 'mediumblue'])
                fq.set_xticklabels(['insert', 'deletion', 'both'],  fontsize=12)
                fq.grid(color='lightgray', linestyle='solid', linewidth=1)
            plt.savefig(directory + name + '_ID-bias_abundance_plot.pdf', figsize=(10, 20), dpi=330, bbox_inches='tight', pad_inches=0.25)
            plt.close()
              
            plt.figure(figsize=(16, 20))
            if single_end == False:
                fig, fq = plt.subplots(1, 2,  sharey=True)
            else:
                fig, fq = plt.subplots(1, 1)
            fig.suptitle('Sequencing ID-bias: relative indel abundance in 10bp from a modified cytosine', fontsize=14)
            if single_end == False:
                plt.subplots_adjust(wspace=0.4)
                plt.subplots_adjust(right=1)
                fq[0].set_ylabel('Relative abundance, %', fontsize=12)
                fq[0].set_xlabel('First in pair', fontsize=12)
                fq[0].bar([0,1,2,3,4,5], [non_zero_division(mod_counts_1_CpG_insert, umod_counts_1_CpG_insert, 0)*100, non_zero_division(mod_counts_1_CHG_insert, umod_counts_1_CHG_insert, 0)*100, non_zero_division(mod_counts_1_CHH_insert, umod_counts_1_CHH_insert, 0)*100, non_zero_division(mod_counts_1_CpG_deletion, umod_counts_1_CpG_deletion, 0)*100, non_zero_division(mod_counts_1_CHG_deletion, umod_counts_1_CHG_deletion, 0)*100, non_zero_division(mod_counts_1_CHH_deletion, umod_counts_1_CHH_deletion, 0)*100], color=colors)
                fq[0].set_xticks([0,1,2,3,4,5])
                fq[0].set_xticklabels(['CpG insert', 'CHG insert', 'CHH insert', 'CpG deletion', 'CHG deletion', 'CHH deletion'],  fontsize=12, rotation=90)
                fq[0].grid(color='lightgray', linestyle='solid', linewidth=1)
                fq[1].set_ylabel('Relative abundance, %', fontsize=12)
                fq[1].set_xlabel('Second in pair', fontsize=12)
                print(mod_counts_2_CHG_insert, umod_counts_2_CHG_insert)
                fq[1].bar([0,1,2,3,4,5], [non_zero_division(mod_counts_2_CpG_insert, umod_counts_2_CpG_insert, 0)*100, non_zero_division(mod_counts_2_CHG_insert, umod_counts_2_CHG_insert, 0)*100, non_zero_division(mod_counts_2_CHH_insert, umod_counts_2_CHH_insert, 0)*100, non_zero_division(mod_counts_2_CpG_deletion, umod_counts_2_CpG_deletion, 0)*100, non_zero_division(mod_counts_2_CHG_deletion, umod_counts_2_CHG_deletion, 0)*100, non_zero_division(mod_counts_2_CHH_deletion, umod_counts_2_CHH_deletion, 0)*100], color=colors)
                fq[1].set_xticks([0,1,2,3,4,5])
                fq[1].set_xticklabels(['CpG insert', 'CHG insert', 'CHH insert', 'CpG deletion', 'CHG deletion', 'CHH deletion'],  fontsize=12, rotation=90)
                fq[1].grid(color='lightgray', linestyle='solid', linewidth=1)
            else:
                fq.set_ylabel('Relative abundance, %', fontsize=12)
                fq.set_xlabel('Indel type', fontsize=12)
                fq.bar([0,1,2,3,4,5], [non_zero_division(mod_counts_1_CpG_insert, umod_counts_1_CpG_insert, 0)*100, non_zero_division(mod_counts_1_CHG_insert, umod_counts_1_CHG_insert, 0)*100, non_zero_division(mod_counts_1_CHH_insert, umod_counts_1_CHH_insert, 0)*100, non_zero_division(mod_counts_1_CpG_deletion, umod_counts_1_CpG_deletion, 0)*100, non_zero_division(mod_counts_1_CHG_deletion, umod_counts_1_CHG_deletion, 0)*100, non_zero_division(mod_counts_1_CHH_deletion, umod_counts_1_CHH_deletion, 0)*100], color=colors)
                fq.set_xticks([0,1,2,3,4,5])
                fq.set_xticklabels(['CpG insert', 'CHG insert', 'CHH insert', 'CpG deletion', 'CHG deletion', 'CHH deletion'],  fontsize=12, rotation=90)
                fq.grid(color='lightgray', linestyle='solid', linewidth=1)
            plt.savefig(directory + name + '_ID-bias_abundance_10bp_mod_site_plot.pdf', figsize=(16, 20), dpi=330, bbox_inches='tight', pad_inches=0.25)
            plt.close()  
              
            y_axis_CpG1_insert, y_axis_CHG1_insert, y_axis_CHH1_insert, y_axis_CpG2_insert, y_axis_CHG2_insert, y_axis_CHH2_insert = list(), list(), list(), list(), list(), list()
            y_axis_CpG1_deletion, y_axis_CHG1_deletion, y_axis_CHH1_deletion, y_axis_CpG2_deletion, y_axis_CHG2_deletion, y_axis_CHH2_deletion = list(), list(), list(), list(), list(), list()
            for row in values_1_CpG_insert:
                y_axis_CpG1_insert.append(row[1])
            for row in values_1_CHG_insert:
                y_axis_CHG1_insert.append(row[1])
            for row in values_2_CHH_insert:
                y_axis_CHH1_insert.append(row[1])
            for row in values_2_CpG_insert:
                y_axis_CpG2_insert.append(row[1])
            for row in values_2_CHG_insert:
                y_axis_CHG2_insert.append(row[1])
            for row in values_2_CHH_insert:
                y_axis_CHH2_insert.append(row[1])
            for row in values_1_CpG_deletion:
                y_axis_CpG1_deletion.append(row[1])
            for row in values_1_CHG_deletion:
                y_axis_CHG1_deletion.append(row[1])
            for row in values_2_CHH_deletion:
                y_axis_CHH1_deletion.append(row[1])
            for row in values_2_CpG_deletion:
                y_axis_CpG2_deletion.append(row[1])
            for row in values_2_CHG_deletion:
                y_axis_CHG2_deletion.append(row[1])
            for row in values_2_CHH_deletion:
                y_axis_CHH2_deletion.append(row[1])
            x_axis = [x for x in range(1,read_length+1)]
            plt.figure(figsize=(16, 12))
            if single_end == False:
                fig, fq = plt.subplots(2, 1,  sharey=True)
            else:
                fig, fq = plt.subplots(1, 1)
            fig.suptitle('Sequencing ID-bias: indel co-localisation at 10bp from a modified position', fontsize=14)
            if single_end == False:
                plt.subplots_adjust(hspace=0.4)
                plt.subplots_adjust(right=1)
                fq[0].set_ylabel('Indel rate, %', fontsize=12)
                fq[0].set_xlabel('First in pair base positions', fontsize=12)
                fq[0].plot(x_axis, y_axis_CpG1_insert, linewidth=1.5, linestyle='-', color=colors[0])
                fq[0].plot(x_axis, y_axis_CHG1_insert, linewidth=1.5, linestyle='-', color=colors[1])
                fq[0].plot(x_axis, y_axis_CHH1_insert, linewidth=1.5, linestyle='-', color=colors[2])
                fq[0].plot(x_axis, y_axis_CpG1_deletion, linewidth=1.5, linestyle='-', color=colors[3])
                fq[0].plot(x_axis, y_axis_CHG1_deletion, linewidth=1.5, linestyle='-', color=colors[4])
                fq[0].plot(x_axis, y_axis_CHH1_deletion, linewidth=1.5, linestyle='-', color=colors[5])
                fq[0].xaxis.set_ticks(numpy.arange(0, read_length + 1, step=ceil(read_length / 10)))
                fq[0].grid(color='lightgray', linestyle='solid', linewidth=1)
                fq[1].set_ylabel('Indel rate, %', fontsize=12)
                fq[1].set_xlabel('Second in pair base positions', fontsize=12)
                fq[1].plot(x_axis, y_axis_CpG2_insert, linewidth=1.5, linestyle='-', color=colors[0])
                fq[1].plot(x_axis, y_axis_CHG2_insert, linewidth=1.5, linestyle='-', color=colors[1])
                fq[1].plot(x_axis, y_axis_CHH2_insert, linewidth=1.5, linestyle='-', color=colors[2])
                fq[1].plot(x_axis, y_axis_CpG2_deletion, linewidth=1.5, linestyle='-', color=colors[3])
                fq[1].plot(x_axis, y_axis_CHG2_deletion, linewidth=1.5, linestyle='-', color=colors[4])
                fq[1].plot(x_axis, y_axis_CHH2_deletion, linewidth=1.5, linestyle='-', color=colors[5])
                fq[1].xaxis.set_ticks(numpy.arange(0, read_length + 1, step=ceil(read_length/10)))
                fq[1].grid(color='lightgray', linestyle='solid', linewidth=1)
                plt.figlegend(['CpG insert', 'CHG insert', 'CHH insert', 'CpG deletion', 'CHG deletion', 'CHH deletion'], loc='center left', bbox_to_anchor=(1.1, 0.5))
            else:
                fq.set_ylabel('Indel rate, %', fontsize=12)
                fq.set_xlabel('Base positions', fontsize=12)
                fq.plot(x_axis, y_axis_CpG1_insert, linewidth=1.5, linestyle='-', color=colors[0])
                fq.plot(x_axis, y_axis_CHG1_insert, linewidth=1.5, linestyle='-', color=colors[1])
                fq.plot(x_axis, y_axis_CHH1_insert, linewidth=1.5, linestyle='-', color=colors[2])
                fq.plot(x_axis, y_axis_CpG1_deletion, linewidth=1.5, linestyle='-', color=colors[3])
                fq.plot(x_axis, y_axis_CHG1_deletion, linewidth=1.5, linestyle='-', color=colors[4])
                fq.plot(x_axis, y_axis_CHH1_deletion, linewidth=1.5, linestyle='-', color=colors[5])
                fq.xaxis.set_ticks(numpy.arange(0, read_length + 1, step=ceil(read_length / 10)))
                fq.grid(color='lightgray', linestyle='solid', linewidth=1)
                plt.figlegend(['CpG insert', 'CHG insert', 'CHH insert', 'CpG deletion', 'CHG deletion', 'CHH deletion'], loc='center left', bbox_to_anchor=(1, 0.5))
            plt.savefig(directory + name + '_ID-bias_modification_colocalisation_plot.pdf', figsize=(16, 12), dpi=330, bbox_inches='tight', pad_inches=0.25)
            plt.close()
            
            plt.figure(figsize=(16, 12))
            if single_end == False:
                fig, fq = plt.subplots(2, 1,  sharey=True)
            else:
                fig, fq = plt.subplots(1, 1)
            fig.suptitle('Sequencing ID-bias: indel rate per read length', fontsize=14)
            if single_end == False:
                plt.subplots_adjust(hspace=0.4)
                plt.subplots_adjust(right=1)
                fq[0].set_ylabel('Indel rate, %', fontsize=12)
                fq[0].set_xlabel('First in pair base positions', fontsize=12)
                fq[0].plot(x_axis, read1_normalised_insert_values, linewidth=1.5, linestyle='-', color=colors[1])
                fq[0].plot(x_axis, read1_normalised_deletion_values, linewidth=1.5, linestyle='-', color=colors[3])
                fq[0].xaxis.set_ticks(numpy.arange(0, read_length + 1, step=ceil(read_length / 10)))
                fq[0].grid(color='lightgray', linestyle='solid', linewidth=1)
                fq[1].set_ylabel('Indel rate, %', fontsize=12)
                fq[1].set_xlabel('Second in pair base positions', fontsize=12)
                fq[1].plot(x_axis, read2_normalised_insert_values, linewidth=1.5, linestyle='-', color=colors[1])
                fq[1].plot(x_axis, read2_normalised_deletion_values, linewidth=1.5, linestyle='-', color=colors[3])
                fq[1].xaxis.set_ticks(numpy.arange(0, read_length + 1, step=ceil(read_length/10)))
                fq[1].grid(color='lightgray', linestyle='solid', linewidth=1)
                plt.figlegend(['insert', 'deletion'], loc='center left', bbox_to_anchor=(1.05, 0.5))
            else:
                fq.set_ylabel('Indel rate, %', fontsize=12)
                fq.set_xlabel('Base positions', fontsize=12)
                fq.plot(x_axis, read1_normalised_insert_values, linewidth=1.5, linestyle='-', color=colors[1])
                fq.plot(x_axis, read1_normalised_deletion_values, linewidth=1.5, linestyle='-', color=colors[3])
                fq.xaxis.set_ticks(numpy.arange(0, read_length + 1, step=ceil(read_length / 10)))
                fq.grid(color='lightgray', linestyle='solid', linewidth=1)
                plt.figlegend(['insert', 'deletion'], loc='center left', bbox_to_anchor=(0.9, 0.5))
            plt.savefig(directory + name + '_ID-bias_indel_rate_plot.pdf', figsize=(16, 12), dpi=330, bbox_inches='tight', pad_inches=0.25)
            plt.close()
    except Exception:
        logs.error('asTair cannot output the IDbias plot.', exc_info=True)
    else:
        pass
    time_m = datetime.now()
    logs.info("asTair's ID-bias summary function finished running. {} seconds".format((
    time_m - time_b).total_seconds()))


if __name__ == '__main__':
    idbias()




