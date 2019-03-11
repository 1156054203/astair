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
import pkg_resources
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
    import matplotlib.pyplot as pyp
    import matplotlib.ticker as ticker
    pyp.style.use('seaborn-whitegrid')
    pyp.ioff()
except Exception:
    warnings.warn("Matplotlib was not found, visualisation output will not be supported.", ImportWarning)

from astair.safe_division import non_zero_division
from astair.bam_file_parser import bam_file_opener
from astair.DNA_sequences_operations import complementary


@click.command()
@click.option('input_file', '--input_file', '-i', required=True, help='BAM|CRAM format file containing sequencing reads.')
@click.option('directory', '--directory', '-d', required=True, help='Output directory to save files.')
@click.option('read_length', '--read_length', '-l', type=int, required=True, help='The read length is needed to calculate the M-bias.')
@click.option('method', '--method', '-m',  required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like).')
@click.option('plot', '--plot', '-p', required=False, is_flag=True, help='Phred scores will be visualised and output as a pdf file. Requires installed matplotlib.')
@click.option('colors', '--colors', '-c', default=['teal', 'gray', 'maroon'], type=list, required=False, help="List of color values used for visualistion of CpG, CHG and CHH modification levels per read, which are given as color1,color2,color3. Accepts valid matplotlib color names, RGB and RGBA hex strings and  single letters denoting color {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}. (Default 'teal','gray','maroon')")
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
def mbias(input_file, directory, read_length, method, plot, colors, N_threads):
    """Generate modification per read length information (Mbias). This is a quality-control measure."""
    Mbias_plotting(input_file, directory, read_length, method, plot, colors, N_threads)


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

#logging.basicConfig(level=logging.WARNING)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def initialise_data_counters(read_length):
    """Clean initialisation of empty dictionaries used for counters."""
    all_read_data = list(({}, {}, {}))
    for read_data in all_read_data:
        for i in range(0, read_length):
            read_data[i] = 0
    return all_read_data[0], all_read_data[1], all_read_data[2]


def strand_and_method(flag, ref_sequence, read_sequence, method):
    """Takes the positions of interest in the read given the flag and the method."""
    try:
        if flag == 99 or flag == 147:
            cytosines_reference = [m.start() for m in re.finditer(r'C', ref_sequence, re.IGNORECASE)]
            if method == 'mCtoT':
                thymines_read = [m.start() for m in re.finditer(r'T', read_sequence, re.IGNORECASE)]
                positions = list(set(thymines_read).intersection(set(cytosines_reference)))
            elif method == 'CtoT':
                thymines_read = [m.start() for m in re.finditer(r'T', read_sequence, re.IGNORECASE)]
                positions = list(set(cytosines_reference).difference(set(thymines_read)))
        elif flag == 163 or flag == 83:
            guanines_reference = [m.start() for m in re.finditer(r'G', ref_sequence, re.IGNORECASE)]
            if method == 'mCtoT':
                adenines_read = [m.start() for m in re.finditer(r'A', read_sequence, re.IGNORECASE)]
                positions = list(set(adenines_read).intersection(set(guanines_reference)))
            elif method == 'CtoT':
                adenines_read = [m.start() for m in re.finditer(r'A', read_sequence, re.IGNORECASE)]
                positions = list(set(guanines_reference).difference(set(adenines_read)))
        return positions
    except (IndexError, TypeError, ValueError):
            logs.error('The input file does not contain a MD tag column.', exc_info=True)
            sys.exit(1)


def mbias_calculator(flag, ref_sequence, read_sequence, read_length, read_mods_CpG, read_mods_CHG, read_mods_CHH, read_umod_CpG, read_umod_CHG, read_umod_CHH, method):
    """Calculates the modification level per read position, pair orientation and cytosine context."""
    positions = strand_and_method(flag, ref_sequence, read_sequence, method)
    if flag == 99 or flag == 147:
        cpg_all = [m.start() for m in re.finditer(r'CG', ref_sequence, re.IGNORECASE)]
        chg_all = [m.start() for m in re.finditer(r'C(A|C|T)G', ref_sequence, re.IGNORECASE)]
        chh_all = [m.start() for m in re.finditer(r'C(A|C|T)(A|T|C)', ref_sequence, re.IGNORECASE)]
    elif flag == 163 or flag == 83:
        cpg_all = [m.start() + 1 for m in re.finditer(r'CG', ref_sequence, re.IGNORECASE)]
        chg_all = [m.start() for m in re.finditer(r'G(A|G|T)C', complementary(ref_sequence), re.IGNORECASE)]
        chh_all = [m.start() for m in re.finditer(r'G(A|G|T)(A|T|G)', complementary(ref_sequence), re.IGNORECASE)]
    if len(positions) >= 1:
        cpg_mods = [x for x in positions if x in cpg_all]
        chg_mods = [x for x in positions if x in chg_all]
        chh_mods = [x for x in positions if x in chh_all]
        if len(read_sequence) <= read_length:
            for i in range(0, len(read_sequence)):
                if i in chh_mods:
                    read_mods_CHH[i] += 1
                elif i in chg_mods:
                    read_mods_CHG[i] += 1
                elif i in cpg_mods:
                    read_mods_CpG[i] += 1
                elif i in chh_all:
                    read_umod_CHH[i] += 1
                elif i in chg_all:
                    read_umod_CHG[i] += 1
                elif i in cpg_all:
                    read_umod_CpG[i] += 1
    else:
        if len(read_sequence) <= read_length:
            for i in range(0, len(read_sequence)):
                if i in chh_all:
                    read_umod_CHH[i] += 1
                if i in chg_all:
                    read_umod_CHG[i] += 1
                if i in cpg_all:
                    read_umod_CpG[i] += 1
    return read_mods_CpG, read_mods_CHG, read_mods_CHH, read_umod_CpG, read_umod_CHG, read_umod_CHH


def mbias_evaluater(input_file, read_length, method, N_threads):
    """Outputs the modification levels per read position, pair orientation and cytosine context."""
    read1_mods_CHH, read1_mods_CHG, read1_mods_CpG = initialise_data_counters(read_length)
    read1_umod_CHH, read1_umod_CHG, read1_umod_CpG = initialise_data_counters(read_length)
    read2_mods_CHH, read2_mods_CHG, read2_mods_CpG = initialise_data_counters(read_length)
    read2_umod_CHH, read2_umod_CHG, read2_umod_CpG = initialise_data_counters(read_length)
    for read in bam_file_opener(input_file, 'fetch', N_threads):
        if read.flag == 83 or read.flag == 99:
            mbias_calculator(read.flag, read.get_reference_sequence(), read.query_sequence, read_length,
                             read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH, method)
        elif read.flag == 163 or read.flag == 147:
            mbias_calculator(read.flag, read.get_reference_sequence(), read.query_sequence, read_length,
                             read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH, method)
    return read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH,\
           read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH


def context_calculator(i, read_mods, read_umods, read_values):
    """Calculates summary statistics per context and read orientation."""
    read_values[i] = non_zero_division(read_mods[i], read_umods[i] + read_mods[i]) * 100
    values = [(keys + 1, round(vals[0], 3)) if isinstance(vals, list) else (keys + 1, round(vals, 3)) for keys, vals in read_values.items()]
    umod_counts = [(keys + 1, round(vals[0], 3)) if isinstance(vals, list) else (keys + 1, round(vals, 3)) for keys, vals in read_umods.items()]
    mod_counts = [(keys + 1, round(vals[0], 3)) if isinstance(vals, list) else (keys + 1, round(vals, 3)) for keys, vals in read_mods.items()]
    return read_values, values, umod_counts, mod_counts


def mbias_statistics_calculator(input_file, name, directory, read_length, method, N_threads):
    """Creates a summary statistics of the modification levels per read position, pair orientation and cytosine context,
    and then writes them as a text file that can be used for independent visualisation."""
    read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH,\
    read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH \
        = mbias_evaluater(input_file, read_length, method, N_threads)
    read_values_1_CpG, read_values_1_CHG, read_values_1_CHH = initialise_data_counters(read_length)
    read_values_2_CpG, read_values_2_CHG, read_values_2_CHH = initialise_data_counters(read_length)
    for i in range(0, read_length):
        read_values_1_CHH, values_1_CHH, umod_counts_1_CHH, mod_counts_1_CHH = context_calculator(i, read1_mods_CHH, read1_umod_CHH, read_values_1_CHH)
        read_values_1_CHG, values_1_CHG, umod_counts_1_CHG, mod_counts_1_CHG = context_calculator(i, read1_mods_CHG, read1_umod_CHG, read_values_1_CHG)
        read_values_1_CpG, values_1_CpG, umod_counts_1_CpG, mod_counts_1_CpG = context_calculator(i, read1_mods_CpG, read1_umod_CpG, read_values_1_CpG)
        read_values_2_CHH, values_2_CHH, umod_counts_2_CHH, mod_counts_2_CHH = context_calculator(i, read2_mods_CHH, read2_umod_CHH, read_values_2_CHH)
        read_values_2_CHG, values_2_CHG, umod_counts_2_CHG, mod_counts_2_CHG = context_calculator(i, read2_mods_CHG, read2_umod_CHG, read_values_2_CHG)
        read_values_2_CpG, values_2_CpG, umod_counts_2_CpG, mod_counts_2_CpG = context_calculator(i, read2_mods_CpG, read2_umod_CpG, read_values_2_CpG)
    all_values = [(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18) for
                  a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18 in
                  zip_longest(values_1_CpG, umod_counts_1_CpG, mod_counts_1_CpG, values_2_CpG,
                                            umod_counts_2_CpG, mod_counts_2_CpG,values_1_CHG, umod_counts_1_CHG,
                                            mod_counts_1_CHG,  values_2_CHG, umod_counts_2_CHG, mod_counts_2_CHG,
                                            values_1_CHH, umod_counts_1_CHH, mod_counts_1_CHH,  values_2_CHH, umod_counts_2_CHH, mod_counts_2_CHH)]
    all_values = [(all_values[i][0][0], all_values[i][0][1], all_values[i][1][1], all_values[i][2][1], all_values[i][3][1],
                   all_values[i][4][1], all_values[i][5][1], all_values[i][6][1], all_values[i][7][1], all_values[i][8][1],
                   all_values[i][9][1], all_values[i][10][1], all_values[i][11][1], all_values[i][12][1], all_values[i][13][1],
                   all_values[i][14][1], all_values[i][15][1], all_values[i][16][1], all_values[i][17][1]) for i in range(0, len(all_values))]
    try:
        with open(directory + name + "_Mbias.txt", 'w') as stats_file:
            line = csv.writer(stats_file, delimiter='\t', lineterminator='\n')
            line.writerow(['POSITION_(bp)', 'MOD_LVL_CpG_READ_1', 'UNMOD_COUNT_CpG_READ_1', 'MOD_COUNT_CpG_READ_1',
                           'MOD_LVL_CpG_READ_2', 'UNMOD_COUNT_CpG_READ_2', 'MOD_COUNT_CpG_READ_2',
                           'MOD_LVL_CHG_READ_1', 'UNMOD_COUNT_CHG_READ_1', 'MOD_COUNT_CHG_READ_1',
                           'MOD_LVL_CHG_READ_2', 'UNMOD_COUNT_CHG_READ_2', 'MOD_COUNT_CHG_READ_2',
                           'MOD_LVL_CHH_READ_1', 'UNMOD_COUNT_CHH_READ_1', 'MOD_COUNT_CHH_READ_1',
                           'MOD_LVL_CHH_READ_2', 'UNMOD_COUNT_CHH_READ_2', 'MOD_COUNT_CHH_READ_2'])
            for row in all_values:
                line.writerow(row)
    except IOError:
        logs.error('asTair cannot write to Mbias file.', exc_info=True)
    return values_1_CpG, values_2_CpG, values_1_CHG, values_2_CHG, values_1_CHH, values_2_CHH

def Mbias_plotting(input_file, directory, read_length, method, plot, colors, N_threads):
    """The general M-bias calculation and statistics output function, which might be also visualised if the plotting module is enabled."""
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    values_1_CpG, values_2_CpG, values_1_CHG, values_2_CHG, values_1_CHH, values_2_CHH = mbias_statistics_calculator(input_file, name, directory, read_length, method, N_threads)
    try:
        if plot:
            if colors != ['teal', 'gray', 'maroon']:
                colors = "".join(colors).split(',')
            y_axis_CpG1, y_axis_CHG1, y_axis_CHH1, y_axis_CpG2, y_axis_CHG2, y_axis_CHH2 = list(), list(), list(), list(), list(), list()
            for row in values_1_CpG:
                y_axis_CpG1.append(row[1])
            for row in values_1_CHG:
                y_axis_CHG1.append(row[1])
            for row in values_1_CHH:
                y_axis_CHH1.append(row[1])
            for row in values_2_CpG:
                y_axis_CpG2.append(row[1])
            for row in values_2_CHG:
                y_axis_CHG2.append(row[1])
            for row in values_2_CHH:
                y_axis_CHH2.append(row[1])
            x_axis = [x for x in range(1,read_length+1)]
            pyp.figure()
            fig, fq = pyp.subplots(2, 1)
            fig.suptitle('Sequencing M-bias', fontsize=14)
            pyp.subplots_adjust(hspace=0.4)
            pyp.subplots_adjust(right=1)
            fq[0].set_ylabel('Modification level, %', fontsize=12)
            fq[0].set_xlabel('First in pair base positions', fontsize=12)
            fq[0].plot(x_axis, y_axis_CpG1, linewidth=1.0, linestyle='-', color=colors[0])
            fq[0].plot(x_axis, y_axis_CHG1, linewidth=1.0, linestyle='-', color=colors[1])
            fq[0].plot(x_axis, y_axis_CHH1, linewidth=1.0, linestyle='-', color=colors[2])
            fq[0].xaxis.set_ticks(numpy.arange(0, read_length + 1, step=ceil(read_length/10)))
            fq[0].yaxis.set_ticks(numpy.arange(0, 101, step=10))
            fq[0].grid(color='lightgray', linestyle='solid', linewidth=1)
            fq[1].set_ylabel('Modification level, %', fontsize=12)
            fq[1].set_xlabel('Second in pair base positions', fontsize=12)
            fq[1].plot(x_axis, y_axis_CpG2, linewidth=1.0, linestyle='-', color=colors[0])
            fq[1].plot(x_axis, y_axis_CHG2, linewidth=1.0, linestyle='-', color=colors[1])
            fq[1].plot(x_axis, y_axis_CHH2, linewidth=1.0, linestyle='-', color=colors[2])
            fq[1].xaxis.set_ticks(numpy.arange(0, read_length + 1, step=ceil(read_length/10)))
            fq[1].yaxis.set_ticks(numpy.arange(0, 101, step=10))
            fq[1].grid(color='lightgray', linestyle='solid', linewidth=1)
            pyp.figlegend(['CpG', 'CHG', 'CHH'], loc='center left', bbox_to_anchor=(1, 0.5))
            pyp.savefig(directory + name + '_M-bias_plot.pdf', figsize=(16, 12), dpi=330, bbox_inches='tight', pad_inches=0.15)
            pyp.close()
    except Exception:
        logs.error('asTair cannot output the Mbias plot.', exc_info=True)
    else:
        pass
    time_m = datetime.now()
    logs.info("asTair's M-bias summary function has finished running. {} seconds".format((
    time_m - time_b).total_seconds()))


if __name__ == '__main__':
    mbias()




