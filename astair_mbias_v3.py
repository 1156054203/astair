#!/usr/bin/env python3

import click
import pysam
import re
import sys
import itertools
import csv
from math import ceil
import numpy
from datetime import datetime
import logging
import warnings
from os import path
try:
    import matplotlib as mplot
    mplot.use('Agg')
    import matplotlib.pyplot as pyp
    import matplotlib.ticker as ticker
    pyp.style.use('seaborn-whitegrid')
    pyp.ioff()
except ImportError:
    pass

from safe_division import non_zero_division



@click.command()
@click.option('input_file', '--input_file', required=True, help='BAM format file containing sequencing reads.')
@click.option('directory', '--directory', required=True, help='Output directory to save files.')
@click.option('read_length', '--read_length', type=int, required=True, help='The read length is needed to calculate the M-bias.')
@click.option('plot', '--plot', required=False, is_flag=True, help='Phred scores will be visualised and output as a pdf file. Requires installed matplotlib.')
def Mbias_exec(input_file, directory, read_length, plot):
    Mbias_plotting(input_file, directory, read_length, plot)


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()

def initialise_data_counters(read_length):
    """Clean initialisation of empty dictionaries used for counters."""
    all_read_data = list(({}, {}, {}))
    for read_data in all_read_data:
        for i in range(0, read_length):
            read_data[i] = 0
    return all_read_data[0], all_read_data[1], all_read_data[2]


def bam_file_opener(input_file):
    """Opens neatly and separately the bam file as an iterator."""
    try:
        open(input_file, 'rb')
    except (SystemExit, KeyboardInterrupt, IOError, FileNotFoundError):
        logs.error('The input file does not exist.', exc_info=True)
        sys.exit(1)
    inbam = pysam.AlignmentFile(input_file, "rb")
    bam_fetch = inbam.fetch(until_eof=True)
    return bam_fetch

def check_read_info(read):
    """Checks if the MD column exist in the input bam file and takes its string."""
    try:
        if isinstance(read.tags[0][1], str) and read.tags[0][0] == 'MD':
            data = read.tags[0][1]
        elif isinstance(read.tags[1][1], str) and read.tags[1][0] == 'MD':
            data = read.tags[1][1]
        return data
    except (IndexError, TypeError):
            logs.error('The input file does not contain a MD tag column.', exc_info=True)
            sys.exit(1)


def mbias_calculator(read, read_length, read_mods_CpG, read_mods_CHG, read_mods_CHH, read_all_CpG, read_all_CHG, read_all_CHH):
    """Calculates the modification level per read position, pair orientation and cytosine context."""
    read_data = check_read_info(read)
    if re.search("(?:.*C.*)", read_data, re.IGNORECASE):
        changes = [int(s) for s in re.findall(r'\d+', read_data)]
        non_overlap = [x + 1 if x == 0 else x for x in changes]
        names = list(re.findall(r'[^\W\d_]+', read_data))
        positions = [x + 1 for x in list(itertools.accumulate(non_overlap))]
        sequence = list(read.query_sequence)
        positions = [x for x in positions if x <= len(sequence)]
        for k in range(len(positions) - 1):
            if len(names[k]) == 1:
                element = positions[k]
                sequence[element] = names[k]
        reads = "".join(sequence)
        cpg_all = [m.start() for m in re.finditer(r'CG', reads, re.IGNORECASE)]
        chg_all = [m.start() for m in re.finditer(r'C(A|C|T)G', reads, re.IGNORECASE)]
        chh_all = [m.start() for m in re.finditer(r'C(A|C|T)(A|T|C)', reads, re.IGNORECASE)]
        cpg_mods = [x for x in positions if x in cpg_all]
        chg_mods = [x for x in positions if x in chg_all]
        chh_mods = [x for x in positions if x in chh_all]
        if len(reads) <= read_length:
            for i in range(0, len(reads)):
                if i in chh_mods:
                    read_mods_CHH[i] += 1
                elif i in chg_mods:
                    read_mods_CHG[i] += 1
                elif i in cpg_mods:
                    read_mods_CpG[i] += 1
                elif i in chh_all:
                    read_all_CHH[i] += 1
                elif i in chg_all:
                    read_all_CHG[i] += 1
                elif i in cpg_all:
                    read_all_CpG[i] += 1
    else:
        sequence = list(read.query_sequence)
        reads = "".join(sequence)
        cpg_all = [m.start() for m in re.finditer(r'CG', reads, re.IGNORECASE)]
        chg_all = [m.start() for m in re.finditer(r'C(A|C|T)G', reads, re.IGNORECASE)]
        chh_all = [m.start() for m in re.finditer(r'C(A|C|T)(A|T|C)', reads, re.IGNORECASE)]
        if len(reads) <= read_length:
            for i in range(0, len(reads)):
                if i in chh_all:
                    read_all_CHH[i] += 1
                if i in chg_all:
                    read_all_CHG[i] += 1
                if i in cpg_all:
                    read_all_CpG[i] += 1
    return read_mods_CpG, read_mods_CHG, read_mods_CHH, read_all_CpG, read_all_CHG, read_all_CHH


def mbias_evaluater(input_file, read_length):
    """Outputs the modification levels per read position, pair orientation and cytosine context."""
    read1_mods_CHH, read1_mods_CHG, read1_mods_CpG = initialise_data_counters(read_length)
    read1_all_CHH, read1_all_CHG, read1_all_CpG = initialise_data_counters(read_length)
    read2_mods_CHH, read2_mods_CHG, read2_mods_CpG = initialise_data_counters(read_length)
    read2_all_CHH, read2_all_CHG, read2_all_CpG = initialise_data_counters(read_length)
    for read in bam_file_opener(input_file):
        if read.flag == 99 or read.flag == 83:
            mbias_calculator(read, read_length, read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_all_CpG, read1_all_CHG, read1_all_CHH)
        if read.flag == 147 or read.flag == 163:
            mbias_calculator(read, read_length, read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_all_CpG, read2_all_CHG, read2_all_CHH)
    return read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_all_CpG, read1_all_CHG, read1_all_CHH,\
           read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_all_CpG, read2_all_CHG, read2_all_CHH


def mbias_statistics_calculator(input_file, name, directory, read_length):
    """Creates a summary statistics of the modification levels per read position, pair orientation and cytosine context,
    and then writes them as a text file that can be used for independent visualisation."""
    read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_all_CpG, read1_all_CHG, read1_all_CHH,\
    read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_all_CpG, read2_all_CHG, read2_all_CHH \
        = mbias_evaluater(input_file, read_length)
    read_values_1_CpG, read_values_1_CHG, read_values_1_CHH = initialise_data_counters(read_length)
    read_values_2_CpG, read_values_2_CHG, read_values_2_CHH = initialise_data_counters(read_length)
    for i in range(0, read_length):
        read_values_1_CHH[i] = non_zero_division(read1_mods_CHH[i], read1_all_CHH[i]+read1_mods_CHH[i])*100
        read_values_1_CHG[i] = non_zero_division(read1_mods_CHG[i], read1_all_CHG[i]+read1_mods_CHG[i])*100
        read_values_1_CpG[i] = non_zero_division(read1_mods_CpG[i], read1_all_CpG[i]+read1_mods_CpG[i])*100
        values_1_CHH = [(keys + 1, round(values[0], 3)) if isinstance(values, list) else (keys + 1, round(values, 3)) for keys, values in read_values_1_CHH.items()]
        values_1_CHG = [(keys + 1, round(values[0], 3)) if isinstance(values, list) else (keys + 1, round(values, 3)) for keys, values in read_values_1_CHG.items()]
        values_1_CpG = [(keys + 1, round(values[0], 3)) if isinstance(values, list) else (keys + 1, round(values, 3)) for keys, values in read_values_1_CpG.items()]
        read_values_2_CHH[i] = non_zero_division(read2_mods_CHH[i], read2_all_CHH[i]+read2_mods_CHH[i])*100
        read_values_2_CHG[i] = non_zero_division(read2_mods_CHG[i], read2_all_CHG[i]+read2_mods_CHG[i])*100
        read_values_2_CpG[i] = non_zero_division(read2_mods_CpG[i], read2_all_CpG[i]+read2_mods_CpG[i])*100
        values_2_CHH = [(keys + 1, round(values[0], 3)) if isinstance(values, list) else (keys + 1, round(values, 3)) for keys, values in read_values_2_CHH.items()]
        values_2_CHG = [(keys + 1, round(values[0], 3)) if isinstance(values, list) else (keys + 1, round(values, 3)) for keys, values in read_values_2_CHG.items()]
        values_2_CpG = [(keys + 1, round(values[0], 3)) if isinstance(values, list) else (keys + 1, round(values, 3)) for keys, values in read_values_2_CpG.items()]
    all_values = [(*a1, *a2, *a3, *a4, *a5, *a6) for a1, a2, a3, a4, a5, a6 in itertools.zip_longest(values_1_CpG, values_2_CpG, values_1_CHG, values_2_CHG, values_1_CHH, values_2_CHH)]
    all_values = [(x[0], x[1], x[3], x[5], x[7], x[9], x[11]) for x in all_values]
    with open(directory + name + ".Mbias.txt", 'w', newline='') as stats_file:
        line = csv.writer(stats_file, delimiter='\t', lineterminator='\n')
        line.writerow(['POSITION (bp)', 'CpG READ 1', 'CpG READ 2', 'CHG READ 1', 'CHG READ 2', 'CHH READ 1', 'CHH READ 2'])
        for row in all_values:
            line.writerow(row)
    return values_1_CpG, values_2_CpG, values_1_CHG, values_2_CHG, values_1_CHH, values_2_CHH


def Mbias_plotting(input_file, directory, read_length, plot):
    """The general M-bias calculation and statistics output function, which might be also visualised if the plotting module is enabled."""
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    values_1_CpG, values_2_CpG, values_1_CHG, values_2_CHG, values_1_CHH, values_2_CHH = mbias_statistics_calculator(input_file, name, directory, read_length)
    if plot:
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
        fq[0].plot(x_axis, y_axis_CpG1, linewidth=1.0, linestyle='-', color='teal')
        fq[0].plot(x_axis, y_axis_CHG1, linewidth=1.0, linestyle='-', color='gray')
        fq[0].plot(x_axis, y_axis_CHH1, linewidth=1.0, linestyle='-', color='maroon')
        fq[0].xaxis.set_ticks(numpy.arange(0, read_length + 1, step=ceil(read_length/10)))
        fq[0].yaxis.set_ticks(numpy.arange(0, 101, step=10))
        fq[0].grid(color='lightgray', linestyle='solid', linewidth=1)
        fq[1].set_ylabel('Modification level, %', fontsize=12)
        fq[1].set_xlabel('Second in pair base positions', fontsize=12)
        fq[1].plot(x_axis, y_axis_CpG2, linewidth=1.0, linestyle='-', color='teal')
        fq[1].plot(x_axis, y_axis_CHG2, linewidth=1.0, linestyle='-', color='gray')
        fq[1].plot(x_axis, y_axis_CHH2, linewidth=1.0, linestyle='-', color='maroon')
        fq[1].xaxis.set_ticks(numpy.arange(0, read_length + 1, step=ceil(read_length/10)))
        fq[1].yaxis.set_ticks(numpy.arange(0, 101, step=10))
        fq[1].grid(color='lightgray', linestyle='solid', linewidth=1)
        pyp.figlegend(['CpG', 'CHG', 'CHH'], loc='center left', bbox_to_anchor=(1, 0.5))
        pyp.savefig(directory + name + '_M-bias_plot.pdf', figsize=(16, 12), dpi=330, bbox_inches='tight', pad_inches=0.15)
        pyp.close()
    else:
        pass
    time_m = datetime.now()
    logs.info("asTair's M-bias summary function has finished running. {} seconds".format((
    time_m - time_b).total_seconds()))


if __name__ == '__main__':
    Mbias_exec()




