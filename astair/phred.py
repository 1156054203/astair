#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import division

import re
import pdb
import csv
import sys
import gzip
import pysam
import click
import random
import logging
import warnings
from os import path
import pkg_resources
from threading import Thread
from datetime import datetime

try:
    import matplotlib as mplot
    mplot.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    plt.style.use('seaborn-whitegrid')
    plt.ioff()
except Exception:
    warnings.warn("Matplotlib was not found, visualisation output will not be supported.", ImportWarning)

if sys.version[0] == '3':
    from queue import Queue as Queue
    from itertools import zip_longest
elif sys.version[0] == '2':
    from Queue import Queue as Queue
    from itertools import izip_longest as zip_longest
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))

from astair.safe_division import non_zero_division_NA
from astair.statistics_summary import general_statistics_summary


@click.command()
@click.option('fq1', '--fq1', '-1', required=True, help='First in pair (R1) sequencing reads file in fastq.gz format.')
@click.option('fq2', '--fq2', '-2', required=False, help='Second in pair (R2) sequencing reads file in fastq.gz format.')
@click.option('calculation_mode', '--calculation_mode', '-cm', required=False, default='means', type=click.Choice(['means', 'absolute']), help='Gives the mode of computation used for the Phred scores summary, where means runs faster. (Default means)')
@click.option('single_end', '--se', '-se', default=False, is_flag=True, required=False, help='Indicates single-end sequencing reads (Default False).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
@click.option('sample_size', '--sample_size', '-s', default=10000000, type=int, required=False, help='The number of reads to sample for the analysis. (Default 10 000 000).')
@click.option('plot', '--plot', '-p', required=False, is_flag=True, help='Phred scores will be visualised and output as a pdf file. Requires installed matplotlib.')
@click.option('minimum_score', '--minimum_score', '-q', required=False, default=15, type=int, help='Minimum Phred score used for visualisation only. (Default 15).')
@click.option('colors', '--colors', '-c', default=['skyblue', 'mediumaquamarine', 'khaki', 'lightcoral'], type=list, required=False, help="List of color values used for visualistion of A, C, G, T, they are given as color1,color2,color3,color4. Accepts valid matplotlib color names, RGB and RGBA hex strings and  single letters denoting color {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}. (Default skyblue,mediumaquamarine,khaki,lightcoral).")
def phred(fq1, fq2, calculation_mode, directory, sample_size, minimum_score, colors, plot, single_end):
    """Calculate per base (A, C, T, G) Phred scores for each strand."""
    Phred_scores_plotting(fq1, fq2, calculation_mode, directory, sample_size, minimum_score, colors, plot, single_end)


warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

#logging.basicConfig(level=logging.WARNING)
logs = logging.getLogger(__name__)

time_b = datetime.now()


def numeric_Phred_score(score):
    """Converts ASCII fastq sequencing scores to numeric scores."""
    final_score = list(pysam.qualitystring_to_array(score))
    return final_score

def clean_open_file(input_file):
    """Opens neatly and separately the fastq file as an iterator."""
    try:
        fastq_file = gzip.open(input_file, "rt")
        return fastq_file
    except Exception:
        logs.error('The input file does not exist.', exc_info=True)
        sys.exit(1)

def Phred_score_main_body(file_to_load, calculation_mode, sample_size):
    """Calculates the mean or the absolute numeric Phred scores per each base using a desired sample size for random sampling from the fastq file."""
    if sample_size is None:
        cutoff = 10000000
    else:
        cutoff = int(sample_size)
    cycle_count = 0
    read_values = []
    Phred_T_sum, Phred_T_len, Phred_C_sum, Phred_C_len, Phred_A_sum, Phred_A_len, Phred_G_sum, Phred_G_len = [0] * 8
    for line in file_to_load:
        if cycle_count == 1 or cycle_count % 4 == 1:
            read_sequence = re.findall(r"(?<!@)(?!@).*[^+|^\n].*(?!@\n)", line)[0]
        elif cycle_count == 3 or cycle_count % 4 == 3:
            read_base_qualities = re.findall(r"(?<!@)(?!@).*[^+|^\n].*(?!@\n)", line)
            read_base_qualities = numeric_Phred_score(read_base_qualities[0])
            results = [(i, j) for i, j in
                       zip_longest(read_base_qualities, read_sequence, fillvalue='BLANK')]
            thymines = [x[0] for x in results if 'T' in x[:]]
            cytosines = [x[0] for x in results if 'C' in x[:]]
            adenines = [x[0] for x in results if 'A' in x[:]]
            guanines = [x[0] for x in results if 'G' in x[:]]
            Phred_T_sum += sum(thymines)
            Phred_T_len += len(thymines)
            Phred_C_sum += sum(cytosines)
            Phred_C_len += len(cytosines)
            Phred_A_sum += sum(adenines)
            Phred_A_len += len(adenines)
            Phred_G_sum += sum(guanines)
            Phred_G_len += len(guanines)
            if calculation_mode == 'means':
                if cycle_count < cutoff * 4:
                    read_values.append(tuple((non_zero_division_NA(sum(thymines), len(thymines)),
                                          non_zero_division_NA(sum(cytosines), len(cytosines)),
                                          non_zero_division_NA(sum(adenines), len(adenines)),
                                          non_zero_division_NA(sum(guanines), len(guanines)))))
                elif cycle_count >= cutoff * 4 and cycle_count < 10 * cutoff * 4:
                    random_index = random.randint(0, len(read_values) - 1)
                    read_values[random_index] = tuple((non_zero_division_NA(sum(thymines), len(thymines)),
                                                       non_zero_division_NA(sum(cytosines), len(cytosines)),
                                                       non_zero_division_NA(sum(adenines), len(adenines)),
                                                       non_zero_division_NA(sum(guanines), len(guanines))))
                else:
                    break
            elif calculation_mode == 'absolute':
                if cycle_count < cutoff * 4:
                    read_values.append((list((thymines)), list((cytosines)), list((adenines)), list((guanines))))
                elif cycle_count >= cutoff * 4 and cycle_count < 10 * cutoff * 4:
                    random_index = random.randint(0, len(read_values) - 1)
                    read_values[random_index] = tuple(((list((thymines)), list((cytosines)), list((adenines)), list((guanines)))))
                else:
                    break
        cycle_count += 1
    return read_values

def Phred_score_statistics_calculation(input_file, sample_size, calculation_mode):
    """ Wraps the Phred_score_main_body function."""
    read_values = Phred_score_main_body(clean_open_file(input_file), calculation_mode, sample_size)
    return read_values


def Phred_score_value_return(fastq_input, sample_size, calculation_mode, fq1, fq2):
    """Returns the numeric Phred scores per base for the forward and reverse reads in pair-end sequencing."""
    if fastq_input == fq1:
        read_values_fq1 = Phred_score_statistics_calculation(fastq_input, sample_size, calculation_mode)
        return read_values_fq1
    elif fastq_input == fq2:
        read_values_fq2 = Phred_score_statistics_calculation(fastq_input, sample_size, calculation_mode)
        return read_values_fq2


def main_Phred_score_calculation_output(fq1, fq2, sample_size, directory, name, calculation_mode, single_end):
    """Runs the Phred score calculation functions in parallel for the forward and the reverse reads in pair-end and single-end sequencing."""
    threads = []
    queue = [Queue(), Queue()]
    threads.append(Thread(target=lambda que, arg1, arg2, arg3, arg4, arg5: que.put(Phred_score_value_return(arg1, arg2, arg3, arg4, arg5)),
                          args=(queue[0], fq1, sample_size, calculation_mode, fq1, fq2), ))
    if single_end == False or fq2 is not None:
        threads.append(Thread(target=lambda que, arg1, arg2, arg3, arg4, arg5: que.put(Phred_score_value_return(arg1, arg2, arg3, arg4, arg5)),
                          args=(queue[1], fq2, sample_size, calculation_mode, fq1, fq2), ))
    threads[0].start()
    if single_end == False or fq2 is not None:
        threads[1].start()
    for thread in threads:
        thread.join()
    read_values_fq1 = queue[0].get()
    if single_end == False or fq2 is not None:
        read_values_fq2 = queue[1].get()
    else:
        read_values_fq2 = None
    return read_values_fq1, read_values_fq2


def summary_statistics_output(directory, name, statistics_data, read_orientation, single_end):
    """Outputs the Phred score calculation statistics as a text file that can be visualised independently from the plotting module."""
    try:
        with open(directory + name + '_total_Phred.txt', 'a') as new_file:
            data_line = csv.writer(new_file, delimiter='\t', lineterminator='\n')
            if read_orientation == 'F':
                if single_end == False:
                    read_orientation_string = 'First in pair_'
                else:
                    read_orientation_string = '__Single-end__'
            else:
                read_orientation_string = 'Second in pair'
            data_line.writerow(["____________________________________{}____________________________________".format(read_orientation_string)])
            data_line.writerow(["______________________________________________________________________________________"])
            data_line.writerow(["mean ", "adenines: " + str(statistics_data[0][0]), "cytosines: " + str(statistics_data[1][0]),
                                "thymines: " + str(statistics_data[2][0]), "guanines: " + str(statistics_data[3][0])])
            data_line.writerow(["median ", "adenines: " + str(statistics_data[0][1]), "cytosines: " + str(statistics_data[1][1]),
                                "thymines: " + str(statistics_data[2][1]), "guanines: " + str(statistics_data[3][1])])
            data_line.writerow(["q25 ", "adenines: " + str(statistics_data[0][3]), "cytosines: " + str(statistics_data[1][3]),
                                "thymines: " + str(statistics_data[2][3]), "guanines: " + str(statistics_data[3][3])])
            data_line.writerow(["q75 ", "adenines: " + str(statistics_data[0][4]), "cytosines: " + str(statistics_data[1][4]),
                                "thymines: " + str(statistics_data[2][4]), "guanines: " + str(statistics_data[3][4])])
            data_line.writerow(["sd ", "adenines: " + str(statistics_data[0][2]), "cytosines: " + str(statistics_data[1][2]),
                                "thymines: " + str(statistics_data[2][2]), "guanines: " + str(statistics_data[3][2])])
            data_line.writerow(["min ", "adenines: " + str(statistics_data[0][5]), "cytosines: " + str(statistics_data[1][5]),
                                "thymines: " + str(statistics_data[2][5]), "guanines: " + str(statistics_data[3][5])])
            data_line.writerow(["max ", "adenines: " + str(statistics_data[0][6]), "cytosines: " + str(statistics_data[1][6]),
                                "thymines: " + str(statistics_data[2][6]), "guanines: " + str(statistics_data[3][6])])
            data_line.writerow(["______________________________________________________________________________________"])
    except IOError:
        logs.error('asTair cannot write to Phred scores summary file.', exc_info=True)


def Phred_values_return(input_values, read_orientation, directory, name, calculation_mode, single_end):
    """Takes the Phred scores per read and splits them in lists per base, and then calculates the summary statistics."""
    Ts = [x[0] for x in input_values if x[0] != 'NA']
    Cs = [x[1] for x in input_values if x[1] != 'NA']
    As = [x[2] for x in input_values if x[2] != 'NA']
    Gs = [x[3] for x in input_values if x[3] != 'NA']
    if calculation_mode == 'absolute':
        Ts = sum(Ts, [])
        Cs = sum(Cs, [])
        As = sum(As, [])
        Gs = sum(Gs, [])
    As_stats = tuple(general_statistics_summary(As))
    Cs_stats = tuple(general_statistics_summary(Cs))
    Gs_stats = tuple(general_statistics_summary(Gs))
    Ts_stats = tuple(general_statistics_summary(Ts))
    statistics_data = list((As_stats, Cs_stats, Ts_stats, Gs_stats))
    summary_statistics_output(directory, name, statistics_data, read_orientation, single_end)
    data = list((As, Cs, Gs, Ts))
    mx = max(Ts + Cs + As + Gs)
    return data, mx


def Phred_scores_color_change(plot_input, colors):
    """Employs the desired colouring scheme for the plotting module."""
    for plot in plot_input:
        for box in plot['boxes']:
            box.set(color='black', linewidth=1)
        for patch, color in zip(plot['boxes'], colors):
            patch.set_facecolor(color)
        for whisker in plot['whiskers']:
            whisker.set(color='black', linewidth=1)
        for median in plot['medians']:
            median.set(color='black', linewidth=1)
        for flier in plot['fliers']:
            flier.set(marker='', markersize=1, markerfacecolor=None, markeredgecolor=None)


def Phred_scores_plotting(fq1, fq2, calculation_mode, directory, sample_size, minimum_score, colors, plot, single_end):
    """The general function that takes the input file, calculates mean or absolute Phred scores per base, and outputs
    their summary statistics as a text file or as a plot when the plotting module is enabled."""
    time_s = datetime.now()
    logs.info("asTair's Phred scores statistics summary function started running. {} seconds".format((time_s - time_b).total_seconds()))
    try:
        file1 = open(fq1, 'r')
        file1.close()
        if single_end == False or fq2 is not None:
            file2 = open(fq2, 'r')
            file2.close()
    except Exception:
        logs.error('The input fastq files do not exist.', exc_info=True)
        sys.exit(1)
    name = path.splitext(path.basename(fq1))[0]
    if single_end == False or fq2 is not None:
        name = re.sub('_(R1|1).fq', '', name)
    else:
        name = re.sub('.fq', '', name)
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    if path.exists(directory) == False:
        raise Exception("The output directory does not exist.")
        sys.exit(1)
    if colors != ['skyblue', 'mediumaquamarine', 'khaki', 'lightcoral']:
        colors = "".join(colors).split(',')
    read_values_fq1, read_values_fq2 = main_Phred_score_calculation_output(fq1, fq2, sample_size, directory, name, calculation_mode, single_end)
    data_fq1, maxy1 = Phred_values_return(read_values_fq1, 'F', directory, name, calculation_mode, single_end)
    if single_end == False  or fq2 is not None:
        data_fq2, maxy2 = Phred_values_return(read_values_fq2, 'R', directory, name, calculation_mode, single_end)
    try:
        if plot:
            plt.figure()
            if single_end == False  or fq2 is not None:
                fig, fq = plt.subplots(1, 2)
                fig.suptitle('Sequencing base quality', fontsize=14)
                plt.subplots_adjust(wspace=0.4)
                maxy = [max(maxy1, maxy2) + 1 if max(maxy1, maxy2) + 1 > 35 else 35][0]
                box1 = fq[0].boxplot(data_fq1, labels=['A', 'C', 'G', 'T'], patch_artist=True)
                fq[0].set_ylabel('Phred score', fontsize=12)
                fq[0].set_xlabel('First in pair', fontsize=12)
                fq[0].axis([0, 5, minimum_score, maxy])
                fq[0].yaxis.set_major_locator(ticker.MultipleLocator(5))
                fq[0].grid(color='lightgray', linestyle='solid', linewidth=1)
                box2 = fq[1].boxplot(data_fq2, labels=['A', 'C', 'G', 'T'], patch_artist=True)
                fq[1].set_xlabel('Second in pair', fontsize=12)
                fq[1].axis([0, 5, minimum_score, maxy])
                fq[1].yaxis.set_major_locator(ticker.MultipleLocator(5))
                fq[1].grid(color='lightgray', linestyle='solid', linewidth=1)
                Phred_scores_color_change([box1, box2], colors)
                plt.vlines(-1, minimum_score, maxy, alpha=0.3, linewidth=1, linestyle='--', color='gray', clip_on=False)
                plt.savefig(directory + name + '_phred_scores_plot.pdf', figsize=(14, 8), dpi=330, bbox_inches='tight')
            else:
                fig, fq = plt.subplots(1, 1)
                fig.suptitle('Sequencing base quality', fontsize=14)
                maxy = [maxy1 + 1 if maxy1 + 1 > 35 else 35][0]
                box1 = fq.boxplot(data_fq1, labels=['A', 'C', 'G', 'T'], patch_artist=True)
                fq.set_ylabel('Phred score', fontsize=12)
                fq.set_xlabel('Single-end read', fontsize=12)
                fq.axis([0, 5, minimum_score, maxy])
                fq.yaxis.set_major_locator(ticker.MultipleLocator(5))
                fq.grid(color='lightgray', linestyle='solid', linewidth=1)
                Phred_scores_color_change([box1], colors)
                plt.savefig(directory + name + '_phred_scores_plot.pdf', figsize=(14, 8), dpi=330, bbox_inches='tight')
            plt.close()
    except Exception:
        logs.error('asTair cannot output the Phred scores plot.', exc_info=True)
    else:
        pass
    time_m = datetime.now()
    logs.info("asTair's Phred scores statistics summary function finished running. {} seconds".format((
    time_m - time_b).total_seconds()))


if __name__ == '__main__':
    phred()
