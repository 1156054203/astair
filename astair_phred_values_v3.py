#!/usr/bin/env python3

import re
import csv
import sys
import gzip
import click
import random
import logging
import itertools
from os import path
from queue import Queue
from threading import Thread
from datetime import datetime
try:
    import matplotlib as mplot
    mplot.use('Agg')
    import matplotlib.pyplot as pyp
    import matplotlib.ticker as ticker
    pyp.style.use('seaborn-whitegrid')
    pyp.ioff()
except ImportError:
    pass

from safe_division import non_zero_division_NA
from statistics_summary import general_statistics_summary


@click.command()
@click.option('fq1', '--fq1', '-1', required=True, help='First in pair (R1) sequencing reads file in fastq.gz format')
@click.option('fq2', '--fq2', '-2', required=True, help='Second in pair (R2) sequencing reads file in fastq.gz format')
@click.option('calculation_mode', '--calculation_mode', '-cm', required=False, default='means', type=click.Choice(['means', 'absolute']), help='Gives the mode of computation used for the Phred scores summary, where means runs faster. (Default is means)')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
@click.option('sample_size', '--sample_size', '-s', default=10000000, type=int, required=False, help='The number of reads to sample for the analysis. (Default 10 000 000)')
@click.option('plot', '--plot', '-p', required=False, is_flag=True, help='Phred scores will be visualised and output as a pdf file. Requires installed matplotlib.')
@click.option('minimum_score', '--minimum_score', '-q', required=False, default=15, type=int, help='Minimum Phred score used for visualisation only. (Default 15)')
@click.option('colors', '--colors', '-c', default=['skyblue', 'mediumaquamarine', 'khaki', 'lightcoral'], type=list, required=False, help="List of color values used for visualistion of A, C, G, T, they are given as color1,color2,color3,color4. Accepts valid matplotlib color names, RGB and RGBA hex strings and  single letters denoting color {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}. (Default skyblue,mediumaquamarine,khaki,lightcoral)")
def Phred_score_calculation_visualisation_exec(fq1, fq2, calculation_mode, directory, sample_size, minimum_score, colors, plot):
    Phred_scores_plotting(fq1, fq2, calculation_mode, directory, sample_size, minimum_score, colors, plot)


logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()


def numeric_Phred_score(score):
    """Converts ASCII fastq sequencing scores to numeric scores."""
    final_score = score.translate({ord("!"): ' 0 ', ord("\""): ' 1 ', ord("#"): ' 2 ', ord("$"): ' 3 ', ord("%"): ' 4 ',
                                   ord("&"): ' 5 ', ord("'"): ' 6 ', ord("("): ' 7 ', ord(")"): ' 8 ', ord("*"): ' 9 ',
                                   ord("+"): ' 10 ', ord("`"): ' 11 ', ord("-"): ' 12 ', ord("."): ' 13 ', ord("/"): ' 14 ',
                                   ord("0"): ' 15 ', ord("1"): ' 16 ', ord("2"): ' 17 ', ord("3"): ' 18 ', ord("4"): ' 19 ', 
                                   ord("5"): ' 20 ', ord("6"): ' 21 ', ord("7"): ' 22 ', ord("8"): ' 23 ', ord("9"): ' 24 ',
                                   ord(":"): ' 25 ', ord(";"): ' 26 ', ord("<"): ' 27 ', ord("="): ' 28 ', ord(">"): ' 29 ',
                                   ord("?"): ' 30 ', ord("@"): ' 31 ', ord("A"): ' 32 ', ord("B"): ' 33 ', ord("C"): ' 34 ',
                                   ord("D"): ' 35 ', ord("E"): ' 36 ', ord("F"): ' 37 ', ord("G"): ' 38 ', ord("H"): ' 39 ',
                                   ord("I"): ' 40 '})
    final_score = [int(x) for x in re.findall(r"(?!'')[0-9][0-9]|[0-9](?=" ")", final_score)]
    return final_score

def clean_open_file(input_file):
    """Opens neatly and separately the fastq file as an iterator."""
    try:
        fastq_file = gzip.open(input_file, "rt")
        return fastq_file
    except (SystemExit, KeyboardInterrupt, IOError, FileNotFoundError):
        logs.error('The input file does not exist.', exc_info=True)
        sys.exit(1)



def Phred_score_statistics_calculation(input_file, sample_size, calculation_mode):
    """Calculates the mean or the absolute numeric Phred scores per each base using a desired sample size for random sampling from the fastq file."""
    if sample_size == None:
        cutoff = 10000000
    else:
        cutoff = int(sample_size)
    cycle_count = 0
    read_values = []
    Phred_T_sum, Phred_T_len, Phred_C_sum, Phred_C_len, Phred_A_sum, Phred_A_len, Phred_G_sum, Phred_G_len = [0] * 8
    for line in clean_open_file(input_file):
        if cycle_count == 1 or cycle_count % 4 == 1:
            read_sequence = re.findall(r"(?<!@)(?!@).*[^+|^\n].*(?!@\n)", line)[0]
        elif cycle_count == 3 or cycle_count % 4 == 3:
            read_base_qualities = re.findall(r"(?<!@)(?!@).*[^+|^\n].*(?!@\n)", line)
            read_base_qualities = numeric_Phred_score(read_base_qualities[0])
            results = [(i, *j) for i, j in
                       itertools.zip_longest(read_base_qualities, read_sequence, fillvalue='BLANK')]
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


def Phred_score_value_return(fastq_input, sample_size, calculation_mode, fq1, fq2):
    """Returns the numeric Phred scores per base for the forward and reverse reads in pair-end sequencing."""
    if fastq_input == fq1:
        read_values_fq1 = Phred_score_statistics_calculation(fastq_input, sample_size, calculation_mode)
        return read_values_fq1
    elif fastq_input == fq2:
        read_values_fq2 = Phred_score_statistics_calculation(fastq_input, sample_size, calculation_mode)
        return read_values_fq2


def main_Phred_score_calculation_output(fq1, fq2, sample_size, directory, name, calculation_mode):
    """Runs the Phred score calculation functions in parallel for the forward and the reverse reads in pair-end sequencing."""
    threads = []
    queue = [Queue(), Queue()]
    threads.append(Thread(target=lambda que, arg1, arg2, arg3, arg4, arg5: que.put(Phred_score_value_return(arg1, arg2, arg3, arg4, arg5)),
                          args=(queue[0], fq1, sample_size, calculation_mode, fq1, fq2), ))
    threads.append(Thread(target=lambda que, arg1, arg2, arg3, arg4, arg5: que.put(Phred_score_value_return(arg1, arg2, arg3, arg4, arg5)),
                          args=(queue[1], fq2, sample_size, calculation_mode, fq1, fq2), ))
    threads[0].start(), threads[1].start()
    for thread in threads:
        thread.join()
    read_values_fq1 = queue[0].get()
    read_values_fq2 = queue[1].get()
    return read_values_fq1, read_values_fq2


def summary_statistics_output(directory, name, statistics_data, read_orientation):
    """Outputs the Phred score calculation statistics as a text file that can be visualised independently from the plotting module."""
    with open(directory + name + '_total_Phred.txt', 'a', newline='') as new_file:
        data_line = csv.writer(new_file, delimiter='\t', lineterminator='\n')
        if read_orientation == 'F':
            read_orientation_string = 'First in pair_'
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


def Phred_values_return(input_values, read_orientation, directory, name, calculation_mode):
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
    summary_statistics_output(directory, name, statistics_data, read_orientation)
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


def Phred_scores_plotting(fq1, fq2, calculation_mode, directory, sample_size, minimum_score, colors, plot):
    """The general function that takes the input file, calculates mean or absolute Phred scores per base, and outputs
    their summary statistics as a text file or as a plot when the plotting module is enabled."""
    try:
        open(fq1, 'r')
        open(fq2, 'r')
    except (SystemExit, KeyboardInterrupt, IOError, FileNotFoundError):
        logs.error('The input fastq files do not exist.', exc_info=True)
        sys.exit(1)
    name = path.splitext(path.basename(fq1))[0]
    name = re.sub('_(R1|1).fq', '', name)
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    if colors != ['skyblue', 'mediumaquamarine', 'khaki', 'lightcoral']:
        colors = "".join(colors).split(',')
    read_values_fq1, read_values_fq2 = main_Phred_score_calculation_output(fq1, fq2, sample_size, directory, name, calculation_mode)
    data_fq1, maxy1 = Phred_values_return(read_values_fq1, 'F', directory, name, calculation_mode)
    data_fq2, maxy2 = Phred_values_return(read_values_fq2, 'R', directory, name, calculation_mode)
    if plot:
        pyp.figure()
        fig, fq = pyp.subplots(1, 2)
        fig.suptitle('Sequencing base quality', fontsize=14)
        pyp.subplots_adjust(wspace=0.4)
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
        pyp.vlines(-1, minimum_score, maxy, alpha=0.3, linewidth=1, linestyle='--', color='gray', clip_on=False)
        pyp.savefig(directory + name + '_phred_scores_plot.pdf', figsize=(14, 8), dpi=330, bbox_inches='tight')
        pyp.close()
    else:
        pass
    time_m = datetime.now()
    logs.info("asTair's Phred scores statistics summary function has finished running. {} seconds".format((
    time_m - time_b).total_seconds()))


if __name__ == '__main__':
    Phred_score_calculation_visualisation_exec()
