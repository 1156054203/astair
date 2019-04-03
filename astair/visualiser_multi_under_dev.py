#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import re
import os
import sys
import pdb
import csv
import gzip
import math
import numpy
import click
import random
import logging
import warnings
import itertools
import statistics
import subprocess
from os import path
from datetime import datetime
from collections import defaultdict
from multiprocessing import Process


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
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as mpatches
    pyp.style.use('seaborn-whitegrid')
    pyp.ioff()
except Exception:
    warnings.warn("Matplotlib was not found, visualisation output will not be supported.", ImportWarning)


class listsParamType(click.ParamType):
    name = 'lists'

    def convert(self, value, param, ctx):
        return value.split(',')

lists = listsParamType()

# from astair.safe_division import non_zero_division
# from astair.bam_file_parser import bam_file_opener
# from astair.simple_fasta_parser import fasta_splitting_by_sequence
# from astair.context_search import sequence_context_set_creation
# from astair.context_search import context_sequence_search


from simple_fasta_parser import fasta_splitting_by_sequence
from safe_division import non_zero_division
from context_search import sequence_context_set_creation
from context_search import context_sequence_search


### make to show coverage and modification per genome, per groups and per chromosomes #10 samples max; gives individual plots for mod level distribution and genomic mod line plots, metagenomic feature plots

@click.command()
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format used to find cytosine contexts.')
@click.option('input_file', '--input_file', '-i0', required=True, help='TAB delimited input MODS file.')
@click.option('input_file1', '--input_file1', '-i1', required=False, help='TAB delimited input MODS file.')
@click.option('input_file2', '--input_file2', '-i2', required=False, help='TAB delimited input MODS file.')
@click.option('input_file3', '--input_file3', '-i3', required=False, help='TAB delimited input MODS file.')
@click.option('input_file4', '--input_file4', '-i4', required=False, help='TAB delimited input MODS file.')
@click.option('input_file5', '--input_file5', '-i5', required=False, help='TAB delimited input MODS file.')
@click.option('input_file6', '--input_file6', '-i6', required=False, help='TAB delimited input MODS file.')
@click.option('input_file7', '--input_file7', '-i7', required=False, help='TAB delimited input MODS file.')
@click.option('input_file8', '--input_file8', '-i8', required=False, help='TAB delimited input MODS file.')
@click.option('input_file9', '--input_file9', '-i9', required=False, help='TAB delimited input MODS file.')
@click.option('features', '--features', '-fe', required=False, help='TAB delimited input BED file containing features of interest to make a metaplot.')
@click.option('metaplot_numbers', '--metaplot_numbers', '-mn', required=False, type=lists, help='Provide a list of distances and sizes if a metplot is desired, as a list positions with commas and no spaces e.g. 1000,100. The first number indicates the distance to the centre of the feature, and the second the number of equally sized bins to create.')
@click.option('mode', '--mode', '-mo', default='total', type=click.Choice(['total', 'strand']), help='Defines the coverage counting as total coverage (all bases covering a position), or strand that will only count C+T in OT and G+A in OB. (Default total.)')
@click.option('modification', '--mod', '-m', default=True, is_flag=True, required=False, help='Outputs modification calls distribution per context (Default True).')
@click.option('bin_size', '--bin', '-b', default=1000000, required=False, help='Non-overlapping bin size used for genomic line plots (Default 1 Mb).')
@click.option('maximum_coverage', '--maximum_coverage', '-md', required=False, default=50, type=int, help='The maximum coverage to be displayed on the x-axis. (Default 50).')
@click.option('normalisation', '--normalisation', '-n', type=click.Choice([None, 'GC', 'spikes', 'coverage', 'optimal']), help='Indicates what kind of normalisation to be used for group comparisons. (Default None).')
@click.option('context', '--context', '-co', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis.')
@click.option('boxplot', '--boxplot', '-bx', default=True, is_flag=True, help='Optional boxplot of the coverage distributions.')
@click.option('dmr', '--dmr', '-d', default=False, is_flag=True, help='Optional DMR comparison using an 2D discovery algorithm.')
@click.option('per_chromosome', '--per_chromosome', '-chr', default=None, type=str, help='When used, it visualises only per the chromosome given. (Default None')
@click.option('correlation', '--correlation', '-c', default='pearson', type=click.Choice(['pearson', 'spearman']), help='Specifies the similarity metric to use (Default Pearson correlation).')
@click.option('correlation_cutoff', '--correlation_cutoff', '-cc', default=3, type=int, help='Specifies the coverage cutoff used for the per base, per context correlation plots. (Default 3)')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
def visualiser_exec(reference, input_file, input_file1, input_file2, input_file3, input_file4, input_file5, input_file6, input_file7, input_file8, input_file9,
                    features, metaplot_numbers, mode, modification, bin_size, maximum_coverage, normalisation, context, dmr, per_chromosome, correlation, boxplot, correlation_cutoff, directory):
    visualiser(reference, input_file, input_file1, input_file2, input_file3, input_file4, input_file5, input_file6, input_file7, input_file8, input_file9,
                    features, metaplot_numbers, mode, modification, bin_size, maximum_coverage, normalisation, context, dmr, per_chromosome, correlation, boxplot, correlation_cutoff, directory)




#
# def normalisation(reference, input_file, input_file1, input_file2, input_file3, input_file4,
#                                   input_file5, input_file6, input_file7, input_file8, input_file9,
#                                   features, metaplot_numbers, mode, modification, bin_size, maximum_coverage,
#                                   normalisation, context, dmr, per_chromosome, correlation, boxplot, correlation_cutoff,
#                                   directory):
# ### should we have this as a separate module; including cnv and enrichment correction
#


def stats_counter(cyt_context, row, mode, maximum_depth, coverage, modification, total_length, modified_values):
    total_length[cyt_context] += 1
    if mode == 'total':
        if int(row[11]) <= maximum_depth:
            coverage[int(row[11])] += 1
        else:
            coverage[maximum_depth+1] += 1
    elif mode == 'strand':
        if int(row[4]) + int(row[5]) <= maximum_depth:
            coverage[int(row[4]) + int(row[5])] += 1
        else:
            coverage[maximum_depth+1] += 1
    if modification is True:
        modified_values[round(float(row[3]) * 100)] += 1


def reading_mods_files(inputs_available, maximum_depth, modification, mode, modified_values_CpG, modified_values_CHG, modified_values_CHH, covered_CpG, covered_CHG, covered_CHH, total_length):
    with open(inputs_available + '.mods') as lines:
        lines_in_file = csv.reader(lines, delimiter='\t', lineterminator='\n')
        for row in lines_in_file:
            if row[0] != 'CHROM':
                if mode == 'total' or (mode == 'strand' and row[4] + row[5] != 0):
                    if row[9] == 'CpG':
                        stats_counter('CpG', row, mode, maximum_depth, covered_CpG, modification, total_length, modified_values_CpG)
                    elif row[9] == 'CHG':
                        stats_counter('CHG', row, mode, maximum_depth, covered_CHG, modification, total_length, modified_values_CHG)
                    elif row[9] == 'CHH':
                        stats_counter('CHH', row, mode, maximum_depth, covered_CHH, modification, total_length, modified_values_CHH)



def preprocessor(reference, input_file, input_file1, input_file2, input_file3, input_file4, input_file5, input_file6, input_file7, input_file8, input_file9,
                    features, metaplot_numbers, mode, modification, bin_size, maximum_coverage, normalisation, context, dmr, per_chromosome, correlation, boxplot, correlation_cutoff, directory):
    try:
        keys, fastas = fasta_splitting_by_sequence(reference, per_chromosome)
    except Exception:
        sys.exit(1)
    contexts, all_keys = sequence_context_set_creation(context, None)
    context_total_counts = defaultdict(int)
    if per_chromosome == None:
        for i in range(0, len(keys)):
            time_m = datetime.now()
            modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys[i], None, context_total_counts, None)
    else:
        time_m = datetime.now()
        modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys, None, context_total_counts, None)
    inputs_available, all_names = list(), list()
    for inputs in [input_file, input_file1, input_file2, input_file3, input_file4, input_file5, input_file6, input_file7, input_file8, input_file9]:
        if inputs:
            name = path.splitext(path.basename(inputs))[0]
            inputs_available.append(name)
    if context == 'all':
        contexts_to_use = ['CpG', 'CHG', 'CHH']
    else:
        contexts_to_use = [context]
    all_names_coverage, all_names_mods = list(), list()
    for i in range(0, len(inputs_available)):
        modified_values_CpG, modified_values_CHG, modified_values_CHH = dict(), dict(), dict()
        covered_CpG, covered_CHG, covered_CHH = dict(), dict(), dict()
        for modified_values in [modified_values_CpG, modified_values_CHG, modified_values_CHH]:
            for l in range(0, 101):
                modified_values[l] = 0
        for covered in [covered_CpG, covered_CHG, covered_CHH]:
            for k in range(0, maximum_coverage+2):
                covered[k] = 0
        total_length = {'CpG': 0, 'CHG': 0, 'CHH': 0}
        reading_mods_files(inputs_available[i], maximum_coverage, modification, mode, modified_values_CpG, modified_values_CHG, modified_values_CHH, covered_CpG, covered_CHG, covered_CHH, total_length)
        for contexts in contexts_to_use:
            if contexts != 'CpG':
                total_covered = context_total_counts[contexts] + context_total_counts[contexts+'b']
            else:
                total_covered = context_total_counts['CG'] + context_total_counts['CGb']
            if contexts == 'CpG':
                all_values = [(i, covered_CpG[i]) for i in range(1, maximum_coverage+2)]
                all_values.insert(0, (0, abs(total_covered-sum(list(covered_CpG.values())))))
            elif contexts == 'CHG':
                all_values = [(i, covered_CHG[i]) for i in range(1, maximum_coverage+2)]
                all_values.insert(0, (0, abs(total_covered - sum(list(covered_CHG.values())))))
            else:
                all_values = [(i, covered_CHH[i]) for i in range(1, maximum_coverage+2)]
                all_values.insert(0, (0, abs(total_covered - sum(list(covered_CHH.values())))))
            all_names_coverage.append('{}_{}.hist'.format(inputs_available[i], contexts))
            with open('{}_{}.hist'.format(inputs_available[i], contexts), 'w') as lines:
                empty_lines = csv.writer(lines, delimiter='\t', lineterminator='\n')
                for row in all_values:
                    empty_lines.writerow(row)
    if modification is True:
        for i in range(0, len(inputs_available)):
            for contexts in contexts_to_use:
                if contexts == 'CpG':
                    all_values = [(i, modified_values_CpG[i]) for i in range(0, 101)]
                elif contexts == 'CHG':
                    all_values = [(i, modified_values_CHG[i]) for i in range(0, 101)]
                else:
                    all_values = [(i, modified_values_CHH[i]) for i in range(0, 101)]
                all_names_mods.append('{}_{}_mods.hist'.format(inputs_available[i], contexts))
                with open('{}_{}_mods.hist'.format(inputs_available[i], contexts), 'w') as lines:
                    empty_lines = csv.writer(lines, delimiter='\t', lineterminator='\n')
                    for row in all_values:
                        empty_lines.writerow(row)
    return all_names_coverage, all_names_mods



#
#
#
# def text_file_writer(reference, input_file, input_file1, input_file2, input_file3, input_file4,
#                       input_file5, input_file6, input_file7, input_file8, input_file9,
#                       features, metaplot_numbers, mode, modification, bin_size, maximum_coverage,
#                       normalisation, context, dmr, per_chromosome, correlation, boxplot, correlation_cutoff,
#                       directory):
#
# ### should we have this as a separate module; including cnv and enrichment correction
#
#

def summary_file_opener(input, input_bedx, input_bedy, data):
    with open(input) as bedgraph:
        bedgraph_file = csv.reader(bedgraph, delimiter='\t', lineterminator='\n')
        for row in bedgraph_file:
            input_bedx.append(int(row[0]))
            input_bedy.append(int(row[1]))
            data.append(int(row[0]) * int(row[1]))


def quartile_finder(input_bedy, known_var):
    coef = 0
    for x in input_bedy:
        if coef < known_var:
            coef += x
        elif coef >= known_var:
            unknown_index = input_bedy.index(x) - 1
            return unknown_index



def histogram_summariser(input_bedx, input_bedy, data):
    if sum(input_bedy) % 2 == 0:
        middle = sum(input_bedy) / 2
    else:
        middle = (math.ceil(sum(input_bedy) / 2))
    quartile = (math.ceil(sum(input_bedy) / 4))
    index = quartile_finder(input_bedy, middle)
    index25 = quartile_finder(input_bedy, quartile)
    index75 = quartile_finder(input_bedy, quartile*3)
    if index75 is None:
        index75 = 0
    if sum(input_bedy) % 2 == 0:
        median = round((input_bedx[index] + input_bedx[index + 1]) / 2, 3)
    else:
        median = round(input_bedx[index], 3)
    mu = round(sum(data) / sum(input_bedy), 3)
    sigma = numpy.std(input_bedx)
    quartile25 = round(input_bedx[index25], 3)
    quartile75 = round(input_bedx[index75], 3)
    proportions = [round((x / sum(input_bedy)) * 100, 3) for x in input_bedy]
    return mu, sigma, quartile25, quartile75, median, proportions


def plots(data, colort, cycle, all_names):
    cycle = cycle
    for box in data['boxes']:
        box.set(color=colort[cycle], linewidth=3)
        box.set(facecolor=colort[cycle], alpha=0.7)
        cycle += 1
    for j in range(0, len(all_names)):
        for element in ['whiskers', 'means', 'medians', 'caps']:
            pyp.setp(data[element][int(j*(non_zero_division(len(data[element]), len(all_names)))):int(j*(non_zero_division(len(data[element]),
                                                                                                                           len(all_names))) + non_zero_division(len(data[element]),
                                                                                                                                                                len(all_names)))], color=colort[j], linewidth=3)


def individual_histogram(input_bedx, proportions, maximum_coverage, name, mu, median, sigma, directory, type):
    pyp.figure()
    pyp.bar(input_bedx, proportions, width=1, facecolor='teal', alpha=0.8)
    props = dict(boxstyle='round', facecolor='white', alpha=0.1)
    if type == 'coverage':
        pyp.xlabel('Coverage', fontsize=12)
        pyp.ylabel('Frequency, %', fontsize=12)
        pyp.title('Histogram of coverage', fontsize=14)
        pyp.axis([0, maximum_coverage, 0, max(proportions) + 0.01 * max(proportions)])
        pyp.xticks(numpy.arange(0, maximum_coverage + 1, step=5))
    else:
        pyp.xlabel('Modification', fontsize=12)
        pyp.ylabel('Frequency, %', fontsize=12)
        pyp.title('Histogram of modification levels', fontsize=14)
        pyp.axis([0, 100, 0, max(proportions) + 0.01 * max(proportions)])
        pyp.xticks(numpy.arange(0, 100 + 1, step=5))
        pyp.xticks(fontsize=8)
        pyp.yticks(fontsize=8)
    pyp.grid(color='lightgray', linestyle='solid', linewidth=1)
    pyp.text(10, max(proportions) - 0.2, '$\mu=%.2f$\n$\mathrm{median}=%.2f$\n$\sigma=%.2f$' % (mu, median, sigma),
             fontsize=8, verticalalignment='top', bbox=props)
    pyp.savefig(directory + name + '_summary_plot.pdf', figsize=(14, 10), dpi=330,
                bbox_inches='tight')
    pyp.close()


def multiple_histograms(colors, all_names, stats, maximum_coverage, plot_specs, plotting_information, input_bedx_all, proportions, directory, name, contexts, type):
    pyp.figure()
    pyp.rcParams.update({'font.size': 18})
    props = dict(boxstyle='round', facecolor='white', alpha=0.1)
    if type == 'coverage':
        fig, axes = pyp.subplots(nrows=2, ncols=1, figsize=(20, 10), gridspec_kw={'height_ratios': [10, 2]})
        for i in range(0, len(all_names)):
            axes[0].bar(input_bedx_all[i], proportions[i], width=1, facecolor=colors[i], alpha=0.7, label=all_names[i][:-9])
        pyp.subplots_adjust(wspace=0, hspace=0)
        if stats is not None:
            positions = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]
            bp0 = axes[1].bxp(stats, positions=positions[:len(all_names)], patch_artist=True, vert=False, showfliers=False, widths=0.3)
            plots(bp0, colors, 0, all_names)
            axes[1].axis('off')
            axes[1].grid(False)
            axes[1].set_position([0.09, 0.78, 0.81, 0.1])
        axes[0].axis([0, maximum_coverage, 0, max(plot_specs) + 0.2 * max(plot_specs)])
        xaxis = [x if x != maximum_coverage else str(maximum_coverage) + '+' for x in
                 numpy.arange(0, maximum_coverage + 1, step=5)]
        axes[0].set_xticks(numpy.arange(0, maximum_coverage + 1, step=5))
        axes[0].ticklabel_format(style='sci', axis='y', scilimits=(0, 3))
        axes[0].set_xticklabels(xaxis)
        axes[0].set_xlabel('Coverage per position', fontsize=18)
        axes[0].set_ylabel('Number of positions', fontsize=18)
        pyp.legend(handles=[item for item in plotting_information.values()], bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.1)
    else:
        for i in range(0, len(all_names)):
            pyp.bar(input_bedx_all[i], proportions[i], facecolor=colors[i], alpha=0.7, label=all_names[i][:-9], width=1)
        pyp.axis([0, 100, 0, max(plot_specs) + 0.2 * max(plot_specs)])
        pyp.xticks(numpy.arange(0, 101, step=5))
        pyp.xticks(fontsize=8)
        pyp.yticks(fontsize=8)
        pyp.xlabel('Modification level per position', fontsize=12)
        pyp.ylabel('Number of positions', fontsize=12)
        pyp.legend(handles=[item for item in plotting_information.values()], bbox_to_anchor=(1.02, 1), loc=2, prop={'size': 6}, borderaxespad=0.1)
    pyp.grid(color='lightgray', linestyle='solid', linewidth=1)
    pyp.savefig(directory + name + '_' + contexts + '_complete_plot_comparison.pdf', figsize=(24, 10), dpi=330,
                bbox_inches='tight')
    pyp.close()


def coverager(reference, input_file, input_file1, input_file2, input_file3, input_file4, input_file5, input_file6, input_file7, input_file8, input_file9,
                    features, metaplot_numbers, mode, modification, bin_size, maximum_coverage, normalisation, context, dmr, per_chromosome, correlation, boxplot, correlation_cutoff, directory):

    proportions, plotting_information = dict(), dict()
    input_bedx_all, input_bedy_all, stats = dict(), dict(), dict()

    if context == 'all':
        contexts_to_use = ['CpG', 'CHG', 'CHH']
    else:
        contexts_to_use = [context]

    for contexts in contexts_to_use:
        input_bedx_all[contexts] = list()
        input_bedy_all[contexts] = list()
        stats[contexts] = list()

    all_names, all_names_mods = preprocessor(reference, input_file, input_file1, input_file2, input_file3, input_file4, input_file5, input_file6,
                 input_file7, input_file8, input_file9, features, metaplot_numbers, mode, modification, bin_size, maximum_coverage, normalisation, context,
                 dmr, per_chromosome, correlation, boxplot, correlation_cutoff, directory)

    for names in all_names:
        input_bedx, input_bedy, data = list(), list(), list()
        summary_file_opener(names, input_bedx, input_bedy, data)
        mu, sigma, quartile25, quartile75, median, proportions = histogram_summariser(input_bedx, input_bedy, data)
        individual_histogram(input_bedx, proportions, maximum_coverage, names, mu, median, sigma, directory, 'coverage')
        for contexts in contexts_to_use:
            if names[-8:-5] == contexts:
                if boxplot is not None:
                    for i in range(0, len(all_names)):
                        item_i = dict()
                        item_i["med"] = median
                        item_i["q1"] = quartile25
                        item_i["q3"] = quartile75
                        item_i["mean"] = mu
                        item_i['iqr'] = quartile75 - quartile25
                        item_i["whislo"] = min(input_bedx)
                        item_i["whishi"] = (quartile75 - quartile25) + 1.5 * quartile75
                        if item_i["whishi"] >= maximum_coverage:
                            item_i["whishi"] = maximum_coverage + 1
                    stats[contexts].append(item_i)
                else:
                    stats[contexts] = None
                input_bedx_all[contexts].append(input_bedx)
                input_bedy_all[contexts].append(input_bedy)

    colors = ['silver', 'lightseagreen', 'slateblue', 'orange', 'gray', 'coral', 'navy', 'deepskyblue', 'purple', 'green']

    proportions = dict()
    for contexts in contexts_to_use:
        context_names = [i[:-9] for i in all_names if i[-8:-5] == contexts]
        proportions[contexts] = list()
        for element in range(0, len(input_bedy_all[contexts])):
            if len(input_bedy_all[contexts][element]) > maximum_coverage:
                input_bedy_all[contexts][element][maximum_coverage] = sum(input_bedy_all[contexts][element][maximum_coverage:-1])
            proportions[contexts].append(input_bedy_all[contexts][element])
        plot_specs = max(input_bedy_all[contexts])
        for i in range(0, len(context_names)):
            plotting_information[i] = mpatches.Patch(color=colors[i], label=context_names[i], alpha=0.7)
        multiple_histograms(colors, context_names, stats[contexts], maximum_coverage, plot_specs, plotting_information, input_bedx_all[contexts], proportions[contexts], directory, all_names[0], contexts, 'coverage')
    return all_names_mods


def modification_distribution(maximum_coverage, normalisation, context, per_chromosome, directory, all_names_mods):

    proportions, plotting_information = dict(), dict()
    input_bedx_all, input_bedy_all = dict(), dict()

    if context == 'all':
        contexts_to_use = ['CpG', 'CHG', 'CHH']
    else:
        contexts_to_use = [context]

    for contexts in contexts_to_use:
        input_bedx_all[contexts] = list()
        input_bedy_all[contexts] = list()

    for names in all_names_mods:
        input_bedx, input_bedy, data = list(), list(), list()
        summary_file_opener(names, input_bedx, input_bedy, data)
        mu, sigma, quartile25, quartile75, median, proportions = histogram_summariser(input_bedx, input_bedy, data)
        individual_histogram(input_bedx, proportions, maximum_coverage, names, mu, median, sigma, directory, 'modification')
        for contexts in contexts_to_use:
            input_bedx_all[contexts].append(input_bedx)
            input_bedy_all[contexts].append(input_bedy)


    colors = ['silver', 'lightseagreen', 'slateblue', 'orange', 'gray', 'coral', 'navy', 'deepskyblue', 'purple', 'green']

    proportions = dict()
    for contexts in contexts_to_use:
        context_names = [i[:-10] for i in all_names_mods if i[-13:-10] == contexts]
        proportions[contexts] = list()
        for element in range(0, len(input_bedy_all[contexts])):
            if len(input_bedy_all[contexts][element]) > maximum_coverage:
                input_bedy_all[contexts][element][maximum_coverage] = sum(input_bedy_all[contexts][element][maximum_coverage:-1])
            proportions[contexts].append(input_bedy_all[contexts][element])
        plot_specs = max(input_bedy_all[contexts])
        for i in range(0, len(context_names)):
            plotting_information[i] = mpatches.Patch(color=colors[i], label=context_names[i], alpha=0.7)
        multiple_histograms(colors, context_names, None, maximum_coverage, plot_specs, plotting_information, input_bedx_all[contexts], proportions[contexts], directory, all_names_mods[0], contexts, 'modification')



#
#
# def modification_genomic_line_plot(reference, input_file, input_file1, input_file2, input_file3, input_file4,
#                                   input_file5, input_file6, input_file7, input_file8, input_file9,
#                                   features, metaplot_numbers, mode, modification, bin_size, maximum_coverage,
#                                   normalisation, context, dmr, per_chromosome, correlation, boxplot, correlation_cutoff,
#                                   directory):
# ### check how many input files are given and calculate coverage; if more options are specified run them meanwhile
#
#
# def metageneplot(reference, input_file, input_file1, input_file2, input_file3, input_file4,
#                                    input_file5, input_file6, input_file7, input_file8, input_file9,
#                                    features, metaplot_numbers, mode, modification, bin_size, maximum_coverage,
#                                    normalisation, context, dmr, per_chromosome, correlation, boxplot,
#                                    correlation_cutoff,
#                                    directory):
# ### check how many input files are given and calculate coverage; if more options are specified run them meanwhile
#
#
# def dmr(reference, input_file, input_file1, input_file2, input_file3, input_file4,
#                  input_file5, input_file6, input_file7, input_file8, input_file9,
#                  features, metaplot_numbers, mode, modification, bin_size, maximum_coverage,
#                  normalisation, context, dmr, per_chromosome, correlation, boxplot,
#                  correlation_cutoff,
#                  directory):
#
# ### check how many input files are given and calculate coverage; if more options are specified run them meanwhile


def visualiser(reference, input_file, input_file1, input_file2, input_file3, input_file4, input_file5, input_file6, input_file7, input_file8, input_file9,
                    features, metaplot_numbers, mode, modification, bin_size, maximum_coverage, normalisation, context, dmr, per_chromosome, correlation, boxplot, correlation_cutoff, directory):
    ### check how many input files are given and calculate coverage; if more options are specified run them meanwhile

    ######### MAIN calls others
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"

    # coverage histograms
    all_names_mods = coverager(reference, input_file, input_file1, input_file2, input_file3, input_file4, input_file5, input_file6,
              input_file7, input_file8, input_file9, features, metaplot_numbers, mode, modification, bin_size, maximum_coverage, normalisation, context, dmr,
              per_chromosome, correlation, boxplot, correlation_cutoff, directory)


    # modification histograms
    modification_distribution(maximum_coverage, normalisation, context, per_chromosome,directory, all_names_mods)




if __name__ == '__main__':
    visualiser_exec()
