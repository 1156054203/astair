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
import click
import pysam
import shutil
import logging
import warnings
import itertools
import subprocess
from os import path
from datetime import datetime
from multiprocessing import Process
from collections import defaultdict

if sys.version[0] == '3':
    from itertools import zip_longest
elif sys.version[0] == '2':
    from itertools import izip_longest as zip_longest
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))

from astair.safe_division import non_zero_division


class listsParamType(click.ParamType):
    name = 'lists'

    def convert(self, value, param, ctx):
        return value.split(',')

lists = listsParamType()


@click.command()
@click.option('input_file', '--input_file', '-i', required=True, help='BAM|CRAM format file containing sequencing reads.')
@click.option('read_length', '--read_length', '-l', type=int, required=True, help='The read length is needed to separate the reads.')
@click.option('method', '--method', '-m', required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT)')
@click.option('modified_positions', '--modified_positions', '-mp', type=lists, required=True, help='Provide a list of positions that were modified in one-based coordinates. List positions with commas and no spaces e.g. 51,111.')
@click.option('modified_positions_orientation', '--modified_positions_orientation', '-mpo', type=lists, required=True, help='Provide a list of additional information for the modified positions with commas and no spaces e.g. OT,OB')
@click.option('output_bam', '-output_bam', '-sb', required=False, default=False, is_flag=True, help='If given, bam files separated by the context surrounding the positions of interest will be output.')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
def separator(input_file, read_length, method, modified_positions, modified_positions_orientation, output_bam, directory):
    """Separates a bam file based on the context at certain positions.
    Suitable for highly variable and covered spike-ins, as it calculates a rough estimate of the modification value at
    few known positions."""
    position_separator(input_file, read_length, method, modified_positions, modified_positions_orientation, output_bam, directory)




def base_counter(method, modified_bases, unmodified_bases, position, context, read_info, base):
    """Does a simplified counting per context certain positions."""
    if method == 'mCtoT':
        if (base == 'T' and read_info == 'OT') or (base == 'A' and read_info == 'OB'):
            modified_bases[(position, context)] += 1
        elif (base == 'C' and read_info == 'OT') or (base == 'G' and read_info == 'OB'):
            unmodified_bases[(position, context)] += 1
    else:
        if (base == 'C' and read_info == 'OT') or (base == 'G' and read_info == 'OB'):
            modified_bases[(position, context)] += 1
        elif (base == 'T' and read_info == 'OT') or (base == 'A' and read_info == 'OB'):
            unmodified_bases[(position, context)] += 1


def mods_like_context_writer(header, file_name, modified_bases, unmodified_bases):
    """Outputs a mods file containing only information about the positions of interest."""
    sorted_keys_mod = [('spike', dict_keys[0]-1, dict_keys[0], dict_keys[1], dict_values) for dict_keys, dict_values in modified_bases.items()]
    sorted_keys_mod.sort()
    sorted_keys_unmod = [('spike', dict_keys[0]-1, dict_keys[0], dict_keys[1], dict_values)for dict_keys, dict_values in unmodified_bases.items()]
    sorted_keys_unmod.sort()
    mod_values = [(x[0], x[1], x[2], x[3], round(non_zero_division(x[4], (x[4]+y[4]))*100,3), x[4], y[4]) for x in sorted_keys_mod for y in sorted_keys_unmod if y[1]==x[1] and y[3]==x[3]]
    with open(file_name, 'w') as contexts_file:
        write_fasta = csv.writer(contexts_file, delimiter='\t', lineterminator='\n')
        if header == True:
            write_fasta.writerow(["CHROM", "START", "END", "CONTEXT", "MEAN_MODIFICATION_RATE_PERCENT", "MOD", "UNMOD"])
        for line in mod_values:
            write_fasta.writerow(line)


def per_context_and_query_output(output_bam, read, outbam_CHH, outbam_CHG, outbam_CpG, whole, start, lens, read_info, three_bases, base, method, modified_bases, unmodified_bases, index1, index2):
    """Gives a generalised output as mods and bam files if bam output is requested."""
    if read_info == "OT":
        strand = +1
    else:
        strand = -1
    if read.query_sequence[index1] in three_bases and read.query_sequence[index1 + strand] in three_bases and lens != 0:
        base_counter(method, modified_bases, unmodified_bases, abs(whole + start), 'CHH', read_info,
                     read.query_sequence[index2 - 1])
        if output_bam is True:
            outbam_CHH.write(read)
    elif read.query_sequence[index1] in three_bases and read.query_sequence[index1 + strand] == base and lens != 0:
        base_counter(method, modified_bases, unmodified_bases, abs(whole + start), 'CHG', read_info,
                     read.query_sequence[index2 - 1])
        if output_bam is True:
            outbam_CHG.write(read)
    elif read.query_sequence[index1] == base and lens != 0:
        base_counter(method, modified_bases, unmodified_bases, abs(whole + start), 'CpG', read_info,
                     read.query_sequence[index2 - 1])
        if output_bam is True:
            outbam_CpG.write(read)


def position_separator(input_file, read_length, method, modified_positions, modified_positions_orientation, output_bam, directory):
    """Does the main bam splitting by context."""
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    if os.path.exists(directory) == False:
        raise Exception("The output directory does not exist.")
        sys.exit(1)
    modified_information = dict([(int(x), y) for x,y in zip(modified_positions, modified_positions_orientation)])
    inbam = pysam.AlignmentFile(input_file, "rb", header=True)
    bam_fetch = inbam.fetch(until_eof=True)
    if output_bam is True:
        outbam_CpG = pysam.AlignmentFile(directory + name + '_' + method + "_CpG" + ".bam", "wb", header=True, template=inbam)
        outbam_CHH = pysam.AlignmentFile(directory + name + '_' + method + "_CHH" + ".bam", "wb", header=True, template=inbam)
        outbam_CHG = pysam.AlignmentFile(directory + name + '_' + method + "_CHG" + ".bam", "wb", header=True, template=inbam)
    else:
        outbam_CpG, outbam_CHG, outbam_CHH = None, None, None
    modified_bases, unmodified_bases = {}, {}
    for position in modified_positions:
        for context in ['CpG', 'CHG', 'CHH']:
            modified_bases[(int(position), context)] = 0
            unmodified_bases[(int(position), context)] = 0
    for read in bam_fetch:
        if len(read.tags) > 0:
            if isinstance(read.tags[0][1],str):
                read_data = read.tags[0][1]
            elif isinstance(read.tags[1][1],str):
                read_data = read.tags[1][1]
        else:
            read_data = list()
        start = read.reference_start
        end = read.reference_start + read.query_length
        if read.flag == 99 or read.flag == 147:
            read_info = 'OT'
        elif read.flag == 163 or read.flag == 83:
            read_info = 'OB'
        else:
            read_info = None
        if (len(read_data)==0 or not re.search("\^", read_data,re.IGNORECASE)) and len(read.query_sequence) <= read_length:
            for whole in modified_information.keys():
                if whole in range(start, end+1):
                    if (modified_information[whole] == 'OT' and whole in range(start, end+1) and read_info == 'OT') or (modified_information[whole] == 'OT' and whole in range(start, end+1) and read_info == 'OB'):
                        whole = abs(whole - start)
                        string_of_interest = read.query_sequence[whole-1:whole+2]
                        lens = len(string_of_interest)
                        if (read.qstart == 0 or read.qstart == read_length) and len(read.query_sequence) > whole + 1:
                            per_context_and_query_output(output_bam, read, outbam_CHH, outbam_CHG, outbam_CpG, whole, start, lens, read_info, ["A", 'C', "T"], 'G', method, modified_bases, unmodified_bases, whole, whole)
                        elif (read.qstart != 0 or read.qstart != read_length) and len(read.query_sequence) > abs(whole - read.qstart) + 1 > 0:
                             if read.qstart > 0 and read.qstart < read_length:
                                 per_context_and_query_output(output_bam, read, outbam_CHH, outbam_CHG, outbam_CpG, whole, start, lens, read_info, ["A", 'C', "T"], 'G', method, modified_bases, unmodified_bases, abs(whole - read.qstart), abs(whole - read.qstart))
                             elif read.qstart > read_length:
                                 per_context_and_query_output(output_bam, read, outbam_CHH, outbam_CHG, outbam_CpG, whole, start, lens, read_info, ["A", 'C', "T"], 'G', method, modified_bases, unmodified_bases, abs(whole - abs(read_length - read.qstart)), abs(whole - abs(read_length - read.qstart)))
                    elif (modified_information[whole] == 'OB' and whole in range(start, end+1) and read_info == 'OT') or (modified_information[whole] == 'OB' and whole in range(start, end+1) and read_info == 'OB'):
                        whole = abs(whole - start)
                        string_of_interest = read.query_sequence[whole-1:whole+2]
                        lens = len(string_of_interest)
                        if (read.qstart == 0 or read.qstart == read_length) and len(read.query_sequence) > whole + 1:
                            per_context_and_query_output(output_bam, read, outbam_CHH, outbam_CHG, outbam_CpG, whole, start, lens, read_info, ["A", 'G', "T"], 'C', method, modified_bases, unmodified_bases, whole-2, whole)
                        elif (read.qstart != 0 or read.qstart != read_length) and len(read.query_sequence) > abs(whole - read.qstart) + 1 > 0:
                            if read.qstart > 0 and read.qstart < read_length:
                                per_context_and_query_output(output_bam, read, outbam_CHH, outbam_CHG, outbam_CpG, whole, start, lens, read_info, ["A", 'G', "T"], 'C', method, modified_bases, unmodified_bases, abs(whole - read.qstart)-2, abs(whole - read.qstart))
                            elif read.qstart > read_length:
                                per_context_and_query_output(output_bam, read, outbam_CHH, outbam_CHG, outbam_CpG, whole, start, lens, read_info, ["A", 'G', "T"], 'C', method, modified_bases, unmodified_bases, abs(whole - abs(read_length - read.qstart))-2, abs(whole - abs(read_length - read.qstart)))
    file_name = path.join(directory, name + '_' + method + "_positions" + ".mods")
    mods_like_context_writer(True, file_name, modified_bases, unmodified_bases)


if __name__ == '__main__':
    separator()






  
