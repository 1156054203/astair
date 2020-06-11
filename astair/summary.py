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
import numpy
import click
import pysam
import pickle
import shutil
import logging
import warnings
import subprocess
from os import path
from datetime import datetime
from collections import defaultdict


from astair.vcf_reader import read_vcf
from astair.safe_division import non_zero_division
from astair.bam_file_parser import bam_file_opener
from astair.context_search import context_sequence_search
from astair.context_search import sequence_context_set_creation
from astair.simple_fasta_parser import fasta_splitting_by_sequence


@click.command()
@click.option('input_file', '--input_file', '-i', required=True, help='BAM|CRAM format file containing sequencing reads.')
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format used for aligning of the sequencing reads and for pileup.')
@click.option('known_snp', '--known_snp', '-ks', default=None, required=False, help='VCF format file containing genotyped WGS high quality variants or known common variants in VCF format (dbSNP, 1000 genomes, etc.).')
@click.option('context', '--context', '-co', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis. (Default all).')
@click.option('user_defined_context', '--user_defined_context', '-uc', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
@click.option('library', '--library', '-li', required=False, default='directional',  type=click.Choice(['directional']), help='Provides information for the library preparation protocol (Default directional).')
@click.option('method', '--method', '-m', required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT).')
@click.option('single_end', '--single_end', '-se', default=False, is_flag=True, required=False, help='Indicates single-end sequencing reads (Default False).')
@click.option('region', '--region', '-r', nargs=3, type=click.Tuple([str, int, int]), default=(None, None, None), required=False, help='The one-based genomic coordinates of the specific region of interest given in the form chromosome, start position, end position, e.g. chr1 100 2000.')
@click.option('minimum_base_quality', '--minimum_base_quality', '-bq', required=False, type=int, default=0, help='Set the minimum base quality for a read base to be used in the pileup (Default 20).')
@click.option('minimum_mapping_quality', '--minimum_mapping_quality', '-mq', required=False, type=int, default=0, help='Set the minimum mapping quality for a read to be used in the pileup (Default 0).')
@click.option('per_chromosome', '--per_chromosome', '-chr', default=None, type=str, help='When used, it calculates the modification rates only per the chromosome given. (Default None).')
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
@click.option('add_underscores', '--add_underscores', '-au', default=False, is_flag=True, required=False, help='Indicates outputting a new reference fasta file with added underscores in the sequence names that is afterwards used for calling. (Default False).')
@click.option('no_information', '--no_information', '-ni', default='0', type=click.Choice(['.', '0', '*', 'NA']), required=False, help='What symbol should be used for a value where no enough quantative information is used. (Default *).')
def summarise(input_file, reference, known_snp,  context, user_defined_context, library,  method, region, minimum_base_quality, minimum_mapping_quality,per_chromosome, N_threads, directory, single_end, add_underscores, no_information):
    """Collects and outputs modification information per read."""
    read_summariser(input_file, reference, known_snp, context, user_defined_context, library,  method, region, minimum_base_quality, minimum_mapping_quality,per_chromosome, N_threads, directory,single_end, add_underscores, no_information)


warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)


logs = logging.getLogger(__name__)

time_b = datetime.now()

def read_summariser(input_file, reference, known_snp, context, user_defined_context, library,  method, region, minimum_base_quality, minimum_mapping_quality, per_chromosome, N_threads, directory, single_end, add_underscores, no_information):
    """Looks for cytosine contexts and their modification status and outpust read information."""
    time_s = datetime.now()
    logs.info("asTair's read information summary function started running. {} seconds".format((time_s - time_b).total_seconds()))
    name = path.splitext(path.basename(input_file))[0]
    if list(directory)[-1]!="/":
        directory = directory + "/"
    if path.exists(directory) == False:
        raise Exception("The output directory does not exist.")
        sys.exit(1)
    if per_chromosome is None:
        file_name = path.join(directory, name + "_" + method + "_" + context + "_read_summary.txt.gz")
    else:
        file_name = path.join(directory, name + "_" + method + "_" + per_chromosome + "_" + context + "read_summary.txt.gz")
    if no_information == '0':
        no_information = int(no_information)
    if not os.path.isfile(file_name):
        try:
            if per_chromosome==None:
                keys = fasta_splitting_by_sequence(reference, 'keys_only', None, add_underscores, None)
            else:
                keys = [per_chromosome]
            all_chrom = pysam.FastaFile(reference)
        except Exception:
            sys.exit(1)
        if single_end:
            flags_expectation_top, flags_expectation_bottom = [0], [16]
        else:
            flags_expectation_top, flags_expectation_bottom = [99,147], [83, 163]
        contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
        context_total_counts, true_variants, possible_mods = defaultdict(int), [], []
        if sys.version[0] == '3':
            data_line = gzip.open(file_name, 'wt', compresslevel=9, encoding='utf8', newline='\n')
        else:
            data_line = gzip.open(file_name, 'wt', compresslevel=9)    
        data_line.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format("#CHROM", "START", "END", "READ_NAME", "MODIFICATION_STATUS", 'FLAG', 'CONTEXT', "SPECIFIC_CONTEXT", "BQ", "MAPQ", "STRAND", "FRAGMENT_LENGTH", "AS", "XS", "EDIT", "KNOWN_SNP", "INFO"))
        if region != (None, None, None):
            keys, start, end = [region[0]], region[1], region[2]
        for i in range(0, len(keys)):
            time_m = datetime.now()
            logs.info("Starting read information extraction on {} chromosome (sequence). {} seconds".format(keys[i], (time_m - time_b).total_seconds()))
            if i == len(keys)-1:
                numbered = 'last'
            else:
                numbered = None
            fastas = fasta_splitting_by_sequence(reference, keys[i], numbered, add_underscores, None)
            modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys[i], user_defined_context, context_total_counts, None, None)
            if region == (None, None, None):
                start, end = 0, all_chrom.get_reference_length(keys[i])
            if known_snp is not None:
                time_s = datetime.now()
                logs.info("Starting reading SNP information on {} chromosome (sequence). {} seconds".format(keys[i], (time_s - time_m).total_seconds()))
                true_variants, possible_mods = read_vcf(known_snp, keys[i], fastas[keys[i]], N_threads, start, end)
                time_sf = datetime.now()
                logs.info("Reading SNP information on {} chromosome (sequence) has finished. {} seconds".format(keys[i], (time_sf - time_s).total_seconds()))
            try:
                inbam = bam_file_opener(input_file, (keys[i], start, end), N_threads)
            except Exception:
                sys.exit(1)
            for read in inbam:
                read_info = [(read.reference_name, i[1], i[1]+1) for i in read.get_aligned_pairs() if (i[1]!=None and i[0]!=None) and  (read.reference_name, i[1], i[1]+1) in modification_information_per_position]
                indices = [i[0] for i in read.get_aligned_pairs() if (i[1]!=None and i[0]!=None) and (read.reference_name, i[1], i[1]+1) in modification_information_per_position]
                if len(possible_mods)>0:
                    read_info_ = [(read.reference_name, i[1], i[1]+1) for i in read.get_aligned_pairs() if (i[1]!=None and i[0]!=None) and (read.reference_name, i[1], i[1]+1) in possible_mods]
                    indices_ = [i[0] for i in read.get_aligned_pairs() if (i[1]!=None and i[0]!=None) and (read.reference_name, i[1], i[1]+1) in possible_mods]
                    read_info.extend(read_info_)
                    indices.extend(indices_)
                snp_status = '*'
                if len(indices) > 0:
                    for k in range(0,len(read_info)):
                        position = read_info[k]
                        index = indices[k]
                        if (known_snp==None and position in modification_information_per_position) or (known_snp!=None and (position in modification_information_per_position or position in possible_mods)):
                            if known_snp!=None and position in possible_mods:
                                label = 'possible_modification'
                            else:
                                label = 'true_modification'
                            if  position in modification_information_per_position:
                                if modification_information_per_position[position][2] == 'C':
                                    strand = '+'
                                elif modification_information_per_position[position][2] == 'G':
                                    strand = '-'
                            elif position in possible_mods:
                                if possible_mods[position][2] == 'C':
                                    strand = '+'
                                elif possible_mods[position][2] == 'G':
                                    strand = '-'
                            if known_snp is not None and position in true_variants:
                                snp_status = 'WGS'
                            if position not in possible_mods:
                                if (read.flag in flags_expectation_top and read.query_sequence[index] in ['T','t'] and method == 'mCtoT') or  (read.flag in flags_expectation_top and read.query_sequence[index] in ['C','c'] and method == 'CtoT') or (read.flag in flags_expectation_bottom and read.query_sequence[index] in ['A','a'] and method == 'mCtoT') or  (read.flag in flags_expectation_bottom and read.query_sequence[index] in ['G','g'] and method == 'CtoT'):
                                    data_line.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(read.reference_name, position[1], position[2], read.query_name, 1, read.flag, modification_information_per_position[position][1], modification_information_per_position[position][0], read.query_qualities[index], read.mapq, strand, abs(read.template_length), read.get_tag('AS'), read.get_tag('XS'), read.get_tag('NM'), snp_status, label))
                                else:
                                    data_line.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(read.reference_name, position[1], position[2], read.query_name, 0, read.flag, modification_information_per_position[position][1], modification_information_per_position[position][0], read.query_qualities[index], read.mapq, strand, abs(read.template_length), read.get_tag('AS'), read.get_tag('XS'), read.get_tag('NM'), snp_status, label))
                            else:
                                if (read.flag in flags_expectation_top and read.query_sequence[index] in ['T','t'] and method == 'mCtoT') or  (read.flag in flags_expectation_top and read.query_sequence[index] in ['C','c'] and method == 'CtoT') or (read.flag in flags_expectation_bottom and read.query_sequence[index] in ['A','a'] and method == 'mCtoT') or  (read.flag in flags_expectation_bottom and read.query_sequence[index] in ['G','g'] and method == 'CtoT'):
                                    data_line.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(read.reference_name, position[1], position[2], read.query_name, 1, read.flag, possible_mods[position][1], possible_mods[position][0], read.query_qualities[index], read.mapq, strand, abs(read.template_length), read.get_tag('AS'), read.get_tag('XS'), read.get_tag('NM'), snp_status, label))
                                else:
                                    data_line.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(read.reference_name, position[1], position[2], read.query_name, 0, read.flag, possible_mods[position][1], possible_mods[position][0], read.query_qualities[index], read.mapq, strand, abs(read.template_length),read.get_tag('AS'), read.get_tag('XS'), read.get_tag('NM'), snp_status, label))
            modification_information_per_position, fastas = None, None
        data_line.close()
        all_chrom.close()
        time_e = datetime.now()
        logs.info("asTair read information summary function finished running. {} seconds".format((time_e - time_b).total_seconds()))
    else:
        logs.error('Read information summary file with this name exists. Please rename before rerunning.')
        sys.exit(1)
        
        

if __name__ == '__main__':
    summarise()
        
