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
import subprocess
from os import path
from datetime import datetime
from collections import defaultdict

if sys.version[0] == '3':
    from itertools import zip_longest
elif sys.version[0] == '2':
    from itertools import izip_longest as zip_longest
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))


from astair.simple_fasta_parser import fasta_splitting_by_sequence
from astair.context_search import sequence_context_set_creation
from astair.context_search import context_sequence_search



@click.command()
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format used for aligning of the sequencing reads and for pileup.')
@click.option('context', '--context', '-co', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis.')
@click.option('user_defined_context', '--user_defined_context', '-uc', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
@click.option('per_chromosome', '--per_chromosome', '-chr', default=None, type=str, help='When used, it calculates the modification rates only per the chromosome given. (Default None')
@click.option('compress', '--gz', '-z', default=False, is_flag=True, required=False, help='Indicates whether the mods file output will be compressed with gzip (Default False).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
def find(reference, context, user_defined_context, per_chromosome, compress, directory):
    """Outputs positions of Cs from fasta file per context."""
    find_contexts(reference, context, user_defined_context, per_chromosome, compress, directory)



warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

# logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()


        
def bed_like_context_writer(header, file_name, cytosine_contexts):
    """Writing cytosine genomic coordinates in a bed-like format."""
    sorted_keys = [dict_keys for dict_keys in cytosine_contexts]
    sorted_keys.sort()
    with open(file_name, 'a') as fasta_contexts_file:
        write_fasta = csv.writer(fasta_contexts_file, delimiter='\t', lineterminator='\n')
        if header == True:
            write_fasta.writerow(["CHROM", "START", "END", "SPECIFIC_CONTEXT", "CONTEXT"])
        for line in sorted_keys:
            write_fasta.writerow([line[0], line[1], line[2], cytosine_contexts[line][0], cytosine_contexts[line][1]])



def find_contexts(reference, context, user_defined_context, per_chromosome, compress, directory):
    """Looks for the coordinates of cytosines in different contexts."""
    time_s = datetime.now()
    logs.info("asTair cytosine contexts positions finder started running. {} seconds".format((time_s - time_b).total_seconds()))
    name = path.splitext(path.basename(reference))[0]
    directory = path.abspath(directory)
    if per_chromosome == None:
        file_name = path.join(directory, name + "_"  + context + ".bed")
    else:
        file_name = path.join(directory, name + "_" + per_chromosome + "_" + context + ".bed")
    try:
        keys, fastas = fasta_splitting_by_sequence(reference, per_chromosome)
    except Exception:
        sys.exit(1)
    if not os.path.isfile(file_name) and not os.path.isfile(file_name + '.gz'):
        contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
        cycles = 0
        context_total_counts, context_sample_counts = defaultdict(int), defaultdict(int)
        if per_chromosome == None:
            for i in range(0, len(keys)):
                time_m = datetime.now()
                logs.info("Looking for cytosine positions on {} chromosome (sequence). {} seconds".format(keys[i], (time_m - time_b).total_seconds()))
                cytosine_contexts = context_sequence_search(contexts, all_keys, fastas, keys[i], user_defined_context, context_total_counts, None)
                if i == 0:
                    bed_like_context_writer(True, file_name, cytosine_contexts)
                else:
                    bed_like_context_writer(False, file_name, cytosine_contexts)
        else:
            cytosine_contexts = context_sequence_search(contexts, all_keys, fastas, keys, user_defined_context, context_total_counts, None)
            bed_like_context_writer(True, file_name, cytosine_contexts)
        if compress == True:
            logs.info("Compressing bed file output cytosine genomic positions per context.")
            subprocess.Popen('gzip -f -9 {}'.format(file_name), shell=True)
        time_e = datetime.now()
        logs.info("asTair cytosine contexts positions finder finished running. {} seconds".format((time_e - time_b).total_seconds()))
    else:
        logs.error('Bed file with this name exists. Please rename before rerunning.')
        sys.exit(1)




if __name__ == '__main__':
    find_contexts()
