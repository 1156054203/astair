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
import logging
import warnings
import itertools
import subprocess
from os import path
from datetime import datetime
from multiprocessing import Process


if sys.version[0] == '3':
    from itertools import zip_longest
elif sys.version[0] == '2':
    from itertools import izip_longest as zip_longest
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))

from astair.safe_division import non_zero_division
from astair.bam_file_parser import bam_file_opener
from astair.simple_fasta_parser import fasta_splitting_by_sequence
from astair.context_search import sequence_context_set_creation
from astair.context_search import context_sequence_search


@click.command()
@click.option('input_file', '--input_file', '-i', required=True, help='BAM|CRAM format file containing sequencing reads.')
@click.option('method', '--method', '-m', required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like).')
@click.option('bases_noncpg', '--bases_noncpg', default=3, type=int, help='The number of cytosines conversion events in CpH content to consider the read for removal. Default value is 3.')
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
def filter3T_exec(input_file, method, bases_noncpg, N_threads, directory):
    """Looks for sequencing reads with more than N CpH modifications."""
    removing_mod_err(input_file, method, bases_noncpg, N_threads, directory)

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

# logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()


def removing_mod_err(input_file, method, bases_noncpg, N_threads, directory):
    """ The main function to look for sequencing reads with more than N CpH modifications and remove them."""
    time_s = datetime.now()
    logs.info("asTair's excessive non-CpG modification read removal function started running. {} seconds".format((time_s - time_b).total_seconds()))
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    try:
        inbam = bam_file_opener(input_file, None, N_threads)
    except Exception:
        sys.exit(1)
    outbam3T = pysam.AlignmentFile(directory+name+"_high_CpH_filtered" + ".bam", "wb", template=inbam)
    removed3T = pysam.AlignmentFile(directory+name+"_high_CpH_removed" + ".bam", "wb", template=inbam)
    for read in bam_file_opener(input_file, 'fetch', N_threads):
        if read.flag == 147 or read.flag == 99:
            regs = "(?:.*C.*)" * int(bases_noncpg)
            ref = 'C'
        elif read.flag == 163 or read.flag == 83:
            regs = "(?:.*G.*)" * int(bases_noncpg)
            ref = 'G'
        try:
            read_data = read.get_tag('MD')
            if re.search(regs, read_data,re.IGNORECASE):
                changes = [int(s) for s in re.findall(r'\d+', read_data)]
                non_overlap = [x + 1 if x == 0 else x for x in changes]
                names = list(re.findall(r'[^\W\d_]+', read_data))
                positions = [x+2 for x in [sum(non_overlap[0:i]) for i in range(1, len(non_overlap)+1)]]
                positions[0] = changes[0]
                sequence = list(read.query_sequence)
                for k in range(0, len(names)):
                    element = positions[k]
                    sequence[element] = names[k]
                reads = "".join(sequence)
                if read.flag == 147 or read.flag == 99:
                    cpg = [m.start() for m in re.finditer(r'CG', reads, re.IGNORECASE)]
                elif read.flag == 163 or read.flag == 83:
                    cpg = [m.start() + 1 for m in re.finditer(r'CG', reads, re.IGNORECASE)]
                index = [(i, j) for i, j in zip_longest(positions, names, fillvalue='BLANK')]
                if method == 'mCtoT':
                    mismatch = [x[0] for x in index if x[1] == ref and x[0] not in cpg]
                else:
                    if read.flag == 147 or read.flag == 99:
                        total = [m.start() for m in re.finditer(r'C', reads, re.IGNORECASE)]
                    elif read.flag == 163 or read.flag == 83:
                        total = [m.start() for m in re.finditer(r'G', reads, re.IGNORECASE)]
                    non_cpg = [x for x in total if x not in cpg]
                    mismatch = [x for x in non_cpg if x not in positions]
                if len(mismatch) >= int(bases_noncpg):
                    removed3T.write(read)
                else:
                    outbam3T.write(read)
        except Exception:
            logs.error('The input file does not contain a MD tag column.', exc_info=True)
            sys.exit(1)
    inbam.close()
    time_m = datetime.now()
    logs.info("asTair's excessive non-CpG modification read removal function finished running. {} seconds".format((time_m - time_b).total_seconds()))


if __name__ == '__main__':
    filter3T_exec()







  
