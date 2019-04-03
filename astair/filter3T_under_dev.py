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
    removing_mod_err(input_file, method, bases_noncpg, N_threads, directory)

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

# logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()


def removing_mod_err(input_file, method, bases_noncpg, N_threads, directory):
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    try:
        inbam = bam_file_opener(input_file, 'fetch', N_threads)
    except Exception:
        sys.exit(1)
    outbam3T = pysam.AlignmentFile(name+"3T_"+".bam", "wb", template=inbam)
    removed3T = pysam.AlignmentFile(name+"removed_" + ".bam", "wb", template=inbam)
    for read in inbam(input_file, 'fetch', N_threads):
        if read.flag == 83 or read.flag == 99:
            regs = "(?:.*C.*)" * int(bases_noncpg)
            ch = "(C)(A|C|T)(A|T|C|G)"
            ref = 'C'
        elif read.flag == 163 or read.flag == 147:
            regs = "(?:.*G.*)" * int(bases_noncpg)
            ch = "(G)(A|G|T)(A|T|C|G)"
            ref = 'G'
        if isinstance(read.tags[0][1],str):
            read_data = read.tags[0][1]
        elif isinstance(read.tags[1][1],str):
            read_data = read.tags[1][1]
            if re.search(regs, read_data,re.IGNORECASE):
                changes = [int(s) for s in re.findall(r'\d+', read_data)]
                non_overlap = [x + 1 if x == 0 else x for x in changes]
                names = list(re.findall(r'[^\W\d_]+', read_data))
                positions = [x +1 for x in list(itertools.accumulate(non_overlap))]
                positions = positions[:-1]
                sequence = list(read.query_sequence)
                for k in range(len(positions)-1):
                    element = positions[k]
                    sequence[element] = names[k]
                reads = "".join(sequence)
                non_cpg = [m.start() for m in re.finditer(ch, reads, re.IGNORECASE)]
                index = [(i, j) for i, j in zip_longest(positions, names, fillvalue='BLANK')]
                if non_cpg:
                    mismatch = [x[0] for x in index if x[1] == ref]
                    if len(mismatch) >= int(bases_noncpg):
                       if mismatch not in non_cpg:
                            outbam3T.write(read)
                       else:
                            removed3T.write(read)
                    else:
                        outbam3T.write(read)
                else:
                    outbam3T.write(read)
            else:
                outbam3T.write(read)
    inbam.close()
    time_m = datetime.now()
    logs.info("asTair's excessive non-CpG modification read removal"
              " function has finished running. {} seconds".format((time_m - time_b).total_seconds()))



if __name__ == '__main__':
    filter3T_exec()







  
