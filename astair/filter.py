#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import re
import os
import sys
import pdb
import click
import pysam
import logging
import warnings
from os import path
from datetime import datetime


from astair.cigar_search import cigar_search
from astair.bam_file_parser import bam_file_opener
from astair.simple_fasta_parser import fasta_splitting_by_sequence


@click.command()
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format.')
@click.option('input_file', '--input_file', '-i', required=True, help='BAM|CRAM format file containing sequencing reads.')
@click.option('method', '--method', '-m', required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT).')
@click.option('bases_noncpg', '--bases_noncpg', default=3, type=int, help='The number of cytosines conversion events in CpH content to consider the read for removal. Default value is 3.')
@click.option('per_chromosome', '--per_chromosome', '-chr', default=None, type=str, help='When used, it calculates the modification rates only per the chromosome given. (Default None).')
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('single_end', '--se', '-se', default=False, is_flag=True, required=False, help='Indicates single-end sequencing reads (Default False).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
def filter(reference, input_file, method, bases_noncpg, per_chromosome, N_threads, single_end, directory):
    """Look for sequencing reads with more than N CpH modifications."""
    removing_mod_err(reference, input_file, method, bases_noncpg, per_chromosome, N_threads, single_end, directory)


warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

# logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

time_b = datetime.now()


def position_correction_cigar(read, positions):
    """Uses the CIGAR string information to correct the expected cytosine positions."""
    names, positions_cigar, changes = cigar_search(read.cigarstring)
    index = 0
    for change in names:
        if len(positions) != 0:
            if change == 'D':
                corrected_positions = [x if x < positions_cigar[index] else x - changes[index] for x in positions]
                index += 1
                positions = corrected_positions
            elif change == 'I':
                corrected_positions = [x if x < positions_cigar[index] else x + changes[index] for x in positions]
                index += 1
                positions = corrected_positions
            elif change == 'S' or change == 'H':
                if index == 0:
                    corrected_positions = [x for x in positions if x > positions_cigar[index]]
                else:
                    corrected_positions = [x for x in positions if x < positions_cigar[index]]
                index += 1
                positions = corrected_positions
            else:
                index += 1
    return positions


def removing_mod_err(reference, input_file, method, bases_noncpg, per_chromosome, N_threads, single_end, directory):
    """ The main function to look for sequencing reads with more than N CpH modifications and remove them."""
    time_s = datetime.now()
    logs.info("asTair's excessive non-CpG modification read removal function started running. {} seconds".format((time_s - time_b).total_seconds()))
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if list(directory)[-1]!="/":
        directory = directory + "/"
    if os.path.exists(directory) == False:
        raise Exception("The output directory does not exist.")
        sys.exit(1)
    try:
        inbam = bam_file_opener(input_file, None, N_threads)
    except Exception:
        sys.exit(1)
    keys, fastas = fasta_splitting_by_sequence(reference, per_chromosome, None, False, 'all')
    outbam3T = pysam.AlignmentFile(directory+name+"_high_CpH_filtered" + ".bam", "wb", template=inbam)
    removed3T = pysam.AlignmentFile(directory+name+"_high_CpH_removed" + ".bam", "wb", template=inbam)
    if single_end == False:
        read_flags = [83, 99, 147, 163]
    else:
        read_flags = [0, 16]
    for read in bam_file_opener(input_file, 'fetch', N_threads):
        if read.flag in read_flags:
            if read.flag == 147 or read.flag == 99 or read.flag == 0:
                regs = "(?:.*C.*)" * int(bases_noncpg)
                ref = 'C'
            elif read.flag == 163 or read.flag == 83 or read.flag == 16:
                regs = "(?:.*G.*)" * int(bases_noncpg)
                ref = 'G'
            try:
                read_data = read.get_tag('MD')
                if re.search(regs, read_data,re.IGNORECASE):
                    if read.flag == 147 or read.flag == 99 or read.flag == 0:
                        cpg = [m.start() for m in re.finditer(r'CG', fastas[read.reference_name][
                                                                    read.reference_start:read.reference_start + read.qlen],
                                                            re.IGNORECASE)]
                        ref, alt = 'C', 'T'
                    elif read.flag == 163 or read.flag == 83 or read.flag == 16:
                        cpg = [m.start() + 1 for m in re.finditer(r'CG', fastas[read.reference_name][
                                                                        read.reference_start:read.reference_start + read.qlen],
                                                                re.IGNORECASE)]
                        ref, alt = 'G', 'A'
                    if len(read.tags) != 0 and (re.findall('I', read.cigarstring, re.IGNORECASE) or re.findall('D', read.cigarstring, re.IGNORECASE)) or re.findall('S',read.cigarstring, re.IGNORECASE) or re.findall('H',
                                                                                                                    read.cigarstring,
                                                                                                                    re.IGNORECASE):
                        total = [m.start() for m in re.finditer(ref, fastas[read.reference_name][read.reference_start:read.reference_start+read.qlen], re.IGNORECASE)]
                        total_ref = position_correction_cigar(read, total)
                        cpg = position_correction_cigar(read, cpg)
                    else:
                        total_ref = [m.start() for m in re.finditer(ref, fastas[read.reference_name][read.reference_start:read.reference_start+read.qlen], re.IGNORECASE)]
                    if method == 'mCtoT':
                        mismatch = (set(total_ref).intersection(set([m.start() for m in re.finditer(alt, read.query_sequence, re.IGNORECASE)]))).difference(cpg)
                    else:
                        mismatch = (set(total_ref).intersection((set(total_ref).difference(cpg))).difference([m.start() for m in re.finditer(alt, read.query_sequence, re.IGNORECASE)]))
                    if len(mismatch) >= int(bases_noncpg):
                        removed3T.write(read)
                    else:
                        outbam3T.write(read)
                else:
                    outbam3T.write(read)
            except Exception:
                logs.error('The input file does not contain a MD tag column.', exc_info=True)
                sys.exit(1)
    inbam.close()
    outbam3T.close()
    removed3T.close()
    time_m = datetime.now()
    logs.info("asTair's excessive non-CpG modification read removal function finished running. {} seconds".format((time_m - time_b).total_seconds()))


if __name__ == '__main__':
    filter()

