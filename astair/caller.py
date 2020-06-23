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

if sys.version[0] == '3': 
    from itertools import zip_longest
elif sys.version[0] == '2':
    from itertools import izip_longest as zip_longest
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))


from astair.vcf_reader import read_vcf
from astair.safe_division import safe_rounder
from astair.safe_division import non_zero_division
from astair.bam_file_parser import bam_file_opener
from astair.context_search import context_sequence_search
from astair.context_search import sequence_context_set_creation
from astair.simple_fasta_parser import fasta_splitting_by_sequence

@click.command()
@click.option('input_file', '--input_file', '-i', required=True, help='BAM|CRAM format file containing sequencing reads.')
#@click.option('control_file', '--control_file', '-c', required=False, help='BAM|CRAM format file containing sequencing reads used as a matched control.')
@click.option('known_snp', '--known_snp', '-ks', required=False, help='VCF format file containing genotyped WGS high quality variants or known common variants in VCF format (dbSNP, 1000 genomes, etc.).')
@click.option('model', '--model', '-mo', default='none',  type=click.Choice(['none']), required=False, help='Decide on model for class estimation.')
@click.option('reference', '--reference', '-f', required=True, help='Reference DNA sequence in FASTA format used for aligning of the sequencing reads and for pileup.')
#@click.option('exclude_variants', '--exclude_variants', '-ev', default=False, is_flag=True, help='When set to true does not output variants.)
@click.option('zero_coverage', '--zero_coverage', '-zc', default=False, is_flag=True, help='When set to True, outputs positions not covered in the bam file. Uncovering zero coverage positions takes longer time than using the default option.')
@click.option('context', '--context', '-co', required=False, default='all',  type=click.Choice(['all', 'CpG', 'CHG', 'CHH']), help='Explains which cytosine sequence contexts are to be expected in the output file. Default behaviour is all, which includes CpG, CHG, CHH contexts and their sub-contexts for downstream filtering and analysis. (Default all).')
@click.option('user_defined_context', '--user_defined_context', '-uc', required=False, type=str, help='At least two-letter contexts other than CG, CHH and CHG to be evaluated, will return the genomic coordinates for the first cytosine in the string.')
@click.option('library', '--library', '-li', required=False, default='directional',  type=click.Choice(['directional', 'reverse']), help='Provides information for the library preparation protocol (Default directional).')
@click.option('method', '--method', '-m', required=False, default='mCtoT', type=click.Choice(['CtoT', 'mCtoT']), help='Specify sequencing method, possible options are CtoT (unmodified cytosines are converted to thymines, bisulfite sequencing-like) and mCtoT (modified cytosines are converted to thymines, TAPS-like). (Default mCtoT).')
@click.option('skip_clip_overlap', '--skip_clip_overlap', '-sc', required=False, default=False, type=bool, help='Skipping the random removal of overlapping bases between pair-end reads. Not recommended for pair-end libraries, unless the overlaps are removed prior to calling. (Default False)')
@click.option('single_end', '--single_end', '-se', default=False, is_flag=True, required=False, help='Indicates single-end sequencing reads (Default False).')
@click.option('minimum_base_quality', '--minimum_base_quality', '-bq', required=False, type=int, default=20, help='Set the minimum base quality for a read base to be used in the pileup (Default 20).')
@click.option('minimum_mapping_quality', '--minimum_mapping_quality', '-mq', required=False, type=int, default=0, help='Set the minimum mapping quality for a read to be used in the pileup (Default 0).')
@click.option('adjust_acapq_threshold', '--adjust_capq_threshold', '-amq', required=False, type=int, default=0, help='Used to adjust the mapping quality with default 0 for no adjustment and a recommended value for adjustment 50. (Default 0).')
@click.option('add_indels', '--add_indels', '-ai', required=False, default=True, type=bool, help='Adds inserted bases and Ns for base skipped from the reference (Default True).')
@click.option('redo_baq', '--redo_baq', '-rbq', required=False, default=False, type=bool, help='Re-calculates per-Base Alignment Qualities ignoring existing base qualities (Default False).')
@click.option('compute_baq', '--compute_baq', '-cbq', required=False, default=True, type=bool, help='Performs re-alignment computing of per-Base Alignment Qualities (Default True).')
@click.option('ignore_orphans', '--ignore_orphans', '-io', required=False, default=True, type=bool, help='Ignore reads not in proper pairs (Default True).')
@click.option('max_depth', '--max_depth', '-md', required=False, type=int, default=250, help='Set the maximum read depth for the pileup. Please increase the maximum value for spike-ins and other highly-covered sequences. (Default 250).')
@click.option('per_chromosome', '--per_chromosome', '-chr', type=str, help='When used, it calculates the modification rates only per the chromosome given.')
@click.option('N_threads', '--N_threads', '-t', default=1, required=True, help='The number of threads to spawn (Default 1).')
@click.option('compress', '--gz', '-z', default=False, is_flag=True, required=False, help='Indicates whether the mods file output will be compressed with gzip (Default False).')
@click.option('directory', '--directory', '-d', required=True, type=str, help='Output directory to save files.')
@click.option('add_underscores', '--add_underscores', '-au', default=False, is_flag=True, required=False, help='Indicates outputting a new reference fasta file with added underscores in the sequence names that is afterwards used for calling. (Default False).')
@click.option('no_information', '--no_information', '-ni', default='*', type=click.Choice(['.', '0', '*', 'NA']), required=False, help='What symbol should be used for a value where no enough quantative information is used. (Default *).')
@click.option('start_clip', '--start_clip', '-scl', required=False, type=int, default=0, help='Set the length of the bases in the start of the reads that will not be used for calling. (Default 8).')
@click.option('end_clip', '--end_clip', '-ecl', required=False, type=int, default=0, help='Set the length of the bases in the end of the reads that will not be used for calling. (Default 0).')
def call(input_file, known_snp, model, reference, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, library,  method, minimum_mapping_quality, adjust_acapq_threshold, add_indels, redo_baq, compute_baq, ignore_orphans, max_depth,per_chromosome, N_threads, directory, compress, single_end, add_underscores, no_information, start_clip, end_clip):
    """Call modified cytosines from a bam or cram file. The output consists of two files, one containing modification counts per nucleotide, the other providing genome-wide statistics per context."""
    cytosine_modification_finder(input_file, known_snp, model, reference, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, library,  method, minimum_mapping_quality, adjust_acapq_threshold, add_indels, redo_baq, compute_baq, ignore_orphans, max_depth, per_chromosome, N_threads, directory, compress, single_end, add_underscores, no_information, start_clip, end_clip)


warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)


logs = logging.getLogger(__name__)

time_b = datetime.now()


def modification_calls_writer(data_mods, compress, data_line, header=False):
    """Outputs the modification calls per position in a tab-delimited format."""
    try:
        if compress == False:
            if header:
                data_line.writerow(["#CHROM", "START", "END", "MOD_LEVEL", "MOD", "UNMOD", "REF", "ALT", "SPECIFIC_CONTEXT", "CONTEXT", "SNV", "TOTAL_DEPTH"])
            data_line.writerow(data_mods)
        else:
            if header:
                data_line.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format("#CHROM", "START", "END", "MOD_LEVEL", "MOD", "UNMOD", "REF", "ALT", "SPECIFIC_CONTEXT", "CONTEXT", "SNV", "TOTAL_DEPTH"))
            data_line.write(
                    '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(data_mods[0], data_mods[1], data_mods[2], data_mods[3],
                                                                 data_mods[4], data_mods[5], data_mods[6], data_mods[7],
                                                                 data_mods[8], data_mods[9], data_mods[10], data_mods[11]))
    except IOError:
        logs.error('asTair cannot write to modification calls file.', exc_info=True)


def statistics_calculator(mean_mod, mean_unmod, data_mod, user_defined_context, context_sample_counts):
    """Calculates the summary statistics of the cytosine modificaton levels."""
    context_sample_counts[data_mod[9]] += 1
    context_sample_counts[data_mod[8]] += 1
    if data_mod[10] == 'No':
        for context in numpy.array(['CpG', 'CHH', 'CHG']):
            if re.match(context, data_mod[9]):
                mean_mod[context] += int(data_mod[4])
                mean_unmod[context] += int(data_mod[5])
        for context in numpy.array(['CAG', 'CCG', 'CTG', 'CTT', 'CCT', 'CAT', 'CTA', 'CTC', 'CAC', 'CAA', 'CCA', 'CCC', 'CGA', 'CGT', 'CGC', 'CGG']):
            if re.match(context, data_mod[8]):
                mean_mod[context] += int(data_mod[4])
                mean_unmod[context] += int(data_mod[5])
        if re.match(r"CN", data_mod[9]):
            mean_mod['CNN'] += int(data_mod[4])
            mean_unmod['CNN'] += int(data_mod[5])
        if user_defined_context and re.match('user defined context', data_mod[9]):
            mean_mod['user defined context'] += int(data_mod[4])
            mean_unmod['user defined context'] += int(data_mod[5])


def context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context, total_contexts, sub_contexts, header, no_information):
    """Writes the summary statistics of the cytosine modificaton levels."""
    with open(file_name, 'a') as statistics_output:
            write_file = csv.writer(statistics_output, delimiter='\t', lineterminator='\n')
            if header == True:
                write_file.writerow(["#CONTEXT", "SPECIFIC_CONTEXT", "MEAN_MODIFICATION_RATE_PERCENT", "TOTAL_POSITIONS", "COVERED_POSITIONS", 'MODIFIED', 'UNMODIFIED'])
                if user_defined_context:
                    write_file.writerow([user_defined_context, "*", safe_rounder(non_zero_division(mean_mod['user defined context'], mean_mod['user defined context'] + mean_unmod['user defined context'], no_information), 3, True), context_total_counts['user defined context'], context_sample_counts['user defined context'], mean_mod['user defined context'], mean_unmod['user defined context']])
            write_file.writerow([context, "*", safe_rounder(non_zero_division(mean_mod[context], mean_mod[context] + mean_unmod[context], no_information), 3, True),
                        context_total_counts[total_contexts]+context_total_counts[total_contexts + 'b'], context_sample_counts[context], mean_mod[context], mean_unmod[context]])
            if len(sub_contexts) >= 1:
                for subcontext in sub_contexts:
                    write_file.writerow(["*", subcontext, safe_rounder(non_zero_division(mean_mod[subcontext], mean_mod[subcontext] + mean_unmod[subcontext], no_information), 3, True), context_total_counts[subcontext], context_sample_counts[subcontext], mean_mod[subcontext], mean_unmod[subcontext]])
            

def final_statistics_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context, no_information):
    """Writes the summary statistics of the cytosine modificaton levels.
    Cytosine modification rate given as the percentage total modified cytosines
    divided by the total number of cytosines covered."""
    try:
        if context == 'all':
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, 'CpG', 'CG', numpy.array(['CGA','CGC', 'CGG', 'CGT']), True, no_information)
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, 'CHG', 'CHG',  numpy.array(['CAG','CCG', 'CTG']), False, no_information)
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, 'CHH', 'CHH',  numpy.array(['CTT', 'CAT', 'CCT', 'CTA', 'CAA', 'CCA', 'CTC', 'CAC', 'CCC']), False, no_information)
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, 'CNN', 'CN',  numpy.array([]), False, no_information)
        elif context == 'CpG':
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context, 'CG',  numpy.array(['CGA','CGC', 'CGG', 'CGT']), True, no_information)
        elif context == 'CHG':
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context, 'CHG',  numpy.array(['CAG','CCG', 'CTG']), True, no_information)
        elif context == 'CHH':
            context_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context, 'CHH',  numpy.array(['CTT', 'CAT', 'CCT', 'CTA', 'CAA', 'CCA', 'CTC', 'CAC', 'CCC']), True, no_information)
    except IOError:
        logs.error('asTair cannot write to statistics summary file.', exc_info=True)


def tuple_handler(read_counts, tuple_with_posinfo, method, no_information):
    if len(tuple_with_posinfo) == 4:
        unmodified_bases = read_counts[tuple_with_posinfo[0]]
        TAAF = read_counts[tuple_with_posinfo[1]]
        RAAF = read_counts[tuple_with_posinfo[2]]
        AAAF = read_counts[tuple_with_posinfo[3]]
    elif len(tuple_with_posinfo) == 8:
        unmodified_bases = read_counts[tuple_with_posinfo[0]] + read_counts[tuple_with_posinfo[4]]
        TAAF = read_counts[tuple_with_posinfo[1]] + read_counts[tuple_with_posinfo[5]]
        RAAF = read_counts[tuple_with_posinfo[2]] + read_counts[tuple_with_posinfo[6]]
        AAAF = read_counts[tuple_with_posinfo[3]] + read_counts[tuple_with_posinfo[7]]
    else:
        unmodified_bases = read_counts[tuple_with_posinfo[0]] + read_counts[tuple_with_posinfo[4]] + read_counts[tuple_with_posinfo[8]] + read_counts[tuple_with_posinfo[12]]
        TAAF = read_counts[tuple_with_posinfo[1]] + read_counts[tuple_with_posinfo[5]] + read_counts[tuple_with_posinfo[9]]
        RAAF = read_counts[tuple_with_posinfo[2]] + read_counts[tuple_with_posinfo[6]] + read_counts[tuple_with_posinfo[10]]
        AAAF = read_counts[tuple_with_posinfo[3]] + read_counts[tuple_with_posinfo[7]] + read_counts[tuple_with_posinfo[11]]
    allele_frequencies = list()
    for AAF in [TAAF, RAAF, AAAF]:
        allele_frequencies.append(non_zero_division(AAF, (AAF + unmodified_bases), no_information))
    if method == 'mCtoT':
        if tuple_with_posinfo[0][1] == 'C':
            ratio_modified = non_zero_division(TAAF,(unmodified_bases+TAAF), no_information)
        elif tuple_with_posinfo[0][1] == 'G':
            ratio_modified = non_zero_division(AAAF, (unmodified_bases+AAAF), no_information)
    elif method == 'CtoT':
        if tuple_with_posinfo[0][1] == 'C':
            ratio_modified = non_zero_division(unmodified_bases,(unmodified_bases+TAAF), no_information)
        elif tuple_with_posinfo[0][1] == 'G':
            ratio_modified = non_zero_division(unmodified_bases, (unmodified_bases+AAAF), no_information)
    return unmodified_bases, TAAF, RAAF, AAAF, allele_frequencies, ratio_modified


def alternative_alele(read_counts, ref):
    if len(read_counts) > 0:
        alt = [key for key in read_counts.keys()][[i for i in read_counts.values()].index(max(read_counts.values()))][1]
    if len(read_counts)==0 or alt not in ['G', 'A', 'C', 'A', 'N'] or alt==ref:
        alt = None
    return alt

def universal_variant_calculation_heuristic(read_counts, unexpected_tuples, ref, method, no_information):
    if len(read_counts)>0 and unexpected_tuples[0][1] != [key for key in read_counts.keys()][[i for i in read_counts.values()].index(max(read_counts.values()))][1]:
        unmodified_bases, TAAF, RAAF, AAAF, allele_frequencies, ratio_modified = tuple_handler(read_counts, unexpected_tuples, method, no_information)
        # add a condition for being a numeric value
        if (non_zero_division(TAAF, (unmodified_bases + TAAF), no_information) >=0.8 and ref=='C') or (non_zero_division(AAAF, (unmodified_bases + AAAF), no_information) >=0.8 and ref=='G'): 
            snp = 'homozygous'
        else:
            snp = None
        alt = alternative_alele(read_counts, ref)
    else:
        snp, alt = None, None
    return snp, alt


def universal_modification_calculation(read_counts, expected_tuples, snp, method, no_information):
    unmodified_bases, TAAF, RAAF, AAAF, allele_frequencies, ratio_modified = tuple_handler(read_counts, expected_tuples, method, no_information)
    if method == 'mCtoT':
        if expected_tuples[0][1] == 'C':
            modified_bases = TAAF
        else:
            modified_bases = AAAF
    else:
        if expected_tuples[0][1] == 'C':
            modified_bases = unmodified_bases
            unmodified_bases = TAAF
        else:
            modified_bases = unmodified_bases
            unmodified_bases = AAAF
    modification_level = safe_rounder(ratio_modified, 3, False)
    return modification_level, modified_bases, unmodified_bases


def context_dictionary(user_defined_context):
    empty_dict = dict()
    for context in ['CHH', 'CHG', 'CpG', 'CNN', 'CAG', 'CCG', 'CTG', 'CTT', 'CCT', 'CAT', 'CTA', 'CTC', 'CAC', 'CAA', 'CCA', 'CCC', 'CGA', 'CGT', 'CGC', 'CGG']:
        empty_dict[context] = 0
    if user_defined_context:
        empty_dict['user defined context'] = 0
    return empty_dict


def clipped_reads(pileups, positions, sequences, start_clip, end_clip):   
    """Simplified read clipping."""
    if positions > start_clip and positions < (pileups.alignment.rlen - end_clip):
        return pileups, sequences   
    else:
        return None, None
    
    
def flags_expectation(modification_information_per_position, position, modification, reference, ignore_orphans, single_end, library):
    """Gives the expected flag-base couples, the reference and the modified base."""
    if reference == 'C':
        expected_references = ['C', 'T', 'G', 'A']
        all_references = ['C', 'T', 'G', 'A']
        if single_end == True:
            expected_flags = [0]
            unexpected_flags = [16]
        else:
            expected_flags = [99, 147]
            unexpected_flags = [83, 163]
            if ignore_orphans == False:
                expected_flags.extend(97, 145)
                unexpected_flags.extend(81, 161)
    elif reference == 'G':
        expected_references = ['G', 'T', 'C', 'A']
        all_references = ['G', 'T', 'C', 'A']
        if single_end == True:
            expected_flags = [16]
            unexpected_flags = [0]
        else:
            expected_flags = [83, 163]
            unexpected_flags = [99, 147]
            if ignore_orphans == False:
                expected_flags.extend(81, 161)
                unexpected_flags.extend(97, 145)
    if library == 'directional':
        desired_tuples = [(flag, ref) for flag in expected_flags for ref in expected_references]
        undesired_tuples = [(flag, ref) for flag in unexpected_flags for ref in all_references]
    elif library == 'reverse':
        desired_tuples = [(flag, ref) for flag in unexpected_flags for ref in expected_references]
        undesired_tuples = [(flag, ref) for flag in expected_flags for ref in all_references]
    return desired_tuples, undesired_tuples


def pileup_summary(modification_information_per_position, position, read_counts, mean_mod, mean_unmod, user_defined_context,
                   header, desired_tuples, undesired_tuples, modification, reference, depth, method, context_sample_counts, ignore_orphans, single_end, compress, data_line, real_snp, model, labels, additional_information, no_information):
    """Creates the modification output per position in the format:
    [chrom, start, end, mod_level, mod, unmod, ref, alt, specific_context, context, snv, total_depth]
    given the strand information and whether the library is pair-end or single-end. The key structure is read_counts
    that contains as dictionary items (read flag, base) tuples.
    Assigns heuristic snv categories of heterozygous and not a snv by using the base ratios of the opposite strand."""
    snp_prob, alt = '*', None
    if real_snp == False:
        snp, alt = universal_variant_calculation_heuristic(read_counts, undesired_tuples, reference, method, no_information)
    else:
        snp = 'WGS_known'
    modification_level, modified_bases, unmodified_bases = universal_modification_calculation(read_counts, desired_tuples, snp, method, no_information)
    if depth < 3:
        label = 'LowCov'
    else:
        label = 'PASS'
    if alt is not None:
        modification = alt
    if snp is None:
        snp = 'No'
    all_data = numpy.array([position[0], position[1], position[1] + 1, modification_level, modified_bases, unmodified_bases, reference, modification, modification_information_per_position[position][0], modification_information_per_position[position][1], snp, depth])
    statistics_calculator(mean_mod, mean_unmod, all_data, user_defined_context, context_sample_counts)
    modification_calls_writer(all_data, compress, data_line, header=header)


def tags_search(read_tags, tag_name, list_):
    list_.extend([i[1] for i in read_tags if i[0] == tag_name])


def clean_pileup(pileups, cycles, modification_information_per_position, mean_mod, mean_unmod, user_defined_context,
                 file_name, method, add_indels, context_sample_counts, ignore_orphans, single_end, compress, data_line, library,
                 true_variants, possible_mods, matched, model, labels, fastas, model_name, no_information, start_clip, end_clip):
    """Takes reads from the piled-up region and calculates modification levels."""
    vclipped_reads = numpy.vectorize(clipped_reads)
    for reads in pileups:
        if cycles == 0:
            header = True
        else:
            header = False
        try:
            if ((reads.reference_name, reads.pos, reads.pos + 1) in modification_information_per_position) or (possible_mods is not None and (reads.reference_name, reads.pos, reads.pos + 1) in possible_mods):
                cycles += 1
                position = (reads.reference_name, reads.pos, reads.pos + 1)
                real_snp = False
                read_counts = defaultdict(int)
                if matched:
                    if true_variants.intersection({position}):
                        real_snp = True
                if modification_information_per_position[position][2] == 'C':
                    modification = 'T'
                    reference = 'C'
                elif modification_information_per_position[position][2] == 'G':
                    modification = 'A'
                    reference = 'G'
                desired_tuples, undesired_tuples = flags_expectation(modification_information_per_position, position,
                                                                            modification, reference, ignore_orphans,
                                                                            single_end, library)
                try:
                    sequences = reads.get_query_sequences(mark_matches=False, mark_ends=False, add_indels=add_indels)
                except AssertionError:
                    logs.exception("Failed getting query sequences (AssertionError, pysam). Please decrease the max_depth parameter.")   
                if start_clip==0 and end_clip==0:
                    clipped_pileups, sequences = reads.pileups, sequences
                else:
                    clipped_pileups, sequences = vclipped_reads(reads.pileups, reads.get_query_positions(), sequences, start_clip, end_clip)
                    clipped_pileups, sequences = clipped_pileups[clipped_pileups != numpy.array(None)], sequences[sequences != numpy.array(None)]
                for pileup, seq in zip_longest(clipped_pileups, sequences, fillvalue='BLANK'):
                    read_counts[(pileup.alignment.flag, seq.upper())] += 1
                if possible_mods is not None and (reads.reference_name, reads.pos, reads.pos + 1) in possible_mods and ((sequences.count(modification) + sequences.count(modification.lower())) != max([sequences.count('A')+sequences.count('a'), sequences.count('C')+sequences.count('c'),
                                                                                                   sequences.count('G')+sequences.count('g'), sequences.count('T')+sequences.count('t')]) or (sequences.count(reference) + sequences.count(reference.lower())) != max([sequences.count('A')+sequences.count('a'), sequences.count('C')+sequences.count('c'),
                                                                                                   sequences.count('G')+sequences.count('g'), sequences.count('T')+sequences.count('t')])):
                    pass
                pileup_summary(modification_information_per_position, position, read_counts, mean_mod, mean_unmod,
                                    user_defined_context, header, desired_tuples, undesired_tuples, modification, reference,
                                    reads.get_num_aligned(), method, context_sample_counts, ignore_orphans, single_end,
                                    compress, data_line, real_snp, model, labels, None, no_information)
                modification_information_per_position.pop(position)
        except Exception:
            continue
    pileups = None


def cytosine_modification_finder(input_file, known_snp, model, reference, context, zero_coverage, skip_clip_overlap, minimum_base_quality, user_defined_context, library, method, minimum_mapping_quality, adjust_acapq_threshold, add_indels, redo_baq, compute_baq, ignore_orphans, max_depth, per_chromosome, N_threads, directory, compress, single_end, add_underscores, no_information, start_clip, end_clip):
    """Searches for cytosine modification positions in the desired contexts and calculates the modificaton levels."""
    time_s = datetime.now()
    logs.info("asTair modification finder started running. {} seconds".format((time_s - time_b).total_seconds()))
    name = path.splitext(path.basename(input_file))[0]
    directory = path.abspath(directory)
    if model == 'none':
        model, labels, model_name = None, None, None
    if os.path.exists(directory) == False:
        raise Exception("The output directory does not exist.")
        sys.exit(1)
    if per_chromosome is None:
        file_name = path.join(directory, name + "_" + method + "_" + context + ".mods")
    else:
        file_name = path.join(directory, name + "_" + method + "_" + per_chromosome + "_" + context + ".mods")
    if no_information == '0':
        no_information = int(no_information)
    if not os.path.isfile(file_name) and not os.path.isfile(file_name + '.gz'):
        try:
            inbam = bam_file_opener(input_file, None, N_threads)
        except Exception:
            sys.exit(1)
        try:
            if per_chromosome is None:
                keys = fasta_splitting_by_sequence(reference, 'keys_only', None, add_underscores, None)
            else:
                keys = [per_chromosome]
        except Exception:
            sys.exit(1)
        contexts, all_keys = sequence_context_set_creation(context, user_defined_context)
        cycles = 0
        mean_mod, mean_unmod = context_dictionary(user_defined_context), context_dictionary(user_defined_context)
        context_total_counts, context_sample_counts = defaultdict(int), defaultdict(int)
        if compress == False:
            calls_output = open(file_name, 'a')
            data_line = csv.writer(calls_output, delimiter='\t', lineterminator='\n')
        else:
            logs.info("Compressing output modification calls file.")
            if sys.version[0] == '3':
                data_line = gzip.open(file_name + '.gz', 'wt', compresslevel=9, encoding='utf8', newline='\n')
            else:
                data_line = gzip.open(file_name + '.gz', 'wt', compresslevel=9)
        true_variants, possible_mods, matched = None, None, False
        for i in range(0, len(keys)):
            time_m = datetime.now()
            logs.info("Starting modification calling on {} chromosome (sequence). {} seconds".format(keys[i], (time_m - time_b).total_seconds()))
            if i == len(keys)-1:
                numbered = 'last'
            else:
                numbered = None
            try:
                fastas = fasta_splitting_by_sequence(reference, keys[i], numbered, add_underscores, None)
                fastas is not None
            except Exception:
                logs.error('The fasta file does not exist or is not indexed.', exc_info=True)
                sys.exit(1)
            modification_information_per_position = context_sequence_search(contexts, all_keys, fastas, keys[i], user_defined_context, context_total_counts, None, None)
            pileups = inbam.pileup(keys[i], ignore_overlaps=skip_clip_overlap, min_base_quality=minimum_base_quality, stepper='samtools',
                                    max_depth=max_depth, redo_baq=redo_baq, ignore_orphans=ignore_orphans, compute_baq=compute_baq,
                                    min_mapping_quality=minimum_mapping_quality, adjust_acapq_threshold=adjust_acapq_threshold)
            if known_snp is not None:
                time_s = datetime.now()
                logs.info("Starting reading SNP information on {} chromosome (sequence). {} seconds".format(keys[i], (time_s - time_m).total_seconds()))
                true_variants, possible_mods = read_vcf(known_snp, keys[i], fastas[keys[i]], N_threads, None, None)
                matched = True
                time_sf = datetime.now()
                logs.info("Reading SNP information on {} chromosome (sequence) has finished. {} seconds".format(keys[i], (time_sf - time_s).total_seconds()))
            clean_pileup(pileups, i, modification_information_per_position, mean_mod, mean_unmod, user_defined_context, file_name, method,
                            add_indels, context_sample_counts, ignore_orphans, single_end, compress, data_line, library, true_variants, possible_mods, matched, model, labels, fastas, model_name, no_information, start_clip, end_clip)
            if zero_coverage:
                for position in modification_information_per_position.keys():
                    if modification_information_per_position[position][2] == 'C':
                        all_data = numpy.array([position[0], position[1], position[1] + 1, no_information, 0, 0, 'C', 'T',
                        modification_information_per_position[position][0], modification_information_per_position[position][1], '*', 0, '*', '*'])
                        modification_calls_writer(all_data, compress, data_line, header=False)
                    elif modification_information_per_position[position][2] == 'G':
                        all_data = numpy.array([position[0], position[1], position[1] + 1, no_information, 0, 0, 'G', 'A',
                        modification_information_per_position[position][0], modification_information_per_position[position][1], '*', 0, '*', '*'])
                        modification_calls_writer(all_data, compress, data_line, header=False)
            modification_information_per_position, true_variants, possible_mods = None, None, None
        inbam.close()
        if per_chromosome is None:
            file_name = path.join(directory, name + "_" + method + "_" + context + ".stats")
        else:
            file_name = path.join(directory, name + "_" + method + "_" + per_chromosome + "_" + context + ".stats")
        final_statistics_output(mean_mod, mean_unmod, user_defined_context, file_name, context_sample_counts, context_total_counts, context, no_information)
        if compress == False:
             calls_output.close()
        else:
            data_line.close()
        time_e = datetime.now()
        logs.info("asTair modification finder finished running. {} seconds".format((time_e - time_b).total_seconds()))
    else:
        logs.error('Mods file with this name exists. Please rename before rerunning.')
        sys.exit(1)

if __name__ == '__main__':
    call()


