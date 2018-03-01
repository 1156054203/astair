import click
import pysam
import tempfile
import os
import sys
import scipy
import allel
import itertools
import csv
import pdb


@click.command()
@click.argument('input_file')
def caller_exec(input_file):
    caller(input_file)

def caller(input_file):
    '''input_file is file name without extension
    '''

    it_OT = splitter_mods(input_file+'OT.vcf', 'OT')
    it_OB = splitter_mods(input_file+'OB.vcf', 'OB')


    first_iteration = True
    for row_OT, row_OB in itertools.zip_longest(it_OT, it_OB):
        data_writer(input_file + 'OT.mods', row_OT, header=first_iteration)
        data_writer(input_file + 'OB.mods', row_OB, header=first_iteration)
        row_OT.extend(row_OB)
        row_OT.sort()
        data_writer(input_file + 'total.mods', row_OT, header=first_iteration)
        bed_graph_row = []
        for element in row_OT:
            bed_graph_row.append(element[0:4])
        data_writer(input_file + 'total.bedGraph', bed_graph_row, header=False)
        first_iteration = False


def splitter_mods(input_file, strand):
    file_fields, file_samples, file_headers, file_iter = allel.iter_vcf_chunks(input_file, fields='*', alt_number=4,
                                                                               buffer_size=16384, chunk_length=65536)
    if strand == 'OT':
        ref_base = 'C'
        alt_base = 'T'
    else:
        ref_base = 'G'
        alt_base = 'A'
    for row in file_iter:
        row_data = row[0]
        data = variants_aggregator(row_data, 'calldata/AD', ref_base, alt_base)
        yield data


def data_writer(input_name, data_mods, header = False):
    with open(input_name, 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        if header:
            wr.writerow(["CHROM", "START", "END", "MOD_LEVEL", "MOD", "NMOD", "ALT", "REF", "DEPTH"])
        for row in data_mods:
            wr.writerow(row)


def variants_aggregator(row_data, alele, REF_base, ALT_base):
    chrom = row_data['variants/CHROM'][:]
    pos = row_data['variants/POS'][:]
    ref = row_data['variants/REF'][:]
    alt = row_data['variants/ALT'][:]
    mods = row_data[alele][:, 0, 1:]
    org = row_data[alele][:, 0, 0]
    depth = row_data['calldata/DP'][:, 0]
    results = [(i, *j) for i, j in itertools.zip_longest(ref, alt, fillvalue='BLANK')]
    indexed = [(i, *j) for i, j in enumerate(results)]
    index = [(x[0], x.index(ALT_base) - 2) for x in indexed if
             ALT_base in x[:] and REF_base in x[1]]  # gives row of c and t and t position + 1
    data = [(chrom[i], pos[i] - 1, pos[i], mods[i, k] / (org[i] + mods[i, k]), mods[i, k], org[i], ALT_base, REF_base,
             depth[i]) for n in range(int(len(index) / 1000) + 1) for i, k in index[n * 1000:n * 1000 + 1000] if
            mods[i, k] > 0]
    return data

if __name__ == '__main__':
    caller_exec()
