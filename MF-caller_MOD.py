import click
import allel
import itertools
import csv
import re
from Bio import SeqIO
from datetime import datetime



@click.command()
@click.argument('input_file')
@click.argument('fasta_file')
def caller_exec(input_file, fasta_file):
    caller(input_file, fasta_file)

def caller(input_file):
    '''input_file is file name without extension
    fasta_file is the reference genome
    '''

    time_b = datetime.now()

    match_name = re.search(r'[^\.]*\.(.*)', input_file)
    if not match_name:
        raise ValueError('File name does not match any chromosome name.')
    chrom_key = match_name.group(1)
    print(chrom_key)

    it_OT = splitter_mods(input_file+'OT.vcf', 'OT')
    it_OB = splitter_mods(input_file+'OB.vcf', 'OB')

    cpg_pat = re.compile("CG", re.IGNORECASE)
    chg_pat = re.compile("C(A|C|T)G", re.IGNORECASE)
    chh_pat = re.compile("C(A|C|T)(A|T|C)", re.IGNORECASE)
    gch_pat = re.compile("GC(A|C|T)", re.IGNORECASE)
    hcg_pat = re.compile("(A|C|T)CG", re.IGNORECASE)

    data_CHH = set()
    data_CHG = set()
    data_CpG = set()
    data_HCG = set()
    data_GCH = set()


    fasta_pos_counter = {'CHH':0, 'CHG':0, 'CpG':0, 'GCH':0, 'HCG':0}
    #pdb.set_trace()
    record_iter = SeqIO.parse(open(fasta_file), "fasta")
    for chrom in record_iter:
        if chrom.id == chrom_key:
            fasta_seq = chrom.seq.tostring()
            for m in chh_pat.finditer(fasta_seq):
                fasta_pos_counter['CHH'] += 1
                data_CHH.add((m.start(), m.start() + 1))
            for m in chg_pat.finditer(fasta_seq):
                fasta_pos_counter['CHG'] += 1
                data_CHG.add((m.start(), m.start() + 1))
            for m in hcg_pat.finditer(fasta_seq):
                fasta_pos_counter['HCG'] += 1
                for pos in [m.start() + 1, m.start() + 2]:
                    data_CpG.add((pos, pos + 1))
            for m in gch_pat.finditer(fasta_seq):
                fasta_pos_counter['GCH'] += 1
                for pos in [m.start(), m.start() + 1]:
                    data_GCH.add((pos, pos + 1))
            for m in cpg_pat.finditer(fasta_seq):
                fasta_pos_counter['CpG'] += 1
                for pos in [m.start(), m.start() + 1]:
                    data_HCG.add((pos, pos + 1))


    print(fasta_pos_counter)

    time_f = datetime.now()
    print((time_f - time_b).total_seconds())

    first_iteration = True
    for row_OT, row_OB in itertools.zip_longest(it_OT, it_OB):
        row_OT.extend(row_OB)
        row_OT.sort()
        data_sample_CHH = list(filter(lambda x: (x[1], x[2]) in data_CHH, row_OT))
        data_sample_CHG = list(filter(lambda x: (x[1], x[2]) in data_CHG, row_OT))
        data_sample_CpG = list(filter(lambda x: (x[1], x[2]) in data_CpG, row_OT))
        data_sample_HCG = list(filter(lambda x: (x[1], x[2]) in data_HCG, row_OT))
        data_sample_GCH = list(filter(lambda x: (x[1], x[2]) in data_GCH, row_OT))
        data_writer(input_file + "_CHH.mods", data_sample_CHH, header=first_iteration)
        data_writer(input_file + "_CHG.mods", data_sample_CHG, header=first_iteration)
        data_writer(input_file + "_CpG.mods", data_sample_CpG, header=first_iteration)
        data_writer(input_file + "_GCH.mods", data_sample_GCH, header=first_iteration)
        data_writer(input_file + "_HCG.mods", data_sample_HCG, header=first_iteration)
        data_writer(input_file + 'total.mods', row_OT, header=first_iteration)
        bed_graph_row = []
        for element in row_OT:
            bed_graph_row.append(element[0:4])
        data_writer(input_file + 'total.bedGraph', bed_graph_row, header=False)
        first_iteration = False

    time_e = datetime.now()
    print((time_e - time_b).total_seconds())

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
    results = [(i, *j, d) for i, j, d in itertools.zip_longest(ref, alt, depth, fillvalue='BLANK')]
    indexed = [(i, *j) for i, j in enumerate(results)]
    index_all = [x[0] for x in indexed if REF_base in x[1] and x[-1] > 0]
    index_mod = [(x[0], x.index(ALT_base) - 2) for x in indexed if ALT_base in x[:] and REF_base in x[1]]
    set_mod = [x[0] for x in index_mod]
    index_nmod = list(filter(lambda x: x not in set_mod, index_all))
    data_mod = [(chrom[i], pos[i] - 1, pos[i], mods[i, k]/(org[i] + mods[i, k]), mods[i, k], org[i], ALT_base, REF_base,
             depth[i]) for n in range(int(len(index_mod)/1000) + 1) for i, k in index_mod[n * 1000:n * 1000 + 1000]]
    data_nmod = [(chrom[i], pos[i] - 1, pos[i], 0, 0, org[i], ALT_base, REF_base,
             depth[i]) for n in range(int(len(index_nmod)/1000) + 1) for i in index_nmod[n * 1000:n * 1000 + 1000]]
    data_mod.extend(data_nmod)
    data_mod.sort()
    return data_mod


if __name__ == '__main__':
    caller_exec()


