import click
import allel
import itertools
import csv
import re
from array import array
from Bio import SeqIO
from datetime import datetime
import pdb
import ahocorasick

@click.command()
@click.argument('input_file')
@click.argument('fasta_file')
def caller_exec(input_file):
    caller(input_file)

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


    it_OT = splitter_mods(input_file+'OT.vcf.gz', 'OT')
    it_OB = splitter_mods(input_file+'OB.vcf.gz', 'OB')


    data_CHH = set()
    data_CHG = set()
    data_CpG = set()

    CHH_top = ['CAA', 'CAC', 'CAT', 'CCC', 'CCT', 'CCA', 'CTA', 'CTC', 'CTT', 'caa', 'cac', 'cat', 'ccc', 'cct', 'cca', 'cta', 'ctc', 'ctt','CAa', 'CAc', 'CAt', 'CCc', 'CCt', 'CCa', 'CTa', 'CTc', 'CTt', 'caA', 'caC', 'caT', 'ccC', 'ccT', 'ccA', 'ctA', 'ctC', 'ctT', 'cAA', 'cAC', 'cAT', 'cCC', 'cCT', 'cCA', 'cTA', 'cTC', 'cTT','Caa', 'Cac', 'Cat', 'Ccc', 'Cct', 'Cca', 'Cta', 'Ctc', 'Ctt', 'CaA', 'CaC', 'CaT', 'CcC', 'CcT', 'CcA', 'CtA', 'CtC', 'CtT','cAa', 'cAc', 'cAt', 'cCc', 'cCt', 'cCa', 'cTa', 'cTc', 'cTt']
    CHH_bottom = ['GAA', 'GAG', 'GAT', 'GGG', 'GGT', 'GGA', 'GTA', 'GTG', 'GTT', 'gaa', 'gag', 'gat', 'ggg', 'ggt', 'gga', 'gta', 'gtg', 'gtt','GAa', 'GAg', 'GAt', 'GGg', 'GGt', 'GGa', 'GTa', 'GTg', 'GTt', 'gaA', 'gaG', 'gaT', 'ggG', 'ggT', 'ggA', 'gtA', 'gtG', 'gtT', 'gAA', 'gAG', 'gAT', 'gGG', 'gGT', 'gGA', 'gTA', 'gTG', 'gTT','Gaa', 'Gag', 'Gat', 'Ggg', 'Ggt', 'Gga', 'Gta', 'Gtg', 'Gtt', 'GaA', 'GaG', 'GaT', 'GgG', 'GgT', 'GgA', 'GtA', 'GtG', 'GtT', 'gAa', 'gAg', 'gAt', 'gGg', 'gGt', 'gGa', 'gTa', 'gTg', 'gTt']
    CHG_top = ['CAG', 'CCG', 'CTG', 'cag', 'ccg', 'ctg','CAg', 'CCg', 'CTg', 'caG', 'ccG', 'ctG', 'CaG', 'CcG', 'CtG', 'cAg', 'cCg', 'cTg', 'cAG', 'cCG', 'cTG', 'Cag', 'Ccg', 'Ctg']
    CHG_bottom = ['GAC', 'GGC', 'GTC', 'gac', 'ggc', 'gtc','GAc', 'GGc', 'GTc', 'gaC', 'ggC', 'gtC', 'GaC', 'GgC', 'GtC', 'gAc', 'gGc', 'gTc','gAC', 'gGC', 'gTC', 'Gac', 'Ggc', 'Gtc']
    CG = ['CG', 'cg', 'Cg', 'cG']

    auto_CG = ahocorasick.Automaton()
    auto_CHG_top = ahocorasick.Automaton()
    auto_CHH_top = ahocorasick.Automaton()
    auto_CHG_bottom = ahocorasick.Automaton()
    auto_CHH_bottom = ahocorasick.Automaton()

    for pattern in CG:
        auto_CG.add_word(pattern, pattern)

    for pattern in CHG_top:
        auto_CHG_top.add_word(pattern, pattern)

    for pattern in CHG_bottom:
        auto_CHG_bottom.add_word(pattern, pattern)

    for pattern in CHH_top:
        auto_CHH_top.add_word(pattern, pattern)

    for pattern in CHH_bottom:
        auto_CHH_bottom.add_word(pattern, pattern)

    auto_CG.make_automaton()
    auto_CHG_top.make_automaton()
    auto_CHH_top.make_automaton()
    auto_CHG_bottom.make_automaton()
    auto_CHH_bottom.make_automaton()


    fasta_pos_counter = {'CHH': 0, 'CHG': 0, 'CpG': 0}


    record_iter = SeqIO.parse(open(fasta_file), "fasta")
    for chrom in record_iter:
        if chrom.id == chrom_key:
            fasta_seq = chrom.seq.tostring()
            for end_ind, found in auto_CG.iter(fasta_seq):
                fasta_pos_counter['CpG'] += 1
                data_CpG.add((end_ind - 1, end_ind))
                data_CpG.add((end_ind, end_ind + 1))
            for end_ind, found in auto_CHG_top.iter(fasta_seq):
                fasta_pos_counter['CHG'] += 1
                data_CHG.add((end_ind - 2, end_ind - 1))
            for end_ind, found in auto_CHG_bottom.iter(fasta_seq):
                fasta_pos_counter['CHG'] += 1
                data_CHG.add((end_ind - 2, end_ind - 1))
            for end_ind, found in auto_CHH_top.iter(fasta_seq):
                fasta_pos_counter['CHH'] += 1
                data_CHH.add((end_ind - 2, end_ind - 1))
            for end_ind, found in auto_CHH_bottom.iter(fasta_seq):
                fasta_pos_counter['CHH'] += 1
                data_CHH.add((end_ind - 2, end_ind - 1))



    data_CHH.difference_update(data_CpG)
    data_CHG.difference_update(data_CpG)

    print(fasta_pos_counter)

    time_f = datetime.now()
    print((time_f - time_b).total_seconds())


    first_iteration = True
    for row_OT, row_OB in itertools.zip_longest(it_OT, it_OB):
        if row_OB is not None and row_OT is not None:
            row_OT.extend(row_OB)
            row_OT.sort()
            data_sample_CpG = list(filter(lambda x: (x[1], x[2]) in data_CpG, row_OT))
            data_sample_CHH = list(filter(lambda x: (x[1], x[2]) in data_CHH, row_OT))
            data_sample_CHG = list(filter(lambda x: (x[1], x[2]) in data_CHG, row_OT))
            data_writer(input_file + "_CHH.mods", data_sample_CHH, header=first_iteration)
            data_writer(input_file + "_CHG.mods", data_sample_CHG, header=first_iteration)
            data_writer(input_file + "_CpG.mods", data_sample_CpG, header=first_iteration)
            data_writer(input_file + 'total.mods', row_OT, header=first_iteration)
            bed_graph_row = []
            for element in row_OT:
                bed_graph_row.append(element[0:4])
            data_writer(input_file + 'total.bedGraph', bed_graph_row, header=False)
            first_iteration = False
        elif row_OT is not None:
            data_sample_CpG = list(filter(lambda x: (x[1], x[2]) in data_CpG, row_OT))
            data_sample_CHH = list(filter(lambda x: (x[1], x[2]) in data_CHH, row_OT))
            data_sample_CHG = list(filter(lambda x: (x[1], x[2]) in data_CHG, row_OT))
            data_writer(input_file + "_CHH.mods", data_sample_CHH, header=first_iteration)
            data_writer(input_file + "_CHG.mods", data_sample_CHG, header=first_iteration)
            data_writer(input_file + "_CpG.mods", data_sample_CpG, header=first_iteration)
            data_writer(input_file + 'total.mods', row_OT, header=first_iteration)
            bed_graph_row = []
            for element in row_OT:
                bed_graph_row.append(element[0:4])
            data_writer(input_file + 'total.bedGraph', bed_graph_row, header=False)
            first_iteration = False
        elif row_OB is not None:
            data_sample_CpG = list(filter(lambda x: (x[1], x[2]) in data_CpG, row_OB))
            data_sample_CHH = list(filter(lambda x: (x[1], x[2]) in data_CHH, row_OB))
            data_sample_CHG = list(filter(lambda x: (x[1], x[2]) in data_CHG, row_OB))
            data_writer(input_file + "_CHH.mods", data_sample_CHH, header=first_iteration)
            data_writer(input_file + "_CHG.mods", data_sample_CHG, header=first_iteration)
            data_writer(input_file + "_CpG.mods", data_sample_CpG, header=first_iteration)
            data_writer(input_file + 'total.mods', row_OB, header=first_iteration)
            bed_graph_row = []
            for element in row_OB:
                bed_graph_row.append(element[0:4])
            data_writer(input_file + 'total.bedGraph', bed_graph_row, header=False)
            first_iteration = False


    time_e = datetime.now()
    print((time_e - time_b).total_seconds())

def splitter_mods(input_file, strand):
    file_fields, file_samples, file_headers, file_iter = allel.iter_vcf_chunks(input_file, fields='*', alt_number=4,
                                                                               buffer_size=50000, chunk_length=100000)
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
    indel = row_data['variants/INDEL'][:]*1
    index_indel = [(i, j) for i, j in enumerate(indel)]
    not_indel = [x[0] for x in index_indel if x[1] == 0]
    chrom = row_data['variants/CHROM'][not_indel]
    pos = row_data['variants/POS'][not_indel]
    ref = row_data['variants/REF'][not_indel]
    alt = row_data['variants/ALT'][not_indel]
    mods = row_data[alele][not_indel, 0, 1:]
    org = row_data[alele][not_indel, 0, 0]
    depth = row_data['calldata/DP'][not_indel, 0]
    results = [(i, *j, d) for i, j, d in itertools.zip_longest(ref, alt, depth, fillvalue='BLANK')]
    indexed = [(i, *j) for i, j in enumerate(results)]
    index_all = [x[0] for x in indexed if REF_base in x[1] and x[-1] > 0]
    index_mod = [(x[0], x.index(ALT_base) - 2) for x in indexed if ALT_base in x[:] and REF_base in x[1]]
    set_mod = [x[0] for x in index_mod]
    index_nmod = list(filter(lambda x: x not in set_mod, index_all))
    data_mod = [(chrom[i], pos[i] - 1, pos[i], mods[i, k]/(org[i] + mods[i, k]), mods[i, k], org[i], ALT_base, REF_base, depth[i]) for n in range(int(len(index_mod)/crows) + 1) for i, k in index_mod[n * crows:n * crows + crows]]
    data_nmod = [(chrom[i], pos[i] - 1, pos[i], 0, 0, org[i], ALT_base, REF_base, depth[i]) for n in range(int(len(index_nmod)/crows) + 1) for i in index_nmod[n * crows:n * crows + crows]]
    data_mod.extend(data_nmod)
    data_mod.sort()
    return data_mod


crows = 50000

if __name__ == '__main__':
    caller_exec()


