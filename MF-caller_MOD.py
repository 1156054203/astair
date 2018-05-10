import click
import allel
import itertools
import csv
import re
import numpy
from array import array
from datetime import datetime
from Bio import SeqIO
import ahocorasick
import warnings


@click.command()
@click.option('input_file', '--input_file', required=True, help='VCF format file containing aggregated variants. Only one file name (OT or OB) is required to run for both strands')
@click.option('fasta_file', '--fasta_file', required=True, help='DNA sequence in fasta format used for aligning the sequencing reads and mpileup.')
@click.option('zero_coverage', '--zero_coverage', default = False, is_flag = True, help='When set to True, outputs positions not covered in the vcf file. Uncovering zero coverage positions takes longer time than using the default option.')
def caller_exec(input_file, fasta_file, zero_coverage):
    caller(input_file, fasta_file, zero_coverage)

def caller(input_file, fasta_file, zero_coverage):

    name = re.search(r'^.*(?=(?:OT|OB).vcf.gz)', input_file).group()

    record_iter = SeqIO.parse(open(fasta_file), "fasta")

    time_b = datetime.now()

    if zero_coverage == True:

        letters_top = ['A', 'C', 'T', 'a', 'c', 't']
        letters_bottom = ['A', 'G', 'T', 'a', 'G', 't']

        CHG_top = [y + x + z for x in letters_top for y in ['C', 'c'] for z in ['G', 'g']]
        CHH_top = [y + x + x for x in letters_top for y in ['C', 'c']]
        CHG_bottom = [y + x + z for x in letters_bottom for y in ['G', 'g'] for z in ['C', 'c']]
        CHH_bottom = [y + x + x for x in letters_bottom for y in ['G', 'g']]
        CG = [y + z + x for x in letters_top for y in ['C', 'c'] for z in ['G', 'g']]

        auto_CG, auto_CHG_top, auto_CHH_top, auto_CHG_bottom, auto_CHH_bottom = ahocorasick.Automaton(), ahocorasick.Automaton(), ahocorasick.Automaton(), ahocorasick.Automaton(), ahocorasick.Automaton()

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

        for chrom in record_iter:
            sequences_key.update({chrom.id: len(chrom.seq.tostring())})
            fasta_seq = chrom.seq.tostring()
            for end_ind, found in auto_CG.iter(fasta_seq):
                data_context.add((chrom.id, end_ind - 2, end_ind - 1, found.upper()))
                data_context.add((chrom.id, end_ind - 1, end_ind, found.upper()))
            for end_ind, found in auto_CHG_top.iter(fasta_seq):
                data_context.add((chrom.id, end_ind - 2, end_ind - 1, found.upper()))
            for end_ind, found in auto_CHG_bottom.iter(fasta_seq):
                data_context.add((chrom.id, end_ind - 2, end_ind - 1, found.upper()))
            for end_ind, found in auto_CHH_top.iter(fasta_seq):
                data_context.add((chrom.id, end_ind - 2, end_ind - 1, found.upper()))
            for end_ind, found in auto_CHH_bottom.iter(fasta_seq):
                data_context.add((chrom.id, end_ind - 2, end_ind - 1, found.upper()))

    time_d = datetime.now()

    print("Cytosine positions are loaded. {}".format((time_d - time_b).total_seconds()))

    def splitter_mods(input_file, strand):
        file_fields, file_samples, file_headers, file_iter = allel.iter_vcf_chunks(input_file, fields='*', alt_number=4,
                                                                                   buffer_size=5000, chunk_length=10000)
        if strand == 'OT':
            ref_base = 'C'
            alt_base = 'T'
        else:
            ref_base = 'G'
            alt_base = 'A'

        previous_data = None

        for row in file_iter:
            current_data = row[0]
            if previous_data is not None:
                current_data = concat_previous_data(current_data, previous_data)
                previous_data = None
            if zero_coverage == True:
                data, coef, data_zero = variants_aggregator(current_data, 'calldata/AD', ref_base, alt_base)
            else:
                data, coef = variants_aggregator(current_data, 'calldata/AD', ref_base, alt_base)
            if coef > 0:
                previous_data = filter_incomplete_pattern(current_data, coef)
            if zero_coverage == True:
                yield data, data_zero
            else:
                    yield data

    def concat_previous_data(current_data, previous_data):
        for key, value in previous_data.items():
            previous_data[key] = numpy.concatenate((value, current_data[key]), axis = 0)
        return previous_data

    def filter_incomplete_pattern(data, coef):
        result = {'variants/INDEL' : None, 'variants/CHROM' : None, 'variants/POS' : None, 'variants/REF' : None, 'variants/ALT' : None, 'calldata/AD': None,
                       'calldata/DP': None}
        for column in ['variants/INDEL', 'variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT', 'calldata/AD',
                       'calldata/DP']:
            index = numpy.searchsorted(data['variants/POS'], coef)

            result[column] = data[column][index:]
        return result

    def data_writer(input_name, data_mods, header = False):
        with open(input_name, 'a', newline='') as myfile:
            wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
            if header:
                wr.writerow(["CHROM", "START", "END", "MOD_LEVEL", "MOD", "NMOD", "ALT", "REF", "DEPTH", "CONTEXT"])
            for row in data_mods:
                wr.writerow(row)

    def variants_aggregator(row_data, alele, REF_base, ALT_base):
        cycle = datetime.now()
        totals_context = {}
        totals_found = {}
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
        zipper = itertools.zip_longest(chrom, pos, fillvalue='BLANK')
        for element in zipper:
            totals_found[element] = 1
        indexed = [(i, *j) for i, j in enumerate(results)]
        index_all = [x[0] for x in indexed if REF_base in x[1]]
        index_mod = [(x[0], x.index(ALT_base) - 2) for x in indexed if ALT_base in x[:] and REF_base in x[1]]
        set_mod = [x[0] for x in index_mod]
        index_nmod = list(filter(lambda x: x not in set_mod, index_all))
        data_mod = [(chrom[i], pos[i] - 1, pos[i], mods[i, k]/(org[i] + mods[i, k]), mods[i, k], org[i], ALT_base, REF_base, depth[i], ref[i]+ref[i+1]+ref[i+2] if i + 2 <= len(ref) - 1 and chrom[i] == chrom[i+2] else ref[i]) for n in range(int(len(index_mod)/crows) + 1) for i, k in index_mod[n * crows:n * crows + crows]]
        data_nmod = [(chrom[i], pos[i] - 1, pos[i], 0, 0, org[i], ALT_base, REF_base, depth[i], ref[i]+ref[i+1]+ref[i+2] if i + 2 <= len(ref) - 1 and chrom[i] == chrom[i+2] else ref[i]) for n in range(int(len(index_nmod)/crows) + 1) for i in index_nmod[n * crows:n * crows + crows]]
        position = max(data_nmod[-1][1], data_mod[-1][1])
        if zero_coverage == True:
            for tuples in data_context:
                context_index = (tuples[0], tuples[1])
                totals_context[context_index] = tuples
            current_data_context = [values for keys, values in totals_context.items() if totals_found.keys() not in keys and max(max(all_coeffs), pos[0]) < values[1] < position]
            current_data_context.sort()
            data_zero = [(x[0], x[1], x[2], 'None', 0, 0, x[3][0].upper(), x[3][0].upper(), 0, x[3]) for x in current_data_context]
        coef = 0
        if len(data_nmod) > 0 and len(data_mod) > 0:
            if data_nmod[-1][0] == data_mod[-1][0]:
                item = [value for key, value in sequences_key.items() if data_mod[-1][0] in key.lower()]
                if item > max(data_mod[-1][1], data_nmod[-1][1]) + context_count - 1:
                    if len(data_mod[-1][9]) > 1 and len(data_nmod[-1][9]) > 1:
                        data_mod.extend(data_nmod)
                        data_mod.sort(key=lambda x: x[0:2])
                        stater(data_mod)
                        index_occ.append(max(data_nmod[-1][1], data_mod[-1][1]))
                    if len(data_mod[-1][9]) > 1 and len(data_nmod[-1][9]) == 1:
                        data_nmod = data_nmod[:-2]
                        coef = data_nmod[-1][1] + 2
                        data_mod = [x for x in data_mod if x[1] <= coef - 2]
                        data_mod.extend(data_nmod)
                        data_mod.sort(key=lambda x: x[0:2])
                        stater(data_mod)
                        index_occ.append(max(data_nmod[-1][1], data_mod[-1][1]))
                    if len(data_mod[-1][9]) == 1 and len(data_nmod[-1][9]) > 1:
                        data_mod = data_mod[:-2]
                        coef = data_mod[-1][1] + 2
                        data_nmod = [x for x in data_nmod if x[1] <= coef - 2]
                        data_mod.extend(data_nmod)
                        data_mod.sort(key=lambda x: x[0:2])
                        stater(data_mod)
                        index_occ.append(max(data_nmod[-1][1], data_mod[-1][1]))
                    if len(data_mod[-1][9]) == 1 and len(data_nmod[-1][9]) == 1:
                        data_mod = data_mod[:-2]
                        data_nmod = data_nmod[:-2]
                        coef = max(data_mod[-1][1], data_nmod[-1][1]) + 2
                        data_mod.extend(data_nmod)
                        data_mod.sort(key=lambda x: x[0:2])
                        stater(data_mod)
                        index_occ.append(max(data_nmod[-1][1], data_mod[-1][1]))
                else:
                    data_mod.extend(data_nmod)
                    data_mod.sort(key=lambda x: x[0:2])
                    stater(data_mod)
        else:
             data_mod.extend(data_nmod)
             data_mod.sort(key=lambda x: x[0:2])
             stater(data_mod)
        if coef != 0:
            all_coeffs.append(coef)
        else:
            all_coeffs.append(data_mod[-1][1])
        if zero_coverage == True:
            return data_mod, coef, data_zero
        else:
            return data_mod, coef

    def stater(data_mod):
        modification_rate_CpG = [x[3] for x in data_mod if re.match(r"CG.", x[9]) or re.match(r"GC.", x[9])]
        modification_rate_CHG = [x[3] for x in data_mod if re.match(r"(?=(C(A|C|T)G))", x[9]) or re.match(r"(?=(G(A|G|T)C))", x[9])]
        modification_rate_CHH = [x[3] for x in data_mod if re.match(r"(?=(C(A|C|T)(A|C|T)))", x[9]) or re.match(r"(?=(G(A|G|T)(A|G|T)))", x[9])]
        modification_rate_unknown = [x[3] for x in data_mod if len(x[9]) == 1]
        mean_mod['CpG'] += round(non_zero_div(sum([x for x in modification_rate_CpG if x is not str]), len([x for x in modification_rate_CpG if x is not str])), 3)
        mean_mod['CHG'] += round(non_zero_div(sum([x for x in modification_rate_CHG if x is not str]), len([x for x in modification_rate_CHG if x is not str])), 3)
        mean_mod['CHH'] += round(non_zero_div(sum([x for x in modification_rate_CHH if x is not str]), len([x for x in modification_rate_CHH if x is not str])), 3)
        mean_mod['Unknown'] += round(non_zero_div(sum([x for x in modification_rate_unknown if x is not str]), len([x for x in modification_rate_unknown if x is not str])), 3)
        mod_len['CpG'] += 1
        mod_len['CHG'] += 1
        mod_len['CHH'] += 1
        mod_len['Unknown'] += 1

    it_OT = splitter_mods(name+'OT.vcf.gz', 'OT')
    it_OB = splitter_mods(name+'OB.vcf.gz', 'OB')

    first_iteration = True
    if zero_coverage == True:
        for row_OT, row_OB in itertools.zip_longest(it_OT, it_OB):
            if row_OB is not None and row_OT is not None:
                row_OT[0].extend(row_OB[0])
                row_OT[1].extend(row_OB[1])
                row_OT[0].sort()
                index_covered = ([(x[0], x[1]) for x in row_OT[0]])
                index_covered.sort()
                if len(index_occ) > 2:
                    unique_zero = set([x for x in row_OT[1] if (x[0], x[1]) not in index_covered and max(index_occ[-2], index_occ[-1]) < x[1] < index_covered[-1][1]])
                if len(index_occ) >= 1 and max(index_occ[-2], index_occ[-1]) < index_covered[0][1] + 1000:
                    unique_zero = set([x for x in row_OT[1] if (x[0], x[1]) not in index_covered and max(index_occ[-2], index_occ[-1]) < x[1] < index_covered[-1][1]])
                else:
                    unique_zero = set([x for x in row_OT[1] if (x[0], x[1]) not in index_covered and index_covered[0][1] < x[1] < index_covered[-1][1]])
                unique_zero = list(unique_zero)
                unique_zero.sort()
                row_OT[0].extend(unique_zero)
                row_OT[0].sort(key=lambda x: x[0:2])
                data_writer(name + 'total.mods', row_OT[0], header=first_iteration)
                first_iteration = False
            elif row_OT is not None:
                if len(index_occ) > 2:
                    unique_zero = set([x for x in row_OT[1] if (x[0], x[1]) not in index_covered and max(index_occ[-2], index_occ[-1]) < x[1] < index_covered[-1][1]])
                if max(index_occ[-2], index_occ[-1]) < index_covered[0][1] + 1000:
                    unique_zero = set([x for x in row_OT[1] if (x[0], x[1]) not in index_covered and max(index_occ[-2], index_occ[-1]) < x[1] < index_covered[-1][1]])
                else:
                    unique_zero = set([x for x in row_OT[1] if (x[0], x[1]) not in index_covered and index_covered[0][1] < x[1] < index_covered[-1][1]])
                unique_zero = list(unique_zero)
                row_OT[0].extend(unique_zero)
                row_OT[0].sort(key=lambda x: x[0:2])
                data_writer(name + 'total.mods', row_OT[0], header=first_iteration)
                first_iteration = False
            elif row_OB is not None:
                if len(index_occ) > 2:
                    unique_zero = set([x for x in row_OB[1] if (x[0], x[1]) not in index_covered and max(index_occ[-2], index_occ[-1]) < x[1] < index_covered[-1][1]])
                if max(index_occ[-2], index_occ[-1]) < index_covered[0][1] + 1000:
                    unique_zero = set([x for x in row_OB[1] if (x[0], x[1]) not in index_covered and max(index_occ[-2], index_occ[-1]) < x[1] < index_covered[-1][1]])
                else:
                    unique_zero = set([x for x in row_OB[1] if (x[0], x[1]) not in index_covered and index_covered[0][1] < x[1] < index_covered[-1][1]])
                unique_zero = list(unique_zero)
                row_OB[0].extend(unique_zero)
                row_OB[0].sort(key=lambda x: x[0:2])
                data_writer(name + 'total.mods', row_OB[0], header=first_iteration)
                first_iteration = False
    else:
        for row_OT, row_OB in itertools.zip_longest(it_OT, it_OB):
            if row_OB is not None and row_OT is not None:
                row_OT.extend(row_OB)
                row_OT.sort(key=lambda x: x[0:2])
                data_writer(name + 'total.mods', row_OT, header=first_iteration)
                first_iteration = False
            elif row_OT is not None:
                data_writer(name + 'total.mods', row_OT, header=first_iteration)
                first_iteration = False
            elif row_OB is not None:
                data_writer(name + 'total.mods', row_OB, header=first_iteration)
                first_iteration = False


    with open(name + ".stats", 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        wr.writerow(["mean CpG modification rate: ", round(non_zero_div(mean_mod['CpG'],mod_len["CpG"]),3)])
        wr.writerow(["mean CHG modification rate: ", round(non_zero_div(mean_mod['CHG'],mod_len["CHG"]),3)])
        wr.writerow(["mean CHH modification rate: ", round(non_zero_div(mean_mod['CHH'],mod_len["CHH"]),3)])
        wr.writerow(["mean CNN modification rate: ", round(non_zero_div(mean_mod['Unknown'],mod_len["Unknown"]),3)])


    time_e = datetime.now()
    print("asTair modification finder finished running. {}".format((time_e - time_b).total_seconds()))

index_occ = []
sequences_key = dict()
data_context = set()
all_coeffs = [0]
crows = 5000
context_count = 3
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

mean_mod = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0}
mod_len = {'CHH': 0, 'CHG': 0, 'CpG': 0, 'Unknown': 0}

def non_zero_div(x, y):
    if y == 0:
        return 0
    else:
        return x / y


if __name__ == '__main__':
    caller_exec()


