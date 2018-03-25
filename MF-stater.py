from Bio import SeqIO
import csv
import itertools
import re
import click


@click.command()
@click.argument('input_file')
@click.argument('fasta_file')
@click.argument('genome')
@click.argument('clip')
def stater_exec(input_file, fasta_file, genome, clip):
    stater(input_file, fasta_file, genome, clip)

def stater(input_file, fasta_file, genome, clip):
    '''input_file is bedgraph file
    fasta_file is the reference genome
    genome is the name of the eukaryotic genome
    clip is number of bases to trim from each genome end
    '''



    cpg_pat = re.compile("CG", re.IGNORECASE)
    chg_pat = re.compile("C(A|C|T)G", re.IGNORECASE)
    chh_pat = re.compile("C(A|C|T)(A|T|C)", re.IGNORECASE)
    gpc_pat = re.compile("GC", re.IGNORECASE)

    chrs = re.compile("chr|CHR", re.IGNORECASE)

    data_CHH = set()
    data_CHG = set()
    data_CpG = set()
    data_GpC = set()


    total_key = set()
    length = []

    record_iter = SeqIO.parse(open(fasta_file), "fasta")
    for chrom in record_iter:
        total_key.add(chrom.id)
        if chrom.id in total_key:
             fasta_seq = chrom.seq.tostring()
             length.append(len(chrom.seq.tostring()))
             for m in chh_pat.finditer(fasta_seq):
                 if m.start() > int(clip) and m.start() + 1 < len(fasta_seq) - int(clip):
                    data_CHH.add((chrom.id, m.start(), m.start() + 1))
             for m in chg_pat.finditer(fasta_seq):
                 if m.start() > int(clip) and m.start() + 1 < len(fasta_seq) - int(clip):
                    data_CHG.add((chrom.id, m.start(), m.start() + 1))
             for m in gpc_pat.finditer(fasta_seq):
                 if m.start() > int(clip) and m.start() + 1 < len(fasta_seq) - int(clip):
                    data_GpC.add((chrom.id, m.start(), m.start() + 1))
                 if m.start() + 1 > int(clip) and m.start() + 2 < len(fasta_seq) - int(clip):
                    data_GpC.add((chrom.id, m.start(), m.start() + 1))
             for m in cpg_pat.finditer(fasta_seq):
                 if m.start() > int(clip) and m.start() + 1 < len(fasta_seq) - int(clip):
                    data_CpG.add((chrom.id, m.start(), m.start() + 1))
                 if m.start() + 1 > int(clip) and m.start() + 2 < len(fasta_seq) - int(clip):
                    data_CpG.add((chrom.id, m.start(), m.start() + 1))




    modification_rate_CHH = []
    modification_rate_CHG = []
    modification_rate_CpG = []
    modification_rate_GpC = []

    key = list(total_key)


    with open(input_file) as handle:
        bedgraph_file = csv.reader(handle, delimiter='\t', lineterminator='\n')
        for row in bedgraph_file:
            data_bedgraph =[(row[0], int(row[1]), int(row[2]), float(row[3]))]
            for k in range(len(key)):
                data_sample_CHH = list(filter(lambda x: (x[0], x[1], x[2]) in data_CHH, data_bedgraph))
                data_sample_CHG = list(filter(lambda x: (x[0], x[1], x[2]) in data_CHG, data_bedgraph))
                data_sample_CpG = list(filter(lambda x: (x[0], x[1], x[2]) in data_CpG, data_bedgraph))
                data_sample_GpC = list(filter(lambda x: (x[0], x[1], x[2]) in data_GpC, data_bedgraph))
                if len(data_sample_CHH) != 0:
                    modification_rate_CHH.append(tuple(itertools.chain([*[x[0] for x in data_sample_CHH], non_zero_div(sum([x[3] for x in data_sample_CHH]),len(data_sample_CHH))])))
                if len(data_sample_CHG) != 0:
                    modification_rate_CHG.append(tuple(itertools.chain([*[x[0] for x in data_sample_CHG], non_zero_div(sum([x[3] for x in data_sample_CHG]),len(data_sample_CHG))])))
                if len(data_sample_CpG) != 0:
                    modification_rate_CpG.append(tuple(itertools.chain([*[x[0] for x in data_sample_CpG], non_zero_div(sum([x[3] for x in data_sample_CpG]), len(data_sample_CpG))])))
                if len(data_sample_GpC) != 0:
                    modification_rate_GpC.append(tuple(itertools.chain([*[x[0] for x in data_sample_GpC], non_zero_div(sum([x[3] for x in data_sample_GpC]),len(data_sample_GpC))])))



    with open(input_file + ".stats", 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        for k in range(len(key)):
            if re.search(chrs, key[k]):
                wr.writerow([genome, "CpG", round(non_zero_div(sum([x[1] for x in modification_rate_CpG]), len([x[1] for x in modification_rate_CpG])),3), len([x[1] for x in modification_rate_CpG])])
                wr.writerow([genome, "CHH", round(non_zero_div(sum([x[1] for x in modification_rate_CHH]), len([x[1] for x in modification_rate_CHH])),3), len([x[1] for x in modification_rate_CHH])])
                wr.writerow([genome, "CHG", round(non_zero_div(sum([x[1] for x in modification_rate_CHG]), len([x[1] for x in modification_rate_CHG])),3), len([x[1] for x in modification_rate_CHG])])
                wr.writerow([genome, "GpC", round(non_zero_div(sum([x[1] for x in modification_rate_GpC]), len([x[1] for x in modification_rate_GpC])),3), len([x[1] for x in modification_rate_GpC])])
            else:
                wr.writerow([key[k], "CpG", round(non_zero_div(sum([x[1] for x in modification_rate_CpG if x[0] == key[k]]), len([x[1] for x in modification_rate_CpG if x[0] == key[k]])),3), len([x[1] for x in modification_rate_CpG if x[0] == key[k]])])
                wr.writerow([key[k], "CHH", round(non_zero_div(sum([x[1] for x in modification_rate_CHH if x[0] == key[k]]), len([x[1] for x in modification_rate_CHH if x[0] == key[k]])),3), len([x[1] for x in modification_rate_CHH if x[0] == key[k]])])
                wr.writerow([key[k], "CHG", round(non_zero_div(sum([x[1] for x in modification_rate_CHG if x[0] == key[k]]), len([x[1] for x in modification_rate_CHG if x[0] == key[k]])),3), len([x[1] for x in modification_rate_CHG if x[0] == key[k]])])
                wr.writerow([key[k], "GpC", round(non_zero_div(sum([x[1] for x in modification_rate_GpC if x[0] == key[k]]), len([x[1] for x in modification_rate_GpC if x[0] == key[k]])),3), len([x[1] for x in modification_rate_GpC if x[0] == key[k]])])


def non_zero_div(x, y):
    if y == 0:
        return 0
    else:
        return x / y


if __name__ == '__main__':
    stater_exec()


