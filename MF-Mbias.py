import click
import pysam
import re
import itertools
import csv
import numpy
import matplotlib.pyplot as pyp
pyp.style.use('seaborn-whitegrid')


@click.command()
@click.option('input_file', '--input_file', required=True, help='BAM format file containing sequencing reads.')
@click.option('read_length', '--read_length', type=int, required=True, help='The read length is needed to calculate the M-bias.')

def Mbias_exec(input_file, read_length):
    Mbias(input_file, read_length)

def Mbias(input_file, read_length):

    name = re.search(r'^.*(?=.bam)', input_file).group()

    for i in range(0, read_length):
        read1_all_CHH[i] = 0
        read1_all_CHG[i] = 0
        read1_all_CpG[i] = 0
        read1_mods_CHH[i] = 0
        read1_mods_CHG[i] = 0
        read1_mods_CpG[i] = 0
        read2_all_CHH[i] = 0
        read2_all_CHG[i] = 0
        read2_all_CpG[i] = 0
        read2_mods_CHH[i] = 0
        read2_mods_CHG[i] = 0
        read2_mods_CpG[i] = 0

    inbam = pysam.AlignmentFile(input_file, "rb")
    bam_fetch = inbam.fetch(until_eof=True)
    for read in bam_fetch:
        if read.flag == 99 or read.flag == 83:
            if isinstance(read.tags[0][1], str):
                read_data = read.tags[0][1]
            elif isinstance(read.tags[1][1], str):
                read_data = read.tags[1][1]
            if re.search(regs, read_data, re.IGNORECASE):
                changes = [int(s) for s in re.findall(r'\d+', read_data)]
                non_overlap = [x + 1 if x == 0 else x for x in changes]
                names = list(re.findall(r'[^\W\d_]+', read_data))
                positions = [x + 1 for x in list(itertools.accumulate(non_overlap))]
                positions = positions[:-1]
                sequence = list(read.query_sequence)
                for k in range(len(positions) - 1):
                    element = positions[k]
                    sequence[element] = names[k]
                reads = "".join(sequence)
                cpg_all = [m.start() for m in re.finditer(r'CG', reads, re.IGNORECASE)]
                chg_all = [m.start() for m in re.finditer(r'C(A|C|T)G', reads, re.IGNORECASE)]
                chh_all = [m.start() for m in re.finditer(r'C(A|C|T)(A|T|C)', reads, re.IGNORECASE)]
                cpg_mods = [x for x in positions if x in cpg_all]
                chg_mods = [x for x in positions if x in chg_all]
                chh_mods = [x for x in positions if x in chh_all]
                for i in range(0,len(reads)):
                    if i in chh_mods:
                        read1_mods_CHH[i] += 1
                    if i in chg_mods:
                        read1_mods_CHG[i] += 1
                    if i in cpg_mods:
                        read1_mods_CpG[i] += 1
                for i in range(0,len(reads)):
                    if i in chh_all:
                        read1_all_CHH[i] += 1
                    if i in chg_all:
                        read1_all_CHG[i] += 1
                    if i in cpg_all:
                        read1_all_CpG[i] += 1
            else:
                sequence = list(read.query_sequence)
                reads = "".join(sequence)
                cpg_all = [m.start() for m in re.finditer(r'CG', reads, re.IGNORECASE)]
                chg_all = [m.start() for m in re.finditer(r'C(A|C|T)G', reads, re.IGNORECASE)]
                chh_all = [m.start()  for m in re.finditer(r'C(A|C|T)(A|T|C)', reads, re.IGNORECASE)]
                for i in range(0,len(reads)):
                    if i in chh_all:
                        read1_all_CHH[i] += 1
                    if i in chg_all:
                        read1_all_CHG[i] += 1
                    if i in cpg_all:
                        read1_all_CpG[i] += 1
        if read.flag == 147 or read.flag == 163:
            if isinstance(read.tags[0][1], str):
                read_data = read.tags[0][1]
            elif isinstance(read.tags[1][1], str):
                read_data = read.tags[1][1]
            if re.search(regs, read_data, re.IGNORECASE):
                changes = [int(s) for s in re.findall(r'\d+', read_data)]
                non_overlap = [x + 1 if x == 0 else x for x in changes]
                names = list(re.findall(r'[^\W\d_]+', read_data))
                positions = [x + 1 for x in list(itertools.accumulate(non_overlap))]
                positions = positions[:-1]
                sequence = list(read.query_sequence)
                for k in range(len(positions) - 1):
                    element = positions[k]
                    sequence[element] = names[k]
                reads = "".join(sequence)
                cpg_all = [m.start() for m in re.finditer(r'CG', reads, re.IGNORECASE)]
                chg_all = [m.start() for m in re.finditer(r'C(A|C|T)G', reads, re.IGNORECASE)]
                chh_all = [m.start() for m in re.finditer(r'C(A|C|T)(A|T|C)', reads, re.IGNORECASE)]
                cpg_mods = [x for x in positions if x in cpg_all]
                chg_mods = [x for x in positions if x in chg_all]
                chh_mods = [x for x in positions if x in chh_all]
                for i in range(0,len(reads)):
                    if i in chh_mods:
                        read2_mods_CHH[i] += 1
                    if i in chg_mods:
                        read2_mods_CHG[i] += 1
                    if i in cpg_mods:
                        read2_mods_CpG[i] += 1
                for i in range(0,len(reads)):
                    if i in chh_all:
                        read2_all_CHH[i] += 1
                    if i in chg_all:
                        read2_all_CHG[i] += 1
                    if i in cpg_all:
                        read2_all_CpG[i] += 1
            else:
                sequence = list(read.query_sequence)
                reads = "".join(sequence)
                cpg_all = [m.start() for m in re.finditer(r'CG', reads, re.IGNORECASE)]
                chg_all = [m.start() for m in re.finditer(r'C(A|C|T)G', reads, re.IGNORECASE)]
                chh_all = [m.start() for m in re.finditer(r'C(A|C|T)(A|T|C)', reads, re.IGNORECASE)]
                for i in range(0,len(reads)):
                    if i in chh_all:
                        read2_all_CHH[i] += 1
                    if i in chg_all:
                        read2_all_CHG[i] += 1
                    if i in cpg_all:
                        read2_all_CpG[i] += 1

    for i in range(0, read_length):
        read_values_1_CHH[i] = [non_zero_div(read1_mods_CHH[i], read1_all_CHH[i])]
        read_values_1_CHG[i] = [non_zero_div(read1_mods_CHG[i], read1_all_CHG[i])]
        read_values_1_CpG[i] = [non_zero_div(read1_mods_CpG[i], read1_all_CpG[i])]
        values_1_CHH = [(keys + 1, round(values[0], 3)) for keys, values in read_values_1_CHH.items()]
        values_1_CHG = [(keys + 1, round(values[0], 3)) for keys, values in read_values_1_CHG.items()]
        values_1_CpG = [(keys + 1, round(values[0], 3)) for keys, values in read_values_1_CpG.items()]
        read_values_2_CHH[i] = [non_zero_div(read2_mods_CHH[i], read2_all_CHH[i])]
        read_values_2_CHG[i] = [non_zero_div(read2_mods_CHG[i], read2_all_CHG[i])]
        read_values_2_CpG[i] = [non_zero_div(read2_mods_CpG[i], read2_all_CpG[i])]
        values_2_CHH = [(keys + 1, round(values[0], 3)) for keys, values in read_values_2_CHH.items()]
        values_2_CHG = [(keys + 1, round(values[0], 3)) for keys, values in read_values_2_CHG.items()]
        values_2_CpG = [(keys + 1, round(values[0], 3)) for keys, values in read_values_2_CpG.items()]

    all_values = [(*a1, *a2, *a3, *a4, *a5, *a6) for a1, a2, a3, a4, a5, a6 in itertools.zip_longest(values_1_CpG, values_2_CpG, values_1_CHG, values_2_CHG, values_1_CHH, values_2_CHH)]
    all_values = [(x[0], x[1], x[3], x[5], x[7], x[9], x[11]) for x in all_values]
    with open(name + ".Mbias.txt", 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        wr.writerow(['POSITION (bp)', 'CpG READ 1', 'CpG READ 2', 'CHG READ 1', 'CHG READ 2', 'CHH READ 1', 'CHH READ 2'])
        for row in all_values:
            wr.writerow(row)

    x_axis = []
    y_axis_CpG1 = []
    y_axis_CHG1 = []
    y_axis_CHH1 = []
    y_axis_CpG2 = []
    y_axis_CHG2 = []
    y_axis_CHH2 = []

    for row in values_1_CpG:
        x_axis.append(row[0])
        y_axis_CpG1.append(row[1])

    for row in values_1_CHG:
        y_axis_CHG1.append(row[1])

    for row in values_1_CHH:
        y_axis_CHH1.append(row[1])

    for row in values_2_CpG:
        y_axis_CpG2.append(row[1])

    for row in values_2_CHG:
        y_axis_CHG2.append(row[1])

    for row in values_2_CHH:
        y_axis_CHH2.append(row[1])

    pyp.figure()
    pyp.plot(x_axis, y_axis_CpG1, linewidth=1.0, linestyle='-', color='teal')
    pyp.plot(x_axis, y_axis_CHG1, linewidth=1.0, linestyle='-', color='gray')
    pyp.plot(x_axis, y_axis_CHH1, linewidth=1.0, linestyle='-', color='maroon')
    pyp.xticks(numpy.arange(0, max(x_axis) + 1, step = max(x_axis)/10))
    pyp.suptitle('Read 1 M-bias', fontsize=14)
    pyp.xlabel('Read position, bp', fontsize=12)
    pyp.ylabel('Cytosine modification level', fontsize=12)
    pyp.legend(['CpG', 'CHG', 'CHH'], loc='center left', bbox_to_anchor=(1, 0.5))
    pyp.savefig(name + '_read1_M-bias_plot.pdf', figsize=(12, 8), dpi=330, bbox_inches='tight')
    pyp.close()

    pyp.figure()
    pyp.plot(x_axis, y_axis_CpG2, linewidth=1.0, linestyle='-', color='teal')
    pyp.plot(x_axis, y_axis_CHG2, linewidth=1.0, linestyle='-', color='gray')
    pyp.plot(x_axis, y_axis_CHH2, linewidth=1.0, linestyle='-', color='maroon')
    pyp.xticks(numpy.arange(0, max(x_axis) + 1, step = max(x_axis)/10))
    pyp.suptitle('Read 2 M-bias', fontsize=14)
    pyp.xlabel('Read position, bp', fontsize=12)
    pyp.ylabel('Cytosine modification level', fontsize=12)
    pyp.legend(['CpG', 'CHG', 'CHH'], loc='center left', bbox_to_anchor=(1, 0.5))
    pyp.savefig(name + '_read2_M-bias_plot.pdf', figsize=(12, 8), dpi=330, bbox_inches='tight')
    pyp.close()

    inbam.close()

read_values_1_CHH = {}
read_values_1_CHG = {}
read_values_1_CpG = {}
read_values_2_CHH = {}
read_values_2_CHG = {}
read_values_2_CpG = {}

read1_all_CHH = {}
read1_all_CHG = {}
read1_all_CpG = {}
read1_mods_CHH = {}
read1_mods_CHG = {}
read1_mods_CpG = {}
read2_all_CHH = {}
read2_all_CHG = {}
read2_all_CpG = {}
read2_mods_CHH = {}
read2_mods_CHG = {}
read2_mods_CpG = {}

regs = "(?:.*C.*)"

def non_zero_div(x, y):
    if y == 0:
        return 0
    else:
        return x / y



if __name__ == '__main__':
    Mbias_exec()







  