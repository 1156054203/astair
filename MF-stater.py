from Bio import SeqIO
import csv
import re


def non_zero_div(x, y):
    if y == 0:
        return 0
    else:
        return x / y




data_bedgraph = []
with open("/home/amygdala/Downloads/trial3/spikes/BS-free_spikes_bwa_L004-stotal.bedGraph") as bedgraph:
    bedgraph_file = csv.reader(bedgraph, delimiter='\t', lineterminator='\n')
    for row in bedgraph_file:
        data_bedgraph.append((row[0], int(row[1]), int(row[2]), float(row[3])))


key = []
test_seq = []
#with open("/media/gergana/DiscB1/Analysis/vcfs/mm9.fa", "rU") as handle:
with open("/home/amygdala/Videos/all-spike-ins-bisulfite-free.fa", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        key.append(record.id)
        test_seq.append(record)



for k in range(len(key)):
    seqs = test_seq[k].seq
    chars = []
    CpG = []
    CHG = []
    CHH = []
    HCG = []
    GCH = []
    for line in seqs:
            chars.append(line)
    fasta_seq = ''.join(chars)
    chh_pat = re.compile("C(A|C|T)(A|T|C)", re.IGNORECASE)
    for m in chh_pat.finditer(fasta_seq):
        CHH.append(m.start())
    chg_pat = re.compile("C(A|C|T)G", re.IGNORECASE)
    for m in chg_pat.finditer(fasta_seq):
        CHG.append(m.start())
    hcg_pat = re.compile("(A|C|T)CG", re.IGNORECASE)
    for m in hcg_pat.finditer(fasta_seq):
        HCG.append(m.start() + 1)
        HCG.append(m.start() + 2)
    gch_pat = re.compile("GC(A|C|T)", re.IGNORECASE)
    for m in gch_pat.finditer(fasta_seq):
        GCH.append(m.start())
        GCH.append(m.start() + 1)
    cpg_pat = re.compile("CG", re.IGNORECASE)
    for m in cpg_pat.finditer(fasta_seq):
        CpG.append(m.start())
        CpG.append(m.start() + 1)
    data_CHH = [(key[k], each, each + 1) for each in CHH]
    data_CHG = [(key[k], each, each + 1) for each in CHG]
    data_CpG = [(key[k], each, each + 1) for each in CpG]
    data_HCG = [(key[k], each, each + 1) for each in HCG]
    data_GCH = [(key[k], each, each + 1) for each in GCH]
    CHH_set = set(data_CHH)
    CHG_set = set(data_CHG)
    CpG_set = set(data_CpG)
    HCG_set = set(data_HCG)
    GCH_set = set(data_GCH)
    data_sample_CHH = list(filter(lambda x: (x[0], x[1], x[2]) in CHH_set, data_bedgraph))
    data_sample_CHG = list(filter(lambda x: (x[0], x[1], x[2]) in CHG_set, data_bedgraph))
    data_sample_CpG = list(filter(lambda x: (x[0], x[1], x[2]) in CpG_set, data_bedgraph))
    data_sample_HCG = list(filter(lambda x: (x[0], x[1], x[2]) in HCG_set, data_bedgraph))
    data_sample_GCH = list(filter(lambda x: (x[0], x[1], x[2]) in GCH_set, data_bedgraph))
    modification_rate_CHH = non_zero_div(sum([x[3] for x in data_sample_CHH]), len(data_sample_CHH))
    modification_rate_CHG = non_zero_div(sum([x[3] for x in data_sample_CHG]), len(data_sample_CHG))
    modification_rate_CpG = non_zero_div(sum([x[3] for x in data_sample_CpG]), len(data_sample_CpG))
    modification_rate_HCG = non_zero_div(sum([x[3] for x in data_sample_HCG]), len(data_sample_HCG))
    modification_rate_GCH = non_zero_div(sum([x[3] for x in data_sample_GCH]), len(data_sample_GCH))
    with open("BSF3_L004_CHH.bedGraph", 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        for row in data_sample_CHH:
            wr.writerow(row)
    with open("BSF3_L004_CHG.bedGraph", 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        for row in data_sample_CHG:
            wr.writerow(row)
    with open("BSF3_L004_CpG.bedGraph", 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        for row in data_sample_CpG:
            wr.writerow(row)
    with open("BSF3_L004_HCG.bedGraph", 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        for row in data_sample_HCG:
            wr.writerow(row)
    with open("BSF3_L004_GCH.bedGraph", 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        for row in data_sample_GCH:
            wr.writerow(row)
    with open("BSF3_L004_GCH.bedGraph", 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        wr.writerow([key[k], "GCH", modification_rate_GCH])
    with open("BSF3_L004.stats", 'a', newline='') as myfile:
        wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
        wr.writerow([key[k], "CpG", round(modification_rate_CpG,3), len(data_CpG)])
        wr.writerow([key[k], "CHH", round(modification_rate_CHH,3), len(data_CHH)])
        wr.writerow([key[k], "CHG", round(modification_rate_CHG,3), len(data_CHG)])
        wr.writerow([key[k], "GCH", round(modification_rate_GCH,3), len(data_GCH)])
        wr.writerow([key[k], "HCG", round(modification_rate_HCG,3), len(data_HCG)])

#modification_rate_CHH = non_zero_div(sum([x[3] for x in data_sample_CHH]), len(data_sample_CHH))
#pos_modification_rate_CHH = non_zero_div(sum([x[3] for x in data_sample_CHH if x[3] >= modification_rate_CHH]), len([x[3] for x in data_sample_CHH if x[3] >= modification_rate_CHH])




                    # with open("allchromCHH.bed", 'a', newline='') as myfile:
    #     wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
    #     for row in data_CHH:
    #         wr.writerow(row)
    # with open("allchromCHG.bed", 'a', newline='') as myfile:
    #     wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
    #     for row in data_CHG:
    #         wr.writerow(row)
    # with open("allchromCpG.bed", 'a', newline='') as myfile:
    #     wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
    #     for row in data_CpG:
    #         wr.writerow(row)
    # with open("allchromHCG.bed", 'a', newline='') as myfile:
    #     wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
    #     for row in data_HCG:
    #         wr.writerow(row)
    # with open("allchromGCH.bed", 'a', newline='') as myfile:
    #     wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
    #     for row in data_GCH:
    #         wr.writerow(row)



