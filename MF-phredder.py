import click
import itertools
import csv
from Bio import SeqIO
import gzip


@click.command()
@click.argument('input_file')
def phredder_exec(input_file):
    phredder(input_file)

def phredder(input_file):
    '''input_file is file name without extension
    '''


    Phred_T_sum = 0
    Phred_T_len = 0
    Phred_C_sum = 0
    Phred_C_len = 0
    Phred_A_sum = 0
    Phred_A_len = 0
    Phred_G_sum = 0
    Phred_G_len = 0


    def non_zero_div(x, y):
        if y == 0:
            return "NA"
        elif x == 0:
            return 0
        else:
            return x / y

    try:


        print("Lets begin.")

        cycle_count = 0
        with gzip.open(input_file + ".fastq.gz", "rt") as handle:
            print("Open file.")
            record_iter = SeqIO.parse(handle, "fastq")
            print("Iterator ready.")
            for line in record_iter:
                results = [(i, *j) for i, j in itertools.zip_longest(line.letter_annotations["phred_quality"], line.seq.tostring(), fillvalue='BLANK')]
                thymines = [x[0] for x in results if 'T' in x[:]]
                cytosines = [x[0] for x in results if 'C' in x[:]]
                adenines = [x[0] for x in results if 'A' in x[:]]
                guanines = [x[0] for x in results if 'G' in x[:]]
                Phred_T_sum = sum(thymines) + Phred_T_sum
                Phred_T_len = len(thymines) + Phred_T_len
                Phred_C_sum = sum(cytosines) + Phred_C_sum
                Phred_C_len = len(cytosines) + Phred_C_len
                Phred_A_sum = sum(adenines) + Phred_A_sum
                Phred_A_len = len(adenines) + Phred_A_len
                Phred_G_sum = sum(guanines) + Phred_G_sum
                Phred_G_len = len(guanines) + Phred_G_len
                read_values = [non_zero_div(sum(thymines), len(thymines)), non_zero_div(sum(cytosines), len(cytosines)), non_zero_div(sum(adenines), len(adenines)), non_zero_div(sum(guanines), len(guanines))]
                with open(input_file + ".Phred.txt", 'a', newline='') as myfile:
                    wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
                    wr.writerow(read_values)
                if cycle_count%100000==0:
                    print("100000 more reads...hold on. Intermediate results: {}".format([non_zero_div(Phred_T_sum, Phred_T_len), non_zero_div(Phred_C_sum, Phred_C_len), non_zero_div(Phred_A_sum, Phred_A_len), non_zero_div(Phred_G_sum, Phred_G_len)]))
                cycle_count += 1

        total_phred = [non_zero_div(Phred_T_sum, Phred_T_len), non_zero_div(Phred_C_sum, Phred_C_len), non_zero_div(Phred_A_sum, Phred_A_len), non_zero_div(Phred_G_sum, Phred_G_len)]

        print(total_phred)

        with open(input_file + '_total_Phred.txt', 'w', newline='') as myfile:
            wr = csv.writer(myfile, delimiter='\t', lineterminator='\n')
            wr.writerow(["Thymine", "Cytosine", "Adenine", "Guanine"])
            wr.writerow(total_phred)



    except OSError:
        print('No fastq file found.')


if __name__ == '__main__':
    phredder_exec()
