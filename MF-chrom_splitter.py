import click
import pysam
from Bio import SeqIO
import re


@click.command()
@click.option('input_file', '--input_file', required=True, help='BAM format file containing sequencing reads.')
@click.option('fasta_file', '--fasta_file', required=True, help='DNA sequence in fasta format used for aligning the sequencing reads and mpileup.')
def filter_exec(input_file, fasta_file):
    filter(input_file, fasta_file)

def filter(input_file, fasta_file):
    
    name = re.search(r'^.*(?=.bam)', input_file).group()
        
    key=[]
    with open(fasta_file, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            key.append(record.id)



    inbam = pysam.AlignmentFile(name + ".bam", "rb")
    for k in range(len(key)):
        bam_fetch = inbam.fetch(key[k])
        data = pysam.AlignmentFile(name +'.'+key[k]+".bam", "wb", template = inbam)
        for read in bam_fetch:
            data.write(read)
        
        
    inbam.close()



if __name__ == '__main__':
    filter_exec()
