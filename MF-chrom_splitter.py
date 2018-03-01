import click
import pysam
from Bio import SeqIO



@click.command()
@click.argument('input_file')
@click.argument('fasta_file')
@click.argument('output_file')
def filter_exec(input_file, fasta_file, output_file):
    filter(input_file, fasta_file, output_file)

def filter(input_file, fasta_file, output_file):
    '''input_file is full name .bam file
    fasta_file is the fasta used
    output_file is file name without extension
    '''
    try:
        
key=[]
with open(fasta_file, "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        key.append(record.id)
   
     s
inbam = pysam.AlignmentFile(input_file, "rb")
for k in range(len(key)):
    bam_fetch = inbam.fetch(key[k])
    data = pysam.AlignmentFile(output_file+key[k]+".bam", "wb", template = inbam)
    for read in bam_fetch:
        data.write(read)
    
    finally:       
        inbam.close()



if __name__ == '__main__':
    filter_exec()
