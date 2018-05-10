import click
import pysam
import re

@click.command()
@click.option('input_file', '--input_file', required=True, help='BAM format file containing sequencing reads.')
def filter_exec(input_file):
    filter(input_file)

def filter(input_file):

    name = re.search(r'^.*(?=.bam)', input_file).group()
    
    try:
        inbam = pysam.AlignmentFile(input_file, "rb")
        bam_fetch = inbam.fetch(until_eof=True)
        outbamOT = pysam.AlignmentFile(name+"OT.bam", "wb", template = inbam)
        outbamOB = pysam.AlignmentFile(name+"OB.bam", "wb", template = inbam) 
        for read in bam_fetch:
            if read.flag == 99 or read.flag == 147:
                outbamOT.write(read)
            elif read.flag == 83 or read.flag == 163: 
                outbamOB.write(read)
    finally:       
        inbam.close()
        outbamOB.close()
        outbamOT.close() 


if __name__ == '__main__':
    filter_exec()
