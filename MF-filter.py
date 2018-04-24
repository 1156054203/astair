import click
import pysam

@click.command()
@click.argument('input_file')
@click.argument('output_file')
def filter_exec(input_file, output_file):
    filter(input_file, output_file)

def filter(input_file, output_file):
    '''input_file is full name .bam file
    output_file is file name without extension
    '''
    try:
        inbam = pysam.AlignmentFile(input_file, "rb")
        bam_fetch = inbam.fetch(until_eof=True)
        outbamOT = pysam.AlignmentFile(output_file+"OT.bam", "wb", template = inbam)
        outbamOB = pysam.AlignmentFile(output_file+"OB.bam", "wb", template = inbam) 
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
