import click
import pysam
import re


@click.command()
@click.argument('input_file')
@click.argument('output_file')
def filter_exec(input_file, output_file):
    filter(input_file, output_file)

def filter(input_file, output_file):
    '''input_file is full name .bam file
    output_filr is file name without extension
    '''

   # try:
        
    regs = re.compile("(?:.*C.*){3}", re.IGNORECASE)
    ch = re.compile("(C)(A|C|T)(A|T|C|G)", re.IGNORECASE)
    numbs = re.compile("(?:.*?(\d+)C.*?)(?:.*?(\d+)C.*?)(?:.*?(\d+)C.*?)", re.IGNORECASE)

    
    inbam = pysam.AlignmentFile(input_file, "rb")
    bam_fetch = inbam.fetch(until_eof=True)
    outbam3T = pysam.AlignmentFile(output_file+"3T.bam", "wb", template = inbam)
    for read in bam_fetch:
        if re.search(regs, read.tags[1][1]):
          non_cpg =  [m.start() for m in re.finditer(ch, read.query_sequence)]
          if non_cpg:
              mismatch = re.search(numbs, read.tags[1][1])
              if mismatch:
                  mismatch_pos = [int(mismatch.groups()[0]), int(mismatch.groups()[1]) + int(mismatch.groups()[0]) + 1, int(mismatch.groups()[2]) + int(mismatch.groups()[1]) + int(mismatch.groups()[0]) + 2]
                  if [x for x in mismatch_pos if x not in non_cpg]:
                    outbam3T.write(read)
              else:
                  outbam3T.write(read)
          else:
               outbam3T.write(read)
        else:
            outbam3T.write(read)
       
       
       
   # finally:       
    inbam.close()


if __name__ == '__main__':
    filter_exec()







  
