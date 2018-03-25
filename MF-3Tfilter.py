import click
import pysam
import re
from multiprocessing import Process
import itertools


@click.command()
@click.argument('input_file')
@click.argument('bases_noncpg')
def filter3T_exec(input_file, bases_noncpg):
    filter3T(input_file, bases_noncpg)

def filter3T(input_file, bases_noncpg):
    '''input_file is file name without extension
    bases_noncpg is the number of non-cpg conversion events to consider the read for removal
    '''

# try:

    p1 = Process(target=single_strander, args=(input_file, "OT", bases_noncpg))
    p2 = Process(target=single_strander, args=(input_file, "OB", bases_noncpg))

    p1.start()
    p2.start()


def single_strander(input_file, strand, bases_noncpg):
    inbam = pysam.AlignmentFile(input_file+strand+".bam", "rb")
    bam_fetch = inbam.fetch(until_eof=True)
    outbam3T = pysam.AlignmentFile(input_file+"3T_"+strand+".bam", "wb", template=inbam)
    removed3T = pysam.AlignmentFile(input_file+"removed_" + strand+".bam", "wb", template=inbam)
    if strand == "OT":
        regs = "(?:.*C.*)" * int(bases_noncpg)
        ch = "(C)(A|C|T)(A|T|C|G)"
        ref = 'C'
    elif strand == "OB":
        regs = "(?:.*G.*)" * int(bases_noncpg)
        ch = "(G)(A|G|T)(A|T|C|G)"
        ref = 'G'
    for read in bam_fetch:
        if isinstance(read.tags[0][1],str):
            read_data = read.tags[0][1]
        elif isinstance(read.tags[1][1],str):
            read_data = read.tags[1][1]
            if re.search(regs, read_data,re.IGNORECASE):
                changes = [int(s) for s in re.findall(r'\d+', read_data)]
                non_overlap = [x + 1 if x == 0 else x for x in changes]
                names = list(re.findall(r'[^\W\d_]+', read_data))
                positions = [x +1 for x in list(itertools.accumulate(non_overlap))]
                positions = positions[:-1]
                sequence = list(read.query_sequence)
                for k in range(len(positions)-1):
                    element = positions[k]
                    sequence[element] = names[k]
                reads = "".join(sequence)
                non_cpg = [m.start() for m in re.finditer(ch, reads, re.IGNORECASE)]
                index = [(i, *j) for i, j in itertools.zip_longest(positions, names, fillvalue='BLANK')]
                if non_cpg:
                    mismatch = [x[0] for x in index if x[1] == ref]
                    if len(mismatch) >= int(bases_noncpg):
                       if mismatch not in non_cpg:
                            outbam3T.write(read)
                       else:
                            removed3T.write(read)
                    else:
                        outbam3T.write(read)
                else:
                    outbam3T.write(read)
            else:
                outbam3T.write(read)


    inbam.close()



if __name__ == '__main__':
    filter3T_exec()







  
