import os
import sys
import pysam
import logging

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

def bam_file_opener(input_file, fetch, threads):
    """Opens neatly and separately the bam file as an iterator."""
    try:
        reads_file = open(input_file, 'rb')
        reads_file.close()
        if (input_file[-4:] == '.bam' and not os.path.isfile(os.path.join(str(input_file)+'.bai'))) \
                or (input_file[-5:] == '.cram' and not os.path.isfile(os.path.join(str(input_file)+'.crai'))):
            try:
                pysam.index(input_file)
                logs.info('Building index for the input file.')
            except Exception:
                logs.error('The input file is not sorted by coordinates. Please sort before running again.', exc_info=True)
                raise
        inbam = pysam.AlignmentFile(input_file, "rb", header=True, threads=threads)
        if isinstance(fetch, str):
            bam_fetch = inbam.fetch(until_eof=True)
            return bam_fetch
        elif isinstance(fetch, tuple):
            bam_fetch = inbam.fetch(fetch[0], fetch[1], fetch[2])
            return bam_fetch
        else:
            return inbam
    except Exception:
        logs.error('The input bam file does not exist or is unsorted.', exc_info=True)
        raise
