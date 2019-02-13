import pysam
import sys
import logging

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

def bam_file_opener(input_file, fetch):
    """Opens neatly and separately the bam file as an iterator."""
    try:
        open(input_file, 'rb')
        inbam = pysam.AlignmentFile(input_file, "rb", header=True, threads=10)
        if isinstance(fetch, str):
            bam_fetch = inbam.fetch(until_eof=True)
            return bam_fetch
        elif isinstance(fetch, tuple):
            bam_fetch = inbam.fetch(fetch[0], fetch[1], fetch[2])
            return bam_fetch
        else:
            return inbam
    except (SystemExit, KeyboardInterrupt, IOError, FileNotFoundError):
        logs.error('The input bam file does not exist.', exc_info=True)
        sys.exit(1)
