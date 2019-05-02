import re
import sys
import pdb
import gzip
import logging

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)


def fasta_splitting_by_sequence(fasta_file, per_chromosome, write):
    """Reads the reference line by line, which enables parsing of fasta files with multiple genomes."""
    fastas = {}
    keys, sequences, sequences_per_chrom = list(), list(), list()
    try:
        if isinstance(fasta_file, str) and fasta_file[-3:] == '.gz':
            fasta_handle = gzip.open(fasta_file, 'rt')
        else:
            fasta_handle = open(fasta_file, 'r')
        spaces = False
        for fasta_sequence in fasta_handle.readlines():
            if per_chromosome == None:
                if re.match(r'^>', fasta_sequence.splitlines()[0]):
                    if fasta_sequence.splitlines()[0][1:].rfind(' ') == -1:
                        keys.append(fasta_sequence.splitlines()[0][1:])
                    else:
                        logs.info("There are spaces in the sequence names of your reference. asTair will replace them with underscores and output a new fasta file recommended for future analyses, or will run analyses with the first word of the reference names.")
                        spaces = True
                        keys.append(fasta_sequence.splitlines()[0][1:].split(' ')[0])
                    sequences.append("".join(sequences_per_chrom))
                    sequences_per_chrom = list()
                else:
                    sequences_per_chrom.append(fasta_sequence.splitlines()[0])
            else:
                if re.match(r'^>', fasta_sequence.splitlines()[0]) and fasta_sequence.splitlines()[0][1:] == per_chromosome:
                    if fasta_sequence.splitlines()[0][1:].rfind(' ') == -1:
                        keys = fasta_sequence.splitlines()[0][1:]
                    else:
                        logs.info("There are spaces in the sequence names of your reference. asTair will replace them with underscores and output a new fasta file recommended for future analyses, or will run analyses with the first word of the reference names.")
                        spaces = True
                        keys.append(fasta_sequence.splitlines()[0][1:].split(' ')[0])
                    chromosome_found = True
                elif re.match(r'^>', fasta_sequence.splitlines()[0]) and fasta_sequence.splitlines()[0][1:] != per_chromosome:
                    chromosome_found = False
                    pass
                else:
                    if chromosome_found == True:
                        sequences_per_chrom.append(fasta_sequence.splitlines()[0])
        if spaces == True and write == 'w':
            if isinstance(fasta_file, str) and fasta_file[-3:] == '.gz':
                fasta_handle = gzip.open(fasta_file, 'rt')
            else:
                fasta_handle = open(fasta_file, 'r')
            data_line = gzip.open(fasta_file[:-3] + '_no_spaces.fa.gz', 'wt')
            for fasta_sequence in fasta_handle.readlines():
                if re.match(r'^>', fasta_sequence.splitlines()[0]):
                        data_line.write('{}\n'.format(fasta_sequence.splitlines()[0].replace(' ', '_')))
                else:
                    data_line.write('{}\n'.format(fasta_sequence.splitlines()[0]))
            data_line.close()
        if per_chromosome == None:
            sequences.append("".join(sequences_per_chrom))
            sequences = sequences[1:]
            for i in range(0, len(keys)):
                fastas[keys[i]] = sequences[i]
        else:
            try:
                sequences = "".join(sequences_per_chrom)
                fastas[keys] = sequences
            except Exception:
                logs.error('The chromosome does not exist in the genome reference fasta file.', exc_info=True)
        fasta_handle.close()
        return keys, fastas
    except Exception:
        logs.error('The genome reference fasta file does not exist.', exc_info=True)
        raise
