import re
import sys
import pdb
import logging

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)


def fasta_splitting_by_sequence(fasta_file, per_chromosome):
    """Reads the reference line by line, which enables parsing of fasta files with multiple genomes."""
    fastas = {}
    keys, sequences, sequences_per_chrom = list(), list(), list()
    try:
        with open(fasta_file, 'r') as fasta_handle:
            for fasta_sequence in fasta_handle.readlines():
                if per_chromosome == None:
                    if re.match(r'^>', fasta_sequence.splitlines()[0]):
                        keys.append(fasta_sequence.splitlines()[0][1:])
                        sequences.append("".join(sequences_per_chrom))
                        sequences_per_chrom = list()
                    else:
                        sequences_per_chrom.append(fasta_sequence.splitlines()[0])
                else:
                    if re.match(r'^>', fasta_sequence.splitlines()[0]) and fasta_sequence.splitlines()[0][1:] == per_chromosome:
                        keys = fasta_sequence.splitlines()[0][1:]
                        chromosome_found = True
                    elif re.match(r'^>', fasta_sequence.splitlines()[0]) and fasta_sequence.splitlines()[0][1:] != per_chromosome:
                        chromosome_found = False
                        pass
                    else:
                        if chromosome_found == True:
                            sequences_per_chrom.append(fasta_sequence.splitlines()[0])
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
        return keys, fastas
    except Exception:
        logs.error('The genome reference fasta file does not exist.', exc_info=True)
        raise
