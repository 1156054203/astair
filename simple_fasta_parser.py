import re
import sys
import logging
import pdb


logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)


def fasta_splitting_by_sequence(fasta_file):
    """Parses fasta files with multiple genomes."""
    fastas = {}
    keys, sequences = list(), list()
    try:
        with open(fasta_file, 'r', newline='') as fasta_handle:
            fasta_sequence = fasta_handle.read()
            if re.match(r".*(?=\r\n)", fasta_sequence):
                keys = re.findall(r"(?<=>).*(?=\r\n)", fasta_sequence)
                sequences = re.findall(r"(?<=\r\n)(?!>).*(?=\r\n)", fasta_sequence)
                if len(keys) < len(sequences):
                    for i in range(0,len(keys_indices)-1):
                        new_strings = re.findall(r"(?<=\n)(?!>).*(?=\n)", fasta_sequence[keys_indices[i]:keys_indices[i+1]])
                        joined_sequence = "".join(new_strings)
                        sequences.append(joined_sequence)
                    new_strings = re.findall(r"(?<=\n)(?!>).*(?=\n)", 
                                             fasta_sequence[keys_indices[-1]:])
                    joined_sequence = "".join(new_strings)
                    sequences.append(joined_sequence)
                    sequences.remove('')
            elif re.match(r".*(?=\n)", fasta_sequence):
                keys = re.findall(r"(?<=>).*(?=\n)", fasta_sequence)
                sequences = re.findall(r"(?<=\n)(?!>).*(?=\n)", fasta_sequence)
                keys_indices = list()
                keys_indices.insert(0,0)
                for index in re.finditer(r"(?<=>).*(?=\n)", fasta_sequence):
                    keys_indices.append(index.end())
                if len(keys) < len(sequences):
                    sequences = list()
                    for i in range(0,len(keys_indices)-1):
                        new_strings = re.findall(r"(?<=\n)(?!>).*(?=\n)", fasta_sequence[keys_indices[i]:keys_indices[i+1]])
                        joined_sequence = "".join(new_strings)
                        sequences.append(joined_sequence)
                    new_strings = re.findall(r"(?<=\n)(?!>).*(?=\n)", 
                                             fasta_sequence[keys_indices[-1]:])
                    joined_sequence = "".join(new_strings)
                    sequences.append(joined_sequence)
                    sequences.remove('')
        for i in range(0, len(keys)):
            fastas[keys[i]] = sequences[i]
        return keys, fastas
    except (SystemExit, KeyboardInterrupt, IOError, FileNotFoundError):
        logs.error('The genome reference fasta file does not exist.', exc_info=True)
        sys.exit(1)
