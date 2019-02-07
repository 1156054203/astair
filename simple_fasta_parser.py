import re
import sys
import logging


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
                    new_strings = [string for string in sequences if len(string) <= 75]
                    sequences = [string for string in sequences if string not in new_strings]
                    joined_sequence = "".join(new_strings)
                    sequences.insert(0, joined_sequence)
            elif re.match(r".*(?=\n)", fasta_sequence):
                keys = re.findall(r"(?<=>).*(?=\n)", fasta_sequence)
                sequences = re.findall(r"(?<=\n)(?!>).*(?=\n)", fasta_sequence)
        for i in range(0, len(keys)):
            fastas[keys[i]] = sequences[i]
        return keys, fastas
    except (SystemExit, KeyboardInterrupt, IOError, FileNotFoundError):
        logs.error('The genome reference fasta file does not exist.', exc_info=True)
        sys.exit(1)
