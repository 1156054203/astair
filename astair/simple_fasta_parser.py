import re
import sys
import logging

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)

def fasta_operator(parser, fasta_sequence, per_chromosome):
    """Parses fasta files with multiple genomes."""
    to_look_for_key = re.compile(r"(?<=>).*(?={})".format(parser))
    to_look_for_string = re.compile(r"(?<={})(?!>).*(?={})".format(parser, parser))
    keys = re.findall(to_look_for_key, fasta_sequence)
    sequences = re.findall(to_look_for_string, fasta_sequence)
    keys_indices = list()
    keys_indices.insert(0,0)
    for index in re.finditer(to_look_for_key, fasta_sequence):
        keys_indices.append(index.end())
    if len(keys) < len(sequences):
        sequences = list()
        for i in range(0,len(keys_indices)-1):
            new_strings = re.findall(to_look_for_string, fasta_sequence[keys_indices[i]:keys_indices[i+1]])
            joined_sequence = "".join(new_strings)
            sequences.append(joined_sequence)
        new_strings = re.findall(to_look_for_string, 
                                    fasta_sequence[keys_indices[-1]:])
        joined_sequence = "".join(new_strings)
        sequences.append(joined_sequence)
        sequences.remove('')
    if per_chromosome != None:
        keys_index = keys.index(per_chromosome)
        keys = keys[keys_index]
        sequences = sequences[keys_index]
    return keys, sequences

def fasta_splitting_by_sequence(fasta_file, per_chromosome):
    """Checks line termination and enables parsing of fasta files with multiple genomes."""
    fastas = {}
    keys, sequences = list(), list()
    try:
        with open(fasta_file, 'r') as fasta_handle:
            fasta_sequence = fasta_handle.read()
            if re.match(r".*(?=\r\n)", fasta_sequence):
                keys, sequences = fasta_operator("\r\n", fasta_sequence, per_chromosome)
            elif re.match(r".*(?=\n\r)", fasta_sequence):
                keys, sequences = fasta_operator("\n\r", fasta_sequence, per_chromosome)
            elif re.match(r".*(?=\n)", fasta_sequence):
                keys, sequences = fasta_operator("\n", fasta_sequence, per_chromosome)
            else:
                logs.error('The new line symbols in the presented fasta file do not match expected use cases. Please, change them to \n.')
        if per_chromosome == None:
            for i in range(0, len(keys)):
                fastas[keys[i]] = sequences[i]
        else:
            fastas[keys] = sequences
        return keys, fastas
    except Exception:
        logs.error('The genome reference fasta file does not exist.', exc_info=True)
        raise
