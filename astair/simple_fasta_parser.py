import re
import os 
import sys
import pdb
import gzip
import numpy
import pysam
import logging
import subprocess

logging.basicConfig(level=logging.DEBUG)
logs = logging.getLogger(__name__)



def gzipped_fasta_read(fasta_file, reference_absolute_name, reference_extension):
    "Reads through a GZIP compressed fasta file."
    if reference_extension == '.gz' and sys.version[0] == '3':
        fasta_handle = gzip.open(fasta_file, 'rt')
    elif reference_extension == '.gz' and sys.version[0] == '2':
        if not os.path.isfile(reference_absolute_name):
            file_ = subprocess.Popen('gunzip {}'.format(fasta_file), shell=True)
            exit_code = file_.wait()
            if exit_code == 0:
                fasta_handle = open(reference_absolute_name, 'r')
        else:
            fasta_handle = open(reference_absolute_name, 'r')
    return fasta_handle


def output_fasta_wuth_underscores(data_line, reference_extension, fasta_file, reference_absolute_name, reference_dir):
    """Outputs a FASTA file with underscores in the names."""
    if reference_extension == '.gz':
        fasta_handle = gzipped_fasta_read(fasta_file, reference_absolute_name, reference_extension)
    else:
        fasta_handle = open(fasta_file, 'r')
    if not os.path.isfile(reference_dir + os.path.splitext(os.path.basename(reference_absolute_name))[0] + '_no_spaces.fa.gz'):
        data_line = open(reference_dir + os.path.splitext(os.path.basename(reference_absolute_name))[0] + '_no_spaces.fa', 'w')
        for fasta_sequence in fasta_handle.readlines():
            if re.match(r'^>', fasta_sequence.splitlines()[0]):
                    data_line.write('{}\n'.format(fasta_sequence.splitlines()[0].replace(' ', '_')))
            else:
                data_line.write('{}\n'.format(fasta_sequence.splitlines()[0]))
        try:
            subprocess.Popen('bgzip {}'.format(reference_dir + os.path.splitext(os.path.basename(reference_absolute_name))[0] + '_no_spaces.fa'), shell=True)
            logs.info("The program will ouput a BGZIP compressed fasta files with underscores in the reference names for future analyses.")
        except Exception:
            subprocess.Popen('gzip {}'.format(reference_dir + os.path.splitext(os.path.basename(reference_absolute_name))[0] + '_no_spaces.fa'), shell=True)
            logs.info("The program will ouput a GZIP compressed fasta files with underscores in the reference names for future analyses.")
        data_line.close()
        fasta_handle.close()


def fasta_splitting_by_sequence(fasta_file, per_chromosome, numbered, add_underscores, all_chromosomes):
    """Reads the reference line by line, which enables parsing of fasta files with multiple genomes."""
    try:
        if (sys.version[0] == '3' and isinstance(fasta_file, str)) or (sys.version[0] == '2' and isinstance(fasta_file, basestring)):
            reference_absolute_name = os.path.splitext(os.path.abspath(fasta_file))[0]
            reference_extension = os.path.splitext(os.path.basename(fasta_file))[1]
            reference_dir = os.path.dirname(fasta_file)
            if list(reference_dir)[-1]!="/":
                reference_dir = reference_dir + "/"
            if reference_extension == '.gz':
                compressed_ = 'gzip'
                try:
                    test = pysam.FastaFile(fasta_file)
                    compressed_ = 'bgzip'
                except Exception:
                    logs.error('The reference FASTA file was not compressed with BGZIP or does not have an index.', exc_info=True)
            if add_underscores:
                output_fasta_wuth_underscores(data_line, reference_extension, fasta_file, reference_absolute_name, reference_dir)
            if reference_extension != '.gz' or compressed_=='bgzip':
                try:
                    keys, fastas, sequences, sequences_per_chrom = numpy.array([]), {}, numpy.array([]), numpy.array([])
                    all_chrom = pysam.FastaFile(fasta_file)
                    keys = all_chrom.references
                    if all_chromosomes is None and per_chromosome == 'keys_only':
                        return keys
                    elif all_chromosomes is None and per_chromosome != 'keys_only':
                        sequences_per_chrom = pysam.FastaFile(fasta_file).fetch(per_chromosome)
                        sequences = "".join(sequences_per_chrom)
                        fastas[per_chromosome] = sequences
                        return fastas
                    else:
                        for sequence_name in keys:
                            fastas[sequence_name] = pysam.FastaFile(fasta_file).fetch(sequence_name)
                        return keys, fastas
                except Exception:
                    logs.error('The chromosome does not exist in the genome reference fasta file or the FASTA file is not indexed.', exc_info=True)
            else:
                chromosome_found = False
                keys, fastas, sequences, sequences_per_chrom = [], {}, [], []
                fasta_handle = gzipped_fasta_read(fasta_file, reference_absolute_name, reference_extension)
                for fasta_sequence in fasta_handle.readlines():
                        if re.match(r'^>', fasta_sequence.splitlines()[0]):
                            if fasta_sequence.splitlines()[0][1:].rfind(' ') == -1:
                                keys.append(fasta_sequence.splitlines()[0][1:])
                            else:
                                logs.info("There are spaces in the sequence names of your reference. Please add them yourself or run asTair with --add_underscores option, which will replace them with underscores and output a new fasta file recommended for future analyses, now it will run analyses with the first word of the reference names, unless you have also set --use_underscores.")
                                keys.append(fasta_sequence.splitlines()[0][1:].split(' ')[0])
                            if all_chromosomes is None and per_chromosome == 'keys_only':
                                fasta_handle.close()
                                if reference_extension == '.gz' and sys.version[0] == '2' and  os.path.isfile(reference_absolute_name):
                                    file_ = subprocess.Popen('gzip {}'.format(reference_absolute_name), shell=True)
                                    exit_code = file_.wait()
                                return keys
                            elif per_chromosome is None:
                                sequences.append("".join(sequences_per_chrom))
                                sequences_per_chrom = list()
                        else:
                            if per_chromosome is not None:
                                if keys[-1] == per_chromosome:
                                    sequences_per_chrom.append(fasta_sequence.splitlines()[0])
                                    sequences = "".join(sequences_per_chrom)
                                    fastas[per_chromosome] = sequences
                                    fasta_handle.close()
                                    if reference_extension == '.gz' and sys.version[0] == '2' and  os.path.isfile(reference_absolute_name):
                                        file_ = subprocess.Popen('gzip {}'.format(reference_absolute_name), shell=True)
                                        exit_code = file_.wait()
                                    if all_chromosomes is None:
                                        return fastas
                                    else:
                                        return keys, fastas
                            else:
                                sequences_per_chrom.append(fasta_sequence.splitlines()[0])
                if per_chromosome is None and all_chromosomes is not None:
                    sequences.append("".join(sequences_per_chrom))
                    sequences = sequences[1:]
                    for i in range(0, len(keys)):
                        fastas[keys[i]] = sequences[i]
                    fasta_handle.close()
                    if reference_extension == '.gz' and sys.version[0] == '2' and  os.path.isfile(reference_absolute_name) and numbered == "last":
                        file_ = subprocess.Popen('gzip {}'.format(reference_absolute_name), shell=True)
                        exit_code = file_.wait()
                    return keys, fastas
        else:
            keys, fastas, sequences_per_chrom, sequences = [], {}, [], []
            fasta_handle = open(fasta_file, 'r')
            for fasta_sequence in fasta_handle.readlines():
                if re.match(r'^>', fasta_sequence.splitlines()[0]):
                    keys.append(fasta_sequence.splitlines()[0][1:])
                    sequences.append("".join(sequences_per_chrom))
                    sequences_per_chrom = list()
                else:
                    sequences_per_chrom.append(fasta_sequence.splitlines()[0])
            sequences.append("".join(sequences_per_chrom))
            sequences = sequences[1:]
            for i in range(0, len(keys)):
                fastas[keys[i]] = sequences[i]
            fasta_handle.close()
            return keys, fastas
    except Exception:
        logs.error('The genome reference fasta file does not exist.', exc_info=True)
        raise
