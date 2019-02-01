import unittest
from io import StringIO
import re
import csv


class FastaParserTest(unittest.TestCase):
    """Tests whether fasta-like strings can be split meaningfully into DNA strings and chromosome names (keys).
     In case \r\n line terminators are discovered together with having more DNA strings than keys, a naive
     concatenation procedure of the multiple short DNA strings is performed."""
    def test_parsing_fasta_files_single(self):
        mock_fasta = ">some_fasta_sequence\nACTGCTCCCTGGaaaTCG\n"
        keys, sequences = list(), list()
        if re.match(r".*(?=\r\n)", mock_fasta):
            keys = re.findall(r"(?<=>).*(?=\r\n)", mock_fasta)
            sequences = re.findall(r"(?<=\r\n)(?!>).*(?=\r\n)", mock_fasta)
            if len(keys) < len(sequences):
                new_strings = [string for string in sequences if len(string) <= 75]
                not_to_merge = [string for string in sequences if string not in new_strings]
                joined_sequence = "".join(new_strings)
                sequences = list()
                sequences.append(joined_sequence)
                sequences.append(not_to_merge)
        elif re.match(r".*(?=\n)", mock_fasta):
            keys = re.findall(r"(?<=>).*(?=\n)", mock_fasta)
            sequences = re.findall(r"(?<=\n)(?!>).*(?=\n)", mock_fasta)
        self.assertEqual(keys, ['some_fasta_sequence'])
        self.assertEqual(sequences, ['ACTGCTCCCTGGaaaTCG'])
    def test_parsing_fasta_files_multiple(self):
        mock_fasta = ">some_fasta_sequence\nACTGCTCCCTGGaaaTCG\n>yet_another_fasta_sequence\nAAACCTGCcctGttug\n>fasta_sequence_again\nAAAAAAACCTGCTAGctaatat\n>and_again_sequence\nCTGATCGTTTAGCAGCA\n"
        keys, sequences = list(), list()
        if re.match(r".*(?=\r\n)", mock_fasta):
            keys = re.findall(r"(?<=>).*(?=\r\n)", mock_fasta)
            sequences = re.findall(r"(?<=\r\n)(?!>).*(?=\r\n)", mock_fasta)
            if len(keys) < len(sequences):
                new_strings = [string for string in sequences if len(string) <= 75]
                not_to_merge = [string for string in sequences if string not in new_strings]
                joined_sequence = "".join(new_strings)
                sequences = list()
                sequences.append(joined_sequence)
                sequences.append(not_to_merge)
        elif re.match(r".*(?=\n)", mock_fasta):
            keys = re.findall(r"(?<=>).*(?=\n)", mock_fasta)
            sequences = re.findall(r"(?<=\n)(?!>).*(?=\n)", mock_fasta)
        self.assertEqual(keys, ['some_fasta_sequence', 'yet_another_fasta_sequence', 'fasta_sequence_again', 'and_again_sequence'])
        self.assertEqual(sequences, ['ACTGCTCCCTGGaaaTCG', 'AAACCTGCcctGttug', 'AAAAAAACCTGCTAGctaatat', 'CTGATCGTTTAGCAGCA'])
    def test_parsing_fasta_files_multiple_old(self):
        mock_fasta = StringIO(">some_fasta_sequence\nACTGCTCCCTGGaaaTCG\n>yet_another_fasta_sequence\nAAACCTGCcctGttug\n>fasta_sequence_again\nAAAAAAACCTGCTAGctaatat\n>and_again_sequence\nCTGATCGTTTAGCAGCA\n")
        keys, sequences = list(), list()
        wr = csv.reader(mock_fasta, delimiter='\t', lineterminator='\n')
        for n in wr:
            if n[0][0] == '>':
                keys.append(n[0][1:])
            else:
                sequences.append(n[0])
        self.assertEqual(keys, ['some_fasta_sequence', 'yet_another_fasta_sequence', 'fasta_sequence_again', 'and_again_sequence'])
        self.assertEqual(sequences, ['ACTGCTCCCTGGaaaTCG', 'AAACCTGCcctGttug', 'AAAAAAACCTGCTAGctaatat', 'CTGATCGTTTAGCAGCA'])
    def test_parsing_fasta_files_multiple_rn_line_terminators(self):
        mock_fasta = ">some_fasta_sequence\r\nGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG\r\nGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG" \
                     "\r\n>fasta_sequence_again\r\nAAAAAAACCTGCTAGctaatatGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG\r\n>and_again_sequence\r\n" \
                     "CTGATCGTTTAGCAGCGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGA\r\n"
        keys, sequences = list(), list()
        if re.match(r".*(?=\r\n)", mock_fasta):
            keys = re.findall(r"(?<=>).*(?=\r\n)", mock_fasta)
            sequences = re.findall(r"(?<=\r\n)(?!>).*(?=\r\n)", mock_fasta)
            if len(keys) < len(sequences):
                new_strings = [string for string in sequences if len(string) <= 75]
                not_to_merge = [string for string in sequences if string not in new_strings]
                joined_sequence = "".join(new_strings)
                sequences = list()
                sequences.append(joined_sequence)
                sequences.extend(not_to_merge)
        elif re.match(r".*(?=\n)", mock_fasta):
            keys = re.findall(r"(?<=>).*(?=\n)", mock_fasta)
            sequences = re.findall(r"(?<=\n)(?!>).*(?=\n)", mock_fasta)
        self.assertEqual(keys, ['some_fasta_sequence', 'fasta_sequence_again', 'and_again_sequence'])
        self.assertEqual(sequences,['GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG',
                                    'AAAAAAACCTGCTAGctaatatGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG',
                                    'CTGATCGTTTAGCAGCGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGA'])


if __name__ == '__main__':
    unittest.main()