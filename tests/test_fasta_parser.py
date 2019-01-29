import unittest
from io import StringIO
import csv



class FastaParserTest(unittest.TestCase):
    def test_parsing_fasta_files_single(self):
        mock_fasta = StringIO(">some_fasta_sequence\nACTGCTCCCTGGaaaTCG\n")
        keys, sequences = list(), list()
        wr = csv.reader(mock_fasta, delimiter='\t', lineterminator='\n')
        for n in wr:
            if n[0][0] == '>':
                keys.append(n[0][1:])
            else:
                sequences.append(n[0])
        self.assertEqual(keys, ['some_fasta_sequence'])
        self.assertEqual(sequences, ['ACTGCTCCCTGGaaaTCG'])
    def test_parsing_fasta_files_multiple(self):
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


if __name__ == '__main__':
    unittest.main()