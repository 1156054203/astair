import unittest
from unittest.mock import patch, mock_open
from simple_fasta_parser import fasta_splitting_by_sequence

class FastaParserTest(unittest.TestCase):
    """Tests whether fasta-like strings can be split meaningfully into DNA strings and chromosome names (keys).
     In case \r\n line terminators are discovered together with having more DNA strings than keys, a naive
     concatenation procedure of the multiple short DNA strings is performed."""

    def test_parsing_fast_files_single(self):
        """Tests whether a short single fasta-like string can be split into DNA sequence and chromosome name."""
        with patch("builtins.open", mock_open(read_data=">some_fasta_sequence\nACTGCTCCCTGGaaaTCG\n")) as mock_file:
            keys, fastas = fasta_splitting_by_sequence(mock_file)
            self.assertEqual(keys, ['some_fasta_sequence'])
            self.assertEqual([fastas[key] for key in keys], ['ACTGCTCCCTGGaaaTCG'])

    def test_parsing_fasta_files_multiple(self):
        """Tests whether several short fasta-like strings can be split into DNA sequences and chromosome names."""
        with patch("builtins.open", mock_open(read_data=">some_fasta_sequence\nACTGCTCCCTGGaaaTCG\n>yet_another_fasta_sequence\nAAACCTGCcctGttug\n"
                                                        ">fasta_sequence_again\nAAAAAAACCTGCTAGctaatat\n>and_again_sequence\nCTGATCGTTTAGCAGCA\n")) as mock_file:
            keys, fastas = fasta_splitting_by_sequence(mock_file)
            self.assertEqual(keys, ['some_fasta_sequence', 'yet_another_fasta_sequence', 'fasta_sequence_again', 'and_again_sequence'])
            self.assertEqual([fastas[key] for key in keys], ['ACTGCTCCCTGGaaaTCG', 'AAACCTGCcctGttug', 'AAAAAAACCTGCTAGctaatat', 'CTGATCGTTTAGCAGCA'])

    def test_parsing_fasta_files_multiple_rn_line_terminators(self):
        """Tests whether several short fasta-like strings can be split into DNA sequences and chromosome names when they have Windows specific line terminators."""
        with patch("builtins.open", mock_open(read_data=">some_fasta_sequence\r\nGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG\r\nGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG" \
                     "\r\n>fasta_sequence_again\r\nAAAAAAACCTGCTAGctaatatGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG\r\n>and_again_sequence\r\n" \
                     "CTGATCGTTTAGCAGCGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGA\r\n")) as mock_file:
            keys, fastas = fasta_splitting_by_sequence(mock_file)
            self.assertEqual(keys, ['some_fasta_sequence', 'fasta_sequence_again', 'and_again_sequence'])
            self.assertEqual([fastas[key] for key in keys],['GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG',
                                        'AAAAAAACCTGCTAGctaatatGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG',
                                        'CTGATCGTTTAGCAGCGGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGA'])


if __name__ == '__main__':
    unittest.main()
