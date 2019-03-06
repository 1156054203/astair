import sys
import unittest

from astair import DNA_sequences_operations as so
    
class DNAOperationsTest(unittest.TestCase):
    """Tests whether the resulting DNA sequences match the expected reverse, complementary and reverse complementary output."""
    
    def test_reverse_DNA(self):
        """Tests whether the resulting DNA sequences match the expected reverse output."""
        self.assertEqual(so.reverse('AAAACTCTGCCGggAaaccttCG'), 'GCttccaaAggGCCGTCTCAAAA')

    def test_complementary_DNA(self):
        """Tests whether the resulting DNA sequences match the expected complementary output."""
        self.assertEqual(so.complementary('AAAACTCTGCCGggAaaccttCG'), 'TTTTGAGACGGCccTttggaaGC')
        
    def test_reverse_complementary_DNA(self):
        """Tests whether the resulting DNA sequences match the expected reverse complementary output."""
        self.assertEqual(so.reverse_complementary('AAAACTCTGCCGggAaaccttCG'), 'CGaaggttTccCGGCAGAGTTTT')


if __name__ == '__main__':
    unittest.main()
