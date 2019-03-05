import sys
import unittest

from DNA_sequences_operations import reverse
from DNA_sequences_operations import complementary
from DNA_sequences_operations import reverse_complementary
    
class DNAOperationsTest(unittest.TestCase):
    """Tests whether the resulting DNA sequences match the expected reverse, complementary and reverse complementary output."""
    
    def test_reverse_DNA(self):
        """Tests whether the resulting DNA sequences match the expected reverse output."""
        self.assertEqual(reverse('AAAACTCTGCCGggAaaccttCG'), 'GCttccaaAggGCCGTCTCAAAA')

    def test_complementary_DNA(self):
        """Tests whether the resulting DNA sequences match the expected complementary output."""
        self.assertEqual(complementary('AAAACTCTGCCGggAaaccttCG'), 'TTTTGAGACGGCccTttggaaGC')
        
    def test_reverse_complementary_DNA(self):
        """Tests whether the resulting DNA sequences match the expected reverse complementary output."""
        self.assertEqual(reverse_complementary('AAAACTCTGCCGggAaaccttCG'), 'CGaaggttTccCGGCAGAGTTTT')


if __name__ == '__main__':
    unittest.main()
