import sys
import unittest

from version_testing import version_testing_builtin
from astair_phred_values_v3 import Phred_score_main_body

class PhredScoreStatisticsTest(unittest.TestCase):
    """Tests whether the resulting DNA sequences match the expected reverse, complementary and reverse complementary output."""
    
    def test_Phred_means(self):
        """Tests whether the resulting DNA sequences match the expected reverse output."""
        read_values = Phred_score_main_body(["@lambda-2/1", "TACTTTCTGTACCAGAAAAATGACGCCTGACTCTGGCTATCTGCTCGTTAAATCTGGCCGTATTCAGATTCAAAT", "+", "AAAAAEEEEEEEAEEE6EEEEEEEEEAEEEEEEEEEEEEAAEEEEEEEEEEEEEEE<EEEEEEEEE0AEEEEEEE"], 'means', 5)
        self.assertEqual(read_values, [(35.333333333333336, 35.111111111111114, 34.85, 33.69230769230769)])

   
if __name__ == '__main__':
    unittest.main()
