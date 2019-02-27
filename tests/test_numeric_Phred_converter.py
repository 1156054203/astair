import unittest

from astair_phred_values_v3 import numeric_Phred_score


class NumericPhredOutputTest(unittest.TestCase):
    """Tests whether a correct translation of ASCII characters encoding Phred scores to numerical scores will be achieved."""

    def test_numeric_Phred_score(self):
        """Tests the translation of Phred scores ASCII string to numerical Phred scores."""
        test_score_string = '!\#$%&"()*+`-./0123456789:;<=>?@ABCDEFGHI'
        test_score_string2 = "!\#$%&'()*+`-./0123456789:;<=>?@ABCDEFGHI"
        self.assertEqual(numeric_Phred_score(test_score_string),
                         [0, 2, 3, 4, 5, 1, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                                                                  27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40])
        self.assertEqual(numeric_Phred_score(test_score_string2),
                         [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                                                                  27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40])



if __name__ == '__main__':
    unittest.main()
