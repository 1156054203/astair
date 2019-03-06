import unittest

from astair import safe_division


class SafeDivisionFunctionsOutputTest(unittest.TestCase):
    """Tests whether the expected 0 or NA will be returned
    as a result of a division by zero operation."""
    def test_non_zero_division(self):
        """Tests whether zero will be returned from division
        by zero operation."""
        self.assertEqual(safe_division.non_zero_division(10,0), 0)
    def test_non_zero_division_NA(self):
        """Tests whether NA will be returned from division
        by zero operation."""
        self.assertEqual(safe_division.non_zero_division_NA(10,0), 'NA')

if __name__ == '__main__':
    unittest.main()
