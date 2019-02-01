import unittest
import xmlrunner

from safe_division import non_zero_division
from safe_division import non_zero_division_NA


class SafeDivisionFunctionsOutputTest(unittest.TestCase):
    """Tests whether the expected 0 or NA will be returned
    as a result of a division by zero operation."""
    def test_non_zero_division(self):
        self.assertEqual(non_zero_division(10,0), 0)
    def test_non_zero_division_NA(self):
        self.assertEqual(non_zero_division_NA(10,0), 'NA')

if __name__ == '__main__':
    unittest.main(testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
