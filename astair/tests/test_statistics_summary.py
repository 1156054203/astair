import unittest

from statistics_summary import general_statistics_summary


class StatisticsSummaryOutputTest(unittest.TestCase):
    """Tests whether the expected mean, median, sd, q25, q75, min
    and max value of a list of numerical values will be returned."""

    def test_statistics_summary(self):
        """Tests whether a correct statistics summary will be returned for a given numerical input."""
        numeric_data = list((1,2,3,40,215,100,5,4,12,16,76,81,9))
        statistics_output = list(general_statistics_summary(numeric_data))
        self.assertEqual(statistics_output, [43.385, 12.0, 62.068, 4.0, 76.0, 1, 215])

if __name__ == '__main__':
    unittest.main()
