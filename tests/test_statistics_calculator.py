import unittest
from collections import defaultdict

from astair import astair_mod_caller_v3 as mod_caller


class StatisticsCalculationTest(unittest.TestCase):
    """"Tests whether mod_caller.statistics_calculator function will count the modified and unmodified reads per cytosine context
    when it is given an ordered list input containing chromosome, start, end, modification ratio, number of
    modified reads, number of unmodified reads, modified base, reference base, specific context, general context,
     SNV or not."""

    def test_modified_position_correct(self):
        """Tests whether a CHH context will be discovered and its modified and unmodified bases counted."""
        mean_mod = {'CHH':0, 'CAT':0}
        mean_unmod = {'CHH':0, 'CAT':0}
        context_sample_counts = defaultdict(int)
        data_mod = ['some_reference_genome', 1, 2, 0.8, 8, 2, 'T', 'C', 'CAT', 'CHH', 'No']
        mod_caller.statistics_calculator(mean_mod, mean_unmod, data_mod, None, context_sample_counts)
        self.assertEqual(mean_mod, {'CHH': 8, 'CAT': 8})
        self.assertEqual(mean_unmod, {'CHH': 2, 'CAT': 2})

    def test_modified_position_incorrect(self):
        """Tests whether a SNV in CHH context will not be counted as a modification."""
        mean_mod = {'CHH':0, 'CAT':0}
        mean_unmod = {'CHH':0, 'CAT':0}
        context_sample_counts = defaultdict(int)
        data_mod = ['some_reference_genome', 1, 2, 0.8, 8, 2, 'T', 'C', 'CAT', 'CHH', 'homozyguous']
        mod_caller.statistics_calculator(mean_mod, mean_unmod, data_mod, None, context_sample_counts)
        self.assertEqual(mean_mod, {'CHH': 0, 'CAT': 0})
        self.assertEqual(mean_unmod, {'CHH': 0, 'CAT': 0})



if __name__ == '__main__':
    unittest.main()
