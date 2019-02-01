import unittest

from astair_mod_caller_v3 import statistics_calculator




class StatisticsCalculationTest(unittest.TestCase):
    """"Tests whether statistics_calculator function will count the modified and unmodified reads per cytosine context
    when it is given an ordered list input containing chromosome, start, end, modification ratio, number of
    modified reads, number of unmodified reads, modified base, reference base, specific context, general context,
     SNV or not."""
    def test_modified_position_correct(self):
        mean_mod = {'CHH':0, 'CAT':0}
        mean_unmod = {'CHH':0, 'CAT':0}
        data_mod = ['some_reference_genome', 1, 2, 0.8, 8, 2, 'T', 'C', 'CAT', 'CHH', 'No']
        statistics_calculator(mean_mod, mean_unmod, data_mod, None)
        self.assertEqual(mean_mod, {'CHH': 8, 'CAT': 8})
        self.assertEqual(mean_unmod, {'CHH': 2, 'CAT': 2})
    def test_modified_position_incorrect(self):
        mean_mod = {'CHH':0, 'CAT':0}
        mean_unmod = {'CHH':0, 'CAT':0}
        data_mod = ['some_reference_genome', 1, 2, 0.8, 8, 2, 'T', 'C', 'CAT', 'CHH', 'homozyguous']
        statistics_calculator(mean_mod, mean_unmod, data_mod, None)
        self.assertNotEqual(mean_mod, {'CHH': 8, 'CAT': 8})
        self.assertNotEqual(mean_unmod, {'CHH': 2, 'CAT': 2})



if __name__ == '__main__':
    unittest.main()
