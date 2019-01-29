import unittest

from astair_mod_caller_v3 import statistics_calculator




class StatisticsCalculationTest(unittest.TestCase):
    def test_modified_position_correct(self):
        mean_mod = {'CHH':0, 'CAT':0}
        mean_nmod = {'CHH':0, 'CAT':0}
        data_mod = ['some_reference_genome', 1, 2, 0.8, 8, 2, 'T', 'C', 'CAT', 'CHH', 'No']
        statistics_calculator(mean_mod, mean_nmod, data_mod, None)
        self.assertEqual(mean_mod, {'CHH': 8, 'CAT': 8})
        self.assertEqual(mean_nmod, {'CHH': 2, 'CAT': 2})
    def test_modified_position_incorrect(self):
        mean_mod = {'CHH':0, 'CAT':0}
        mean_nmod = {'CHH':0, 'CAT':0}
        data_mod = ['some_reference_genome', 1, 2, 0.8, 8, 2, 'T', 'C', 'CAT', 'CHH', 'homozyguous']
        statistics_calculator(mean_mod, mean_nmod, data_mod, None)
        self.assertEqual(mean_mod, {'CHH': 8, 'CAT': 8})
        self.assertEqual(mean_nmod, {'CHH': 2, 'CAT': 2})



if __name__ == '__main__':
    unittest.main()