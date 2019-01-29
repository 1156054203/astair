import unittest

from astair_mod_caller_v3 import sequence_context_set_creation
from astair_mod_caller_v3 import ahocorasick_search



class SequenceSearchOutputTest(unittest.TestCase):
    def test_ahocorasick_search_CHH_top_correct(self):
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CHH', None)
        ahocorasick_search('CHH', contexts, 'AACTTCATCACT', data_context, 'test_string', None)
        self.assertEqual(data_context, {('test_string', 2, 3): ('CTT', 'CHH', 'T', 'C'),
                                        ('test_string', 5, 6): ('CAT', 'CHH', 'T', 'C'),
                                        ('test_string', 8, 9): ('CAC', 'CHH', 'T', 'C')})
    def test_ahocorasick_search_CHH_top_incorrect(self):
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CHH', None)
        ahocorasick_search('CHH', contexts, 'AACTTCATCACT', data_context, 'test_string', None)
        self.assertEqual(data_context, {('test_string', 2, 3): ('CAA', 'CHH', 'T', 'C'),
                                        ('test_string', 5, 6): ('CTT', 'CHH', 'T', 'C'),
                                        ('test_string', 8, 9): ('CTA', 'CHH', 'T', 'C')})
    def test_ahocorasick_search_CHH_bottom_correct(self):
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CHH', None)
        ahocorasick_search('CHHb', contexts, 'AAGGCTTTGccc', data_context, 'test_string', None)
        self.assertEqual(data_context, {('test_string', 2, 3): ('CTT', 'CHH', 'A', 'G'),
                                        ('test_string', 3, 4): ('CCT', 'CHH', 'A', 'G'),
                                        ('test_string', 8, 9): ('CAA', 'CHH', 'A', 'G')})
    def test_ahocorasick_search_CHH_bottom_incorrect(self):
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CHH', None)
        ahocorasick_search('CHHb', contexts, 'AAGGCTTTGccc', data_context, 'test_string', None)
        self.assertEqual(data_context, {})
    def test_ahocorasick_search_CHG_top_correct(self):
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CHG', None)
        ahocorasick_search('CHG', contexts, 'AACTTCAGCACT', data_context, 'test_string', None)
        self.assertEqual(data_context, {('test_string', 5, 6): ('CAG', 'CHG', 'T', 'C')})
    def test_ahocorasick_search_CHG_bottom_correct(self):
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CHG', None)
        ahocorasick_search('CHGb', contexts, 'AACTTCAGCACT', data_context, 'test_string', None)
        self.assertEqual(data_context, {('test_string', 7, 8): ('CTG', 'CHG', 'A', 'G')})
    def test_ahocorasick_search_CpG_correct(self):
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CpG', None)
        ahocorasick_search('CG', contexts, 'AAGCGTTTGccc', data_context, 'test_string', None)
        self.assertEqual(data_context, {('test_string', 3, 4): ('CG', 'CpG', 'T', 'C'),
                                        ('test_string', 4, 5): ('CG', 'CpG', 'A', 'G')})

if __name__ == '__main__':
    unittest.main()
