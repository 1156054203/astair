import unittest

from astair_mod_caller_v3 import sequence_context_set_creation
from reference_context_search_triad import ahocorasick_search



class SequenceSearchOutputTest(unittest.TestCase):
    """Tests whether a particular expected cytosine context would be discovered on the top strand.
    The function sequence_context_set_creation expects desired contexts (CHH, CHG, CpG, all,
    user-defined) as input strings and returns a dictionary containing all possible context
    realisations as a list with key values the desired contexts. Then ahocorasick_search
    uses the  dictionary of contexts, the key values, the fasta DNA string, chromosome
    name and fills in a data_context dictionary. The data_context dictionary contains as
    key values tuples of (chromomes, start, end) of the desired contexts discovered in
    the provided fasta file, and as dictionary items the general and specific context
    discovered at that position. The test battery should have two failures."""

    def test_ahocorasick_search_CHH_top(self):
        """Tests whether the expected CHH positions will be discovered on the top DNA strand."""
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CHH', None)
        ahocorasick_search('CHH', contexts, 'AACTTCATCACT', 'test_string', None, data_context)
        self.assertEqual(data_context, {('test_string', 2, 3): ('CTT', 'CHH', 'T', 'C'),
                                        ('test_string', 5, 6): ('CAT', 'CHH', 'T', 'C'),
                                        ('test_string', 8, 9): ('CAC', 'CHH', 'T', 'C')})

    def test_ahocorasick_search_CHH_bottom(self):
        """Tests whether the expected CHH positions will be discovered on the bottom DNA strand."""
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CHH', None)
        ahocorasick_search('CHHb', contexts, 'AAGGCTTTGccc', 'test_string', None, data_context)
        self.assertEqual(data_context, {('test_string', 2, 3): ('CTT', 'CHH', 'A', 'G'),
                                        ('test_string', 3, 4): ('CCT', 'CHH', 'A', 'G'),
                                        ('test_string', 8, 9): ('CAA', 'CHH', 'A', 'G')})

    def test_ahocorasick_search_CHG_top(self):
        """Tests whether the expected CHG positions will be discovered on the top DNA strand."""
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CHG', None)
        ahocorasick_search('CHG', contexts, 'AACTTCAGCACT', 'test_string', None, data_context)
        self.assertEqual(data_context, {('test_string', 5, 6): ('CAG', 'CHG', 'T', 'C')})

    def test_ahocorasick_search_CHG_bottom(self):
        """Tests whether the expected CHG positions will be discovered on the bottom DNA strand."""
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CHG', None)
        ahocorasick_search('CHGb', contexts, 'AACTTCAGCACT', 'test_string', None, data_context)
        self.assertEqual(data_context, {('test_string', 7, 8): ('CTG', 'CHG', 'A', 'G')})

    def test_ahocorasick_search_CpG(self):
        """Tests whether the expected CpG positions will be discovered."""
        data_context = {}
        contexts, all_keys = sequence_context_set_creation('CpG', None)
        ahocorasick_search('CG', contexts, 'AAGCGTTTGccc', 'test_string', None, data_context)
        self.assertEqual(data_context, {('test_string', 3, 4): ('CG', 'CpG', 'T', 'C'),
                                        ('test_string', 4, 5): ('CG', 'CpG', 'A', 'G')})

if __name__ == '__main__':
    unittest.main()
