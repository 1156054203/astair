import pdb
import csv
import unittest
import subprocess
from os import path

from astair.caller import cytosine_modification_finder

current = path.abspath(path.dirname(__file__))


class CallOutputTest(unittest.TestCase):
    """Tests whether the same cytosine positions will be called as modified for both TAPS and WGBS."""

    def test_call_default_taps(self):
        """Looks for TAPS modified cytosine positions and compares their modification level with
        the expected one by context."""
        cytosine_modification_finder(current + '/test_data/small_lambda.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'mCtoT', 0, 0, True, True, True, False, True, True, 250, None, 1, current + '/test_data/')
        data_generated = list()
        with open(current + '/test_data/small_lambda_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS'), 
                                               ('CpG', '', '63.177', '6225', '46'), ('', 'CGA', '53.409', '1210', '12'), ('', 'CGC', '97.222', '1730', '11'),
                                               ('', 'CGG', '53.556', '1847', '13'), ('', 'CGT', '57.93', '1438', '10'), ('CHG', '', '0.129', '6451', '44'), 
                                               ('', 'CAG', '0.072', '2302', '15'), ('', 'CCG', '0.099', '1847', '13'), ('', 'CTG', '0.201', '2302', '16'), 
                                               ('CHH', '', '25.366', '11503', '81'), ('', 'CTT', '18.301', '1349', '15'), ('', 'CAT', '29.12', '1802', '14'), 
                                               ('', 'CCT', '29.799', '1182', '14'), ('', 'CTA', '0.0', '501', '1'), ('', 'CAA', '41.735', '1432', '7'),
                                               ('', 'CCA', '31.498', '1610', '8'), ('', 'CTC', '8.27', '1116', '10'), ('', 'CAC', '41.093', '1474', '5'), 
                                               ('', 'CCC', '3.939', '1037', '7')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)
        
        

    def test_call_default_wgbs(self):
        """Looks for WGBS modified cytosine positions and compares their modification level with
        the expected one by context."""
        cytosine_modification_finder(current + '/test_data/small_lambda.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'CtoT', 0, 0, True, True, True, False, True, True, 250, None, 1, current + '/test_data/')
        data_generated = list()
        with open(current + '/test_data/small_lambda_CtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS'),
                                               ('CpG', '', '36.823', '6225', '46'), ('', 'CGA', '46.591', '1210', '12'), ('', 'CGC', '2.778', '1730', '11'), 
                                               ('', 'CGG', '46.444', '1847', '13'), ('', 'CGT', '42.07', '1438', '10'), ('CHG', '', '99.871', '6451', '44'),
                                               ('', 'CAG', '99.928', '2302', '15'), ('', 'CCG', '99.901', '1847', '13'), ('', 'CTG', '99.799', '2302', '16'),
                                               ('CHH', '', '74.634', '11503', '81'), ('', 'CTT', '81.699', '1349', '15'), ('', 'CAT', '70.88', '1802', '14'),
                                               ('', 'CCT', '70.201', '1182', '14'), ('', 'CTA', '100.0', '501', '1'), ('', 'CAA', '58.265', '1432', '7'),
                                               ('', 'CCA', '68.502', '1610', '8'), ('', 'CTC', '91.73', '1116', '10'), ('', 'CAC', '58.907', '1474', '5'), 
                                               ('', 'CCC', '96.061', '1037', '7')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_CtoT_all.*')
        subprocess.Popen(remove, shell=True)

if __name__ == '__main__':
    unittest.main()
