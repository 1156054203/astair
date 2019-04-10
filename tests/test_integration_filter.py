import pdb
import csv
import unittest
import subprocess
from os import path

from astair.caller import cytosine_modification_finder
from astair.filter import removing_mod_err

current = path.abspath(path.dirname(__file__))


class CallOutputTest(unittest.TestCase):
    """Tests whether the same cytosine positions will be called as modified for both TAPS and WGBS."""

    def test_filter_default_taps(self):
        """Looks for TAPS sequencing reads with more than N CpH modified positions and filters them out."""
        removing_mod_err(current + '/test_data/small_real_taps_lambda_mCtoT.bam', 'mCtoT', 3, 1, current + '/test_data/')
        cytosine_modification_finder(current + '/test_data/small_real_taps_lambda_mCtoT_high_CpH_filtered.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'mCtoT', 0, 0, True, True, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_high_CpH_filtered_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS'), ('CpG', '', '98.068', '6225', '48'), ('', 'CGA', '99.502', '1210', '12'), ('', 'CGC', '90.789', '1730', '13'), ('', 'CGG', '97.826', '1847', '13'), ('', 'CGT', '100.0', '1438', '10'), ('CHG', '', '0.108', '6451', '45'), ('', 'CAG', '0.238', '2302', '16'), ('', 'CCG', '0.0', '1847', '13'), ('', 'CTG', '0.0', '2302', '16'), ('CHH', '', '2.779', '11503', '81'), ('', 'CTT', '0.0', '1349', '15'), ('', 'CAT', '0.0', '1802', '14'), ('', 'CCT', '0.0', '1182', '14'), ('', 'CTA', '0.0', '501', '1'), ('', 'CAA', '0.0', '1432', '7'), ('', 'CCA', '0.0', '1610', '8'), ('', 'CTC', '25.843', '1116', '10'), ('', 'CAC', '0.676', '1474', '5'), ('', 'CCC', '0.0', '1037', '7')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT_high*')
        subprocess.Popen(remove, shell=True)
        
    
    def test_filter_default_wgbs(self):
        """Looks for WGBS sequencing reads with more than N CpH modified positions and filters them out."""
        removing_mod_err(current + '/test_data/small_real_wgbs_lambda_CtoT.bam', 'CtoT', 3, 1, current + '/test_data/')
        cytosine_modification_finder(current + '/test_data/small_real_wgbs_lambda_CtoT_high_CpH_filtered.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'CtoT', 0, 0, True, True, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_real_wgbs_lambda_CtoT_high_CpH_filtered_CtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS'), ('CpG', '', '95.284', '6225', '48'), ('', 'CGA', '93.05', '1210', '12'), ('', 'CGC', '96.414', '1730', '13'), ('', 'CGG', '93.519', '1847', '13'), ('', 'CGT', '97.4', '1438', '10'), ('CHG', '', '0.733', '6451', '45'), ('', 'CAG', '0.351', '2302', '16'), ('', 'CCG', '0.0', '1847', '13'), ('', 'CTG', '1.728', '2302', '16'), ('CHH', '', '0.596', '11503', '81'), ('', 'CTT', '0.218', '1349', '15'), ('', 'CAT', '0.0', '1802', '14'), ('', 'CCT', '1.111', '1182', '14'), ('', 'CTA', '0.0', '501', '1'), ('', 'CAA', '0.495', '1432', '7'), ('', 'CCA', '0.0', '1610', '8'), ('', 'CTC', '1.938', '1116', '10'), ('', 'CAC', '0.0', '1474', '5'), ('', 'CCC', '1.19', '1037', '7')])
        remove = 'rm {}'.format(current + '/test_data/small_real_wgbs_lambda_CtoT_high*')
        subprocess.Popen(remove, shell=True)

if __name__ == '__main__':
    unittest.main()
