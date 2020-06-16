import pdb
import csv
import unittest
import subprocess
from os import path

from astair.idbiaser import IDbias_plotting


current = path.abspath(path.dirname(__file__))


class IDbiasOutputTest(unittest.TestCase):
    """Tests whether the read information will be used to correctly
    calculate the indel rate per read length for both TAPS and WGBS."""

    def test_idbias_default_taps(self):
        """Tests whether the read information will be used to correctly
        calculate the indel rate per read length in a real TAPS sample."""
        IDbias_plotting(current + '/test_data/lambda_phage.fa', current + '/test_data/small_real_taps_lambda_mCtoT.bam', current + '/test_data/', 80, 'mCtoT', False, False, ['teal', 'deepskyblue', 'mediumblue', 'orange', 'gold', 'sienna'], 1, None, '.')
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_ID-bias.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row[0:5])))
        self.assertEqual(data_generated[0:5], [('#POS', 'R1_total_reads', 'R1_total_insertions', 'R1_total_deletions', 'R1_total_both'), ('.', '193', '0', '52', '0'), ('1', '193', '0', '0', '.'), ('2', '193', '0', '0', '.'), ('3', '193', '0', '0', '.')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT*ID-bias.stats')
        subprocess.Popen(remove, shell=True)
        
    def test_idbias_default_taps_SE(self):
        """Tests whether the read information will be used to correctly
        calculate the indel rate per read length in a real single-end TAPS sample."""
        IDbias_plotting(current + '/test_data/lambda_phage.fa', current + '/test_data/small_real_taps_lambda_mCtoT_SE.bam', current + '/test_data/', 75, 'mCtoT', True, False, ['teal', 'deepskyblue', 'mediumblue', 'orange', 'gold', 'sienna'], 1, None, 0)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_SE_ID-bias.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row[0:5])))
        self.assertEqual(data_generated[0:5], [('#POS', 'R1_total_reads', 'R1_total_insertions', 'R1_total_deletions', 'R1_total_both'), ('0', '10', '0', '0', '0'), ('1', '10', '0', '0', '0'), ('2', '10', '0', '0', '0'), ('3', '10', '0', '0', '0')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT*ID-bias.stats')
        subprocess.Popen(remove, shell=True)
        
    def test_idbias_default_wgbs(self):
        """Tests whether the read information will be used to correctly
        calculate the indel rate per read length in a real WGBS sample."""
        IDbias_plotting(current + '/test_data/lambda_phage.fa', current + '/test_data/small_real_wgbs_lambda_CtoT.bam', current + '/test_data/', 80, 'CtoT', False, False, ['teal', 'deepskyblue', 'mediumblue', 'orange', 'gold', 'sienna'], 1, None, 0)
        data_generated = list()
        with open(current + '/test_data/small_real_wgbs_lambda_CtoT_ID-bias.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)[0:5]))     
        self.assertEqual(data_generated[0:5], [('#POS', 'R1_total_reads', 'R1_total_insertions', 'R1_total_deletions', 'R1_total_both'), ('0', '138', '0', '6', '0'), ('1', '138', '0', '0', '0'), ('2', '138', '0', '0', '0'), ('3', '138', '0', '0', '0')])
        remove = 'rm {}'.format(current + '/test_data/small_real_wgbs_lambda_CtoT*ID-bias.stats')
        subprocess.Popen(remove, shell=True)
        

if __name__ == '__main__':
    unittest.main()
