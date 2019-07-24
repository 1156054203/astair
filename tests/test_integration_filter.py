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
        removing_mod_err(current + '/test_data/lambda_phage.fa', current + '/test_data/small_real_taps_lambda_mCtoT.bam', 'mCtoT', 3, None, 1, False, current + '/test_data/')
        cytosine_modification_finder(current + '/test_data/small_real_taps_lambda_mCtoT_high_CpH_filtered.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'directional', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_high_CpH_filtered_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[0:6], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '87.134', '6225', '48', '1578', '233'), ('*', 'CGA', '96.852', '1210', '12', '523', '17'), ('*', 'CGC', '64.479', '1730', '13', '334', '184'), ('*', 'CGG', '95.473', '1847', '13', '464', '22'), ('*', 'CGT', '96.255', '1438', '10', '257', '10')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT_high*')
        subprocess.Popen(remove, shell=True)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        

    def test_filter_default_taps_gzip(self):
        """Looks for TAPS sequencing reads with more than N CpH modified positions and filters them out with GZIP compressed fasta reference."""
        removing_mod_err(current + '/test_data/lambda_phage_.fa.gz', current + '/test_data/small_real_taps_lambda_mCtoT.bam', 'mCtoT', 3, None, 1, False, current + '/test_data/')
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        cytosine_modification_finder(current + '/test_data/small_real_taps_lambda_mCtoT_high_CpH_filtered.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'directional', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_high_CpH_filtered_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[0:6],[('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '87.134', '6225', '48', '1578', '233'), ('*', 'CGA', '96.852', '1210', '12', '523', '17'), ('*', 'CGC', '64.479', '1730', '13', '334', '184'), ('*', 'CGG', '95.473', '1847', '13', '464', '22'), ('*', 'CGT', '96.255', '1438', '10', '257', '10')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT_high*')
        subprocess.Popen(remove, shell=True)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        
    
    def test_filter_default_wgbs(self):
        """Looks for WGBS sequencing reads with more than N CpH modified positions and filters them out."""
        removing_mod_err(current + '/test_data/lambda_phage.fa', current + '/test_data/small_real_wgbs_lambda_CtoT.bam', 'CtoT', 3, None, 1, False,  current + '/test_data/')
        cytosine_modification_finder(current + '/test_data/small_real_wgbs_lambda_CtoT_high_CpH_filtered.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'directional', 'CtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_real_wgbs_lambda_CtoT_high_CpH_filtered_CtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[0:6], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '94.537', '6225', '48', '1419', '82'), ('*', 'CGA', '91.039', '1210', '12', '254', '25'), ('*', 'CGC', '95.349', '1730', '13', '246', '12'), ('*', 'CGG', '93.468', '1847', '13', '415', '29'), ('*', 'CGT', '96.923', '1438', '10', '504', '16')])
        remove = 'rm {}'.format(current + '/test_data/small_real_wgbs_lambda_CtoT_high*')
        subprocess.Popen(remove, shell=True)
        
        
    def test_filter_default_wgbs_gzip(self):
        """Looks for WGBS sequencing reads with more than N CpH modified positions and filters them out with GZIP compressed fasta reference."""
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        removing_mod_err(current + '/test_data/lambda_phage_.fa.gz', current + '/test_data/small_real_wgbs_lambda_CtoT.bam', 'CtoT', 3, None, 1, False,  current + '/test_data/')
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        cytosine_modification_finder(current + '/test_data/small_real_wgbs_lambda_CtoT_high_CpH_filtered.bam', current + '/test_data/lambda_phage_.fa.gz', 'all', False, False, 13, None, 'directional', 'CtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_real_wgbs_lambda_CtoT_high_CpH_filtered_CtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[0:6], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '94.537', '6225', '48', '1419', '82'), ('*', 'CGA', '91.039', '1210', '12', '254', '25'), ('*', 'CGC', '95.349', '1730', '13', '246', '12'), ('*', 'CGG', '93.468', '1847', '13', '415', '29'), ('*', 'CGT', '96.923', '1438', '10', '504', '16')])
        remove = 'rm {}'.format(current + '/test_data/small_real_wgbs_lambda_CtoT_high*')
        subprocess.Popen(remove, shell=True)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)

if __name__ == '__main__':
    unittest.main()
