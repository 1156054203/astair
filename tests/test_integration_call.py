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
        cytosine_modification_finder(current + '/test_data/small_lambda.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_lambda_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '63.338', '6225', '46', '2023', '1171'), ('*', 'CGA', '53.383', '1210', '12', '434', '379'), ('*', 'CGC', '96.688', '1730', '11', '613', '21'), ('*', 'CGG', '54.018', '1847', '13', '531', '452'), ('*', 'CGT', '58.246', '1438', '10', '445', '319'), ('CHG', '*', '0.125', '6451', '44', '5', '3998'), ('*', 'CAG', '0.07', '2302', '15', '1', '1433'), ('*', 'CCG', '0.097', '1847', '13', '1', '1034'), ('*', 'CTG', '0.196', '2302', '16', '3', '1531'), ('CHH', '*', '25.23', '11503', '81', '1701', '5041'), ('*', 'CTT', '18.36', '1349', '15', '215', '956'), ('*', 'CAT', '29.035', '1802', '14', '349', '853'), ('*', 'CCT', '29.564', '1182', '14', '332', '791'), ('*', 'CTA', '0.0', '501', '1', '0', '23'), ('*', 'CAA', '41.746', '1432', '7', '263', '367'), ('*', 'CCA', '31.18', '1610', '8', '251', '554'), ('*', 'CTC', '8.447', '1116', '10', '68', '737'), ('*', 'CAC', '40.594', '1474', '5', '205', '300'), ('*', 'CCC', '3.766', '1037', '7', '18', '460')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)
        
    
    def test_call_default_taps_gzip(self):
        """Looks for TAPS modified cytosine positions and compares their modification level with
        the expected one by context with GZIP compressed fasta reference."""
        cytosine_modification_finder(current + '/test_data/small_lambda.bam', current + '/test_data/lambda_phage_.fa.gz', 'all', False, False, 13, None, 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_lambda_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '63.338', '6225', '46', '2023', '1171'), ('*', 'CGA', '53.383', '1210', '12', '434', '379'), ('*', 'CGC', '96.688', '1730', '11', '613', '21'), ('*', 'CGG', '54.018', '1847', '13', '531', '452'), ('*', 'CGT', '58.246', '1438', '10', '445', '319'), ('CHG', '*', '0.125', '6451', '44', '5', '3998'), ('*', 'CAG', '0.07', '2302', '15', '1', '1433'), ('*', 'CCG', '0.097', '1847', '13', '1', '1034'), ('*', 'CTG', '0.196', '2302', '16', '3', '1531'), ('CHH', '*', '25.23', '11503', '81', '1701', '5041'), ('*', 'CTT', '18.36', '1349', '15', '215', '956'), ('*', 'CAT', '29.035', '1802', '14', '349', '853'), ('*', 'CCT', '29.564', '1182', '14', '332', '791'), ('*', 'CTA', '0.0', '501', '1', '0', '23'), ('*', 'CAA', '41.746', '1432', '7', '263', '367'), ('*', 'CCA', '31.18', '1610', '8', '251', '554'), ('*', 'CTC', '8.447', '1116', '10', '68', '737'), ('*', 'CAC', '40.594', '1474', '5', '205', '300'), ('*', 'CCC', '3.766', '1037', '7', '18', '460')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        
        
    def test_call_default_taps_SE(self):
        """Looks for TAPS modified cytosine positions and compares their modification level with
        the expected one by context."""
        cytosine_modification_finder(current + '/test_data/small_real_taps_lambda_mCtoT_SE.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, True)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_SE_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:6],[('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '98.039', '6225', '103', '50', '1'), ('*', 'CGA', '100.0', '1210', '24', '12', '0'), ('*', 'CGC', '100.0', '1730', '21', '7', '0'), ('*', 'CGG', '94.118', '1847', '33', '16', '1'), ('*', 'CGT', '100.0', '1438', '25', '15', '0')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT_SE_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)



    def test_call_default_wgbs(self):
        """Looks for WGBS modified cytosine positions and compares their modification level with
        the expected one by context."""
        cytosine_modification_finder(current + '/test_data/small_lambda.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'CtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_lambda_CtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '36.662', '6225', '46', '1171', '2023'), ('*', 'CGA', '46.617', '1210', '12', '379', '434'), ('*', 'CGC', '3.312', '1730', '11', '21', '613'), ('*', 'CGG', '45.982', '1847', '13', '452', '531'), ('*', 'CGT', '41.754', '1438', '10', '319', '445'), ('CHG', '*', '99.875', '6451', '44', '3998', '5'), ('*', 'CAG', '99.93', '2302', '15', '1433', '1'), ('*', 'CCG', '99.903', '1847', '13', '1034', '1'), ('*', 'CTG', '99.804', '2302', '16', '1531', '3'), ('CHH', '*', '74.77', '11503', '81', '5041', '1701'), ('*', 'CTT', '81.64', '1349', '15', '956', '215'), ('*', 'CAT', '70.965', '1802', '14', '853', '349'), ('*', 'CCT', '70.436', '1182', '14', '791', '332'), ('*', 'CTA', '100.0', '501', '1', '23', '0'), ('*', 'CAA', '58.254', '1432', '7', '367', '263'), ('*', 'CCA', '68.82', '1610', '8', '554', '251'), ('*', 'CTC', '91.553', '1116', '10', '737', '68'), ('*', 'CAC', '59.406', '1474', '5', '300', '205'), ('*', 'CCC', '96.234', '1037', '7', '460', '18')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_CtoT_all.*')
        subprocess.Popen(remove, shell=True)
        
   
        
if __name__ == '__main__':
    unittest.main()
