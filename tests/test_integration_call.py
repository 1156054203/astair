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
        cytosine_modification_finder(current + '/test_data/small_lambda.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'mCtoT', 0, 0, True, True, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_lambda_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '63.177', '6225', '46', '1961', '1143'), ('*', 'CGA', '53.409', '1210', '12', '423', '369'), ('*', 'CGC', '97.222', '1730', '11', '595', '17'), ('*', 'CGG', '53.556', '1847', '13', '512', '444'), ('*', 'CGT', '57.93', '1438', '10', '431', '313'), ('CHG', '*', '0.129', '6451', '44', '5', '3885'), ('*', 'CAG', '0.072', '2302', '15', '1', '1388'), ('*', 'CCG', '0.099', '1847', '13', '1', '1005'), ('*', 'CTG', '0.201', '2302', '16', '3', '1492'), ('CHH', '*', '25.366', '11503', '81', '1663', '4893'), ('*', 'CTT', '18.301', '1349', '15', '209', '933'), ('*', 'CAT', '29.12', '1802', '14', '341', '830'), ('*', 'CCT', '29.799', '1182', '14', '326', '768'), ('*', 'CTA', '0.0', '501', '1', '0', '20'), ('*', 'CAA', '41.735', '1432', '7', '255', '356'), ('*', 'CCA', '31.498', '1610', '8', '246', '535'), ('*', 'CTC', '8.27', '1116', '10', '65', '721'), ('*', 'CAC', '41.093', '1474', '5', '203', '291'), ('*', 'CCC', '3.939', '1037', '7', '18', '439')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)
        
        
    def test_call_default_taps_SE(self):
        """Looks for TAPS modified cytosine positions and compares their modification level with
        the expected one by context."""
        cytosine_modification_finder(current + '/test_data/small_real_taps_lambda_mCtoT_SE.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'mCtoT', 0, 0, True, True, True, False, True, True, 250, None, 1, current + '/test_data/', False, True)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_SE_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:6],[('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '100.0', '6225', '103', '42', '0'), ('*', 'CGA', '100.0', '1210', '24', '12', '0'), ('*', 'CGC', '100.0', '1730', '21', '7', '0'), ('*', 'CGG', '100.0', '1847', '33', '8', '0'), ('*', 'CGT', '100.0', '1438', '25', '15', '0')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT_SE_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)
        

    def test_call_default_wgbs(self):
        """Looks for WGBS modified cytosine positions and compares their modification level with
        the expected one by context."""
        cytosine_modification_finder(current + '/test_data/small_lambda.bam', current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'CtoT', 0, 0, True, True, True, False, True, True, 250, None, 1, current + '/test_data/', False, False)
        data_generated = list()
        with open(current + '/test_data/small_lambda_CtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '36.823', '6225', '46', '1143', '1961'), ('*', 'CGA', '46.591', '1210', '12', '369', '423'), ('*', 'CGC', '2.778', '1730', '11', '17', '595'), ('*', 'CGG', '46.444', '1847', '13', '444', '512'), ('*', 'CGT', '42.07', '1438', '10', '313', '431'), ('CHG', '*', '99.871', '6451', '44', '3885', '5'), ('*', 'CAG', '99.928', '2302', '15', '1388', '1'), ('*', 'CCG', '99.901', '1847', '13', '1005', '1'), ('*', 'CTG', '99.799', '2302', '16', '1492', '3'), ('CHH', '*', '74.634', '11503', '81', '4893', '1663'), ('*', 'CTT', '81.699', '1349', '15', '933', '209'), ('*', 'CAT', '70.88', '1802', '14', '830', '341'), ('*', 'CCT', '70.201', '1182', '14', '768', '326'), ('*', 'CTA', '100.0', '501', '1', '20', '0'), ('*', 'CAA', '58.265', '1432', '7', '356', '255'), ('*', 'CCA', '68.502', '1610', '8', '535', '246'), ('*', 'CTC', '91.73', '1116', '10', '721', '65'), ('*', 'CAC', '58.907', '1474', '5', '291', '203'), ('*', 'CCC', '96.061', '1037', '7', '439', '18')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_CtoT_all.*')
        subprocess.Popen(remove, shell=True)

if __name__ == '__main__':
    unittest.main()
