import pdb
import csv
import unittest
import subprocess
from os import path

from astair.phred import Phred_scores_plotting


current = path.abspath(path.dirname(__file__))


class PhredOutputTest(unittest.TestCase):
    """Tests whether the cytosine positions modified at CpG context will be used to correctly
    calculate the CpG modification rate per read length for both TAPS and WGBS."""

    def test_phred_default_taps(self):
        """Tests whether the cytosine positions modified at CpG context will be used to correctly
        calculate the CpG modification rate per read length in a real TAPS sample."""
        Phred_scores_plotting(current + '/test_data/small_real_taps_lambda_1.fq.gz', current + '/test_data/small_real_taps_lambda_2.fq.gz', 'means', current + '/test_data', 772, 0, ['skyblue', 'mediumaquamarine', 'khaki', 'lightcoral'], False)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_total_Phred.txt','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        print(data_generated)
        self.assertEqual(data_generated, [('____________________________________First in pair_____________________________________',), ('______________________________________________________________________________________',), ('mean ', 'adenines: 34.814', 'cytosines: 35.039', 'thymines: 35.082', 'guanines: 34.924'), ('median ', 'adenines: 35.5', 'cytosines: 35.579', 'thymines: 35.652', 'guanines: 35.556'), ('q25 ', 'adenines: 34.654', 'cytosines: 35.2', 'thymines: 35.273', 'guanines: 35.158'), ('q75 ', 'adenines: 35.789', 'cytosines: 36.0', 'thymines: 35.852', 'guanines: 35.765'), ('sd ', 'adenines: 1.718', 'cytosines: 1.808', 'thymines: 1.774', 'guanines: 1.932'), ('min ', 'adenines: 24.947', 'cytosines: 22.889', 'thymines: 23.633', 'guanines: 24.053'), ('max ', 'adenines: 36.0', 'cytosines: 36.0', 'thymines: 36.0', 'guanines: 36.0'), ('______________________________________________________________________________________',), ('____________________________________Second in pair____________________________________',), ('______________________________________________________________________________________',), ('mean ', 'adenines: 34.584', 'cytosines: 34.962', 'thymines: 34.992', 'guanines: 34.246'), ('median ', 'adenines: 35.31', 'cytosines: 35.636', 'thymines: 35.778', 'guanines: 34.667'), ('q25 ', 'adenines: 34.421', 'cytosines: 35.294', 'thymines: 35.267', 'guanines: 34.462'), ('q75 ', 'adenines: 35.789', 'cytosines: 35.778', 'thymines: 36.0', 'guanines: 35.667'), ('sd ', 'adenines: 1.956', 'cytosines: 2.046', 'thymines: 1.889', 'guanines: 2.52'), ('min ', 'adenines: 25.455', 'cytosines: 22.0', 'thymines: 25.056', 'guanines: 17.917'), ('max ', 'adenines: 36.0', 'cytosines: 36.0', 'thymines: 36.0', 'guanines: 36.0'), ('______________________________________________________________________________________',)])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_total_Phred.txt')
        subprocess.Popen(remove, shell=True)
        
    def test_phred_default_wgbs(self):
        """Tests whether the cytosine positions modified at CpG context will be used to correctly
        calculate the CpG modification rate per read length in a real WGBS sample."""
        Phred_scores_plotting(current + '/test_data/small_real_wgbs_lambda_1.fq.gz', current + '/test_data/small_real_wgbs_lambda_2.fq.gz', 'means', current + '/test_data', 552, 0, ['skyblue', 'mediumaquamarine', 'khaki', 'lightcoral'], False)
        data_generated = list()
        with open(current + '/test_data/small_real_wgbs_lambda_total_Phred.txt','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        print(data_generated)
        self.assertEqual(data_generated, [('____________________________________First in pair_____________________________________',), ('______________________________________________________________________________________',), ('mean ', 'adenines: 33.032', 'cytosines: 26.475', 'thymines: 33.415', 'guanines: 34.624'), ('median ', 'adenines: 33.833', 'cytosines: 26.819', 'thymines: 34.35', 'guanines: 35.158'), ('q25 ', 'adenines: 31.479', 'cytosines: 22.917', 'thymines: 32.095', 'guanines: 34.767'), ('q75 ', 'adenines: 35.5', 'cytosines: 30.0', 'thymines: 35.293', 'guanines: 35.158'), ('sd ', 'adenines: 2.743', 'cytosines: 4.843', 'thymines: 2.361', 'guanines: 1.638'), ('min ', 'adenines: 22.833', 'cytosines: 14.0', 'thymines: 26.325', 'guanines: 25.722'), ('max ', 'adenines: 36.0', 'cytosines: 36.0', 'thymines: 36.0', 'guanines: 36.0'), ('______________________________________________________________________________________',), ('____________________________________Second in pair____________________________________',), ('______________________________________________________________________________________',), ('mean ', 'adenines: 29.611', 'cytosines: 33.095', 'thymines: 32.397', 'guanines: 21.906'), ('median ', 'adenines: 30.789', 'cytosines: 33.983', 'thymines: 33.774', 'guanines: 21.225'), ('q25 ', 'adenines: 27.486', 'cytosines: 32.174', 'thymines: 30.481', 'guanines: 18.542'), ('q75 ', 'adenines: 32.615', 'cytosines: 35.052', 'thymines: 34.858', 'guanines: 24.333'), ('sd ', 'adenines: 4.162', 'cytosines: 2.844', 'thymines: 3.543', 'guanines: 4.779'), ('min ', 'adenines: 16.054', 'cytosines: 21.0', 'thymines: 16.75', 'guanines: 14.0'), ('max ', 'adenines: 35.19', 'cytosines: 35.852', 'thymines: 36.0', 'guanines: 34.0'), ('______________________________________________________________________________________',)])
        remove = 'rm {}'.format(current + '/test_data/small_real_wgbs_lambda_total_Phred.txt')
        subprocess.Popen(remove, shell=True)
        

if __name__ == '__main__':
    unittest.main()
