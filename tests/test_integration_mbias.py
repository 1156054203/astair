import pdb
import csv
import unittest
import subprocess
from os import path

from astair.mbias import Mbias_plotting


current = path.abspath(path.dirname(__file__))


class MbiasOutputTest(unittest.TestCase):
    """Tests whether the cytosine positions modified at CpG context will be used to correctly
    calculate the CpG modification rate per read length for both TAPS and WGBS."""

    def test_mbias_default_taps(self):
        """Tests whether the cytosine positions modified at CpG context will be used to correctly
        calculate the CpG modification rate per read length in a real TAPS sample."""
        Mbias_plotting(current + '/test_data/lambda_phage.fa', current + '/test_data/small_real_taps_lambda_mCtoT.bam', current + '/test_data/', 80, 'mCtoT', False, False, ['teal', 'gray', 'maroon'], 1, None, '.')
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_Mbias.txt','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row[0:2])))
        self.assertEqual(data_generated,[('#POSITION_(bp)', 'MOD_LVL_CpG_READ_1'), ('1', '100.0'), ('2', '100.0'), ('3', '33.333'), ('4', '46.667'), ('5', '100.0'), ('6', '66.667'), ('7', '33.333'), ('8', '100.0'), ('9', '100.0'), ('10', '100.0'), ('11', '100.0'), ('12', '100.0'), ('13', '100.0'), ('14', '100.0'), ('15', '100.0'), ('16', '66.667'), ('17', '100.0'), ('18', '100.0'), ('19', '100.0'), ('20', '100.0'), ('21', '66.667'), ('22', '100.0'), ('23', '94.595'), ('24', '100.0'), ('25', '100.0'), ('26', '75.0'), ('27', '100.0'), ('28', '100.0'), ('29', '66.667'), ('30', '100.0'), ('31', '100.0'), ('32', '100.0'), ('33', '100.0'), ('34', '100.0'), ('35', '100.0'), ('36', '100.0'), ('37', '100.0'), ('38', '100.0'), ('39', '75.0'), ('40', '100.0'), ('41', '100.0'), ('42', '100.0'), ('43', '97.561'), ('44', '100.0'), ('45', '100.0'), ('46', '100.0'), ('47', '100.0'), ('48', '100.0'), ('49', '75.0'), ('50', '100.0'), ('51', '100.0'), ('52', '88.889'), ('53', '100.0'), ('54', '100.0'), ('55', '100.0'), ('56', '100.0'), ('57', '100.0'), ('58', '100.0'), ('59', '100.0'), ('60', '100.0'), ('61', '100.0'), ('62', '100.0'), ('63', '100.0'), ('64', '87.5'), ('65', '100.0'), ('66', '.'), ('67', '100.0'), ('68', '100.0'), ('69', '100.0'), ('70', '100.0'), ('71', '80.0'), ('72', '100.0'), ('73', '100.0'), ('74', '100.0'), ('75', '100.0'), ('76', '66.667'), ('77', '100.0'), ('78', '100.0'), ('79', '57.143'), ('80', '33.333')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT*Mbias.txt')
        subprocess.Popen(remove, shell=True)
        
    def test_mbias_default_taps_SE(self):
        """Tests whether the cytosine positions modified at CpG context will be used to correctly
        calculate the CpG modification rate per read length in a real single-end TAPS sample."""
        Mbias_plotting(current + '/test_data/lambda_phage.fa', current + '/test_data/small_real_taps_lambda_mCtoT_SE.bam', current + '/test_data/', 75, 'mCtoT', True, False, ['teal', 'gray', 'maroon'], 1, None, 0)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_SE_Mbias.txt','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row[0:2])))
        self.assertEqual(data_generated[0:5], [('#POSITION_(bp)', 'MOD_LVL_CpG_READ_OT'), ('1', '100.0'), ('2', '100.0'), ('3', '100.0'), ('4', '100.0')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT*Mbias.txt')
        subprocess.Popen(remove, shell=True)
        
    def test_mbias_default_wgbs(self):
        """Tests whether the cytosine positions modified at CpG context will be used to correctly
        calculate the CpG modification rate per read length in a real WGBS sample."""
        Mbias_plotting(current + '/test_data/lambda_phage.fa', current + '/test_data/small_real_wgbs_lambda_CtoT.bam', current + '/test_data/', 80, 'CtoT', False, False, ['teal', 'gray', 'maroon'], 1, None, 0)
        data_generated = list()
        with open(current + '/test_data/small_real_wgbs_lambda_CtoT_Mbias.txt','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)[0:2]))
        self.assertEqual(data_generated[-10:-1], [('71', '92.857'), ('72', '100.0'), ('73', '100.0'), ('74', '100.0'), ('75', '84.615'), ('76', '71.429'), ('77', '100.0'), ('78', '37.5'), ('79', '52.941')])
        remove = 'rm {}'.format(current + '/test_data/small_real_wgbs_lambda_CtoT*Mbias.txt')
        subprocess.Popen(remove, shell=True)
        

if __name__ == '__main__':
    unittest.main()
