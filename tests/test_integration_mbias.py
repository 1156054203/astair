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
        Mbias_plotting(current + '/test_data/small_real_taps_lambda_mCtoT.bam', current + '/test_data/', 80, 'mCtoT', False, False, ['teal', 'gray', 'maroon'], 1)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_Mbias.txt','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row[0:2])))
        self.assertEqual(data_generated, [('POSITION_(bp)', 'MOD_LVL_CpG_READ_1'), ('1', '100.0'), ('2', '100.0'), ('3', '71.429'), ('4', '52.941'), ('5', '100.0'), ('6', '66.667'), ('7', '50.0'), ('8', '6.667'), ('9', '100.0'), ('10', '100.0'), ('11', '100.0'), ('12', '100.0'), ('13', '100.0'), ('14', '100.0'), ('15', '100.0'), ('16', '100.0'), ('17', '100.0'), ('18', '50.0'), ('19', '100.0'), ('20', '100.0'), ('21', '100.0'), ('22', '100.0'), ('23', '100.0'), ('24', '100.0'), ('25', '100.0'), ('26', '75.0'), ('27', '70.0'), ('28', '97.059'), ('29', '100.0'), ('30', '100.0'), ('31', '100.0'), ('32', '85.714'), ('33', '88.889'), ('34', '100.0'), ('35', '100.0'), ('36', '75.0'), ('37', '88.889'), ('38', '97.059'), ('39', '80.0'), ('40', '75.0'), ('41', '42.105'), ('42', '100.0'), ('43', '86.667'), ('44', '12.5'), ('45', '80.0'), ('46', '90.909'), ('47', '50.0'), ('48', '85.714'), ('49', '60.0'), ('50', '53.846'), ('51', '50.0'), ('52', '80.0'), ('53', '100.0'), ('54', '16.129'), ('55', '50.0'), ('56', '80.0'), ('57', '80.0'), ('58', '88.889'), ('59', '94.737'), ('60', '60.0'), ('61', '60.0'), ('62', '92.308'), ('63', '42.857'), ('64', '70.0'), ('65', '100.0'), ('66', '20.0'), ('67', '100.0'), ('68', '55.556'), ('69', '91.304'), ('70', '83.333'), ('71', '75.0'), ('72', '90.0'), ('73', '100.0'), ('74', '72.727'), ('75', '90.0'), ('76', '80.0'), ('77', '80.0'), ('78', '66.667'), ('79', '75.0'), ('80', '60.0')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT*Mbias.txt')
        subprocess.Popen(remove, shell=True)
        
    def test_mbias_default_wgbs(self):
        """Tests whether the cytosine positions modified at CpG context will be used to correctly
        calculate the CpG modification rate per read length in a real WGBS sample."""
        Mbias_plotting(current + '/test_data/small_real_wgbs_lambda_CtoT.bam', current + '/test_data/', 80, 'CtoT', False, False, ['teal', 'gray', 'maroon'], 1)
        data_generated = list()
        with open(current + '/test_data/small_real_wgbs_lambda_CtoT_Mbias.txt','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)[0:2]))
        self.assertEqual(data_generated[-10:], [('71', '100.0'), ('72', '100.0'), ('73', '100.0'), ('74', '100.0'), ('75', '50.0'), ('76', '100.0'), ('77', '100.0'), ('78', '100.0'), ('79', '100.0'), ('80', '100.0')])
        remove = 'rm {}'.format(current + '/test_data/small_real_wgbs_lambda_CtoT*Mbias.txt')
        subprocess.Popen(remove, shell=True)
        

if __name__ == '__main__':
    unittest.main()
