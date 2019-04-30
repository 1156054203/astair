import pdb
import csv
import unittest
import subprocess
from os import path

from astair.separator import position_separator

current = path.abspath(path.dirname(__file__))


class SeparatorOutputTest(unittest.TestCase):
    """Tests whether the modification levels at known positions in a synthetic NCNN spike-in will be equal to the expected."""

    def test_separator_default_taps(self):
        """Looks for TAPS sequencing reads matching CpG, CHG and CHH contexts and calculate their modification levels."""
        position_separator(current + '/test_data/small_real_taps_synthetic.bam', 80, 'mCtoT', [51,111], ['OT','OB'], False, False, current + '/test_data/')
        data_generated = list()
        with open(current + '/test_data/small_real_taps_synthetic_mCtoT_positions.mods','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated, [('CHROM', 'START', 'END', 'CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'MOD', 'UNMOD'), ('5mC_five_prime', '50', '51', 'CHG', '84.836', '800', '143'), ('5mC_five_prime', '50', '51', 'CHH', '84.873', '2003', '357'), ('5mC_five_prime', '50', '51', 'CpG', '96.92', '1762', '56'), ('5mC_five_prime', '110', '111', 'CHG', '81.787', '952', '212'), ('5mC_five_prime', '110', '111', 'CHH', '80.476', '1624', '394'), ('5mC_five_prime', '110', '111', 'CpG', '88.038', '1965', '267')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_synthetic_mCtoT_positions.mods')
        subprocess.Popen(remove, shell=True)
        

if __name__ == '__main__':
    unittest.main()
