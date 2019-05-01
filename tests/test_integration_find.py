import pdb
import csv
import unittest
import subprocess
from os import path

from astair.finder import find_contexts

current = path.abspath(path.dirname(__file__))


class FindOutputTest(unittest.TestCase):
    """Tests whether the same cytosine positions will be called as modified for both TAPS and WGBS."""

    def test_find_default_all(self):
        """Looks for cytosine positions and compares them with the expected ones by context."""
        find_contexts(current + '/test_data/lambda_phage.fa', 'all', None, None, False, current + '/test_data/')
        data_generated = list()
        with open(current + '/test_data/lambda_phage_all.bed','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[0:20], [('CHROM', 'START', 'END', 'STRAND', 'SPECIFIC_CONTEXT', 'CONTEXT'), ('lambda', '2', '3', '-', 'CCC', 'CHH'), ('lambda', '3', '4', '+', 'CGG', 'CpG'), ('lambda', '4', '5', '-', 'CGC', 'CpG'), ('lambda', '5', '6', '-', 'CCG', 'CHG'), ('lambda', '6', '7', '+', 'CGA', 'CpG'), ('lambda', '7', '8', '-', 'CGC', 'CpG'), ('lambda', '9', '10', '+', 'CCT', 'CHH'), ('lambda', '10', '11', '+', 'CTC', 'CHH'), ('lambda', '12', '13', '+', 'CGC', 'CpG'), ('lambda', '13', '14', '-', 'CGA', 'CpG'), ('lambda', '14', '15', '+', 'CGG', 'CpG'), ('lambda', '15', '16', '-', 'CGC', 'CpG'), ('lambda', '16', '17', '-', 'CCG', 'CHG'), ('lambda', '17', '18', '-', 'CCC', 'CHH'), ('lambda', '22', '23', '+', 'CGC', 'CpG'), ('lambda', '23', '24', '-', 'CGA', 'CpG'), ('lambda', '24', '25', '+', 'CTA', 'CHH'), ('lambda', '32', '33', '-', 'CAT', 'CHH'), ('lambda', '41', '42', '+', 'CCG', 'CHG')])
        remove = 'rm {}'.format(current + '/test_data/lambda_phage_all.bed*')
        subprocess.Popen(remove, shell=True)
        
if __name__ == '__main__':
    unittest.main()
