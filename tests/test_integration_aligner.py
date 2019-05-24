import pdb
import csv
import sys
import unittest
import subprocess
from os import path
from tests.version_testing import version_testing_builtin

if sys.version[0] == '2':
    from mock import patch, mock_open
elif sys.version[0] == '3':
    from unittest.mock import patch, mock_open
else:
    raise Exception("This is not the python we're looking for (version {})".format(sys.version[0])) 


from astair.aligner import run_alignment

current = path.abspath(path.dirname(__file__))


class AlignFastaOutputTest(unittest.TestCase):
    """Tests whether the aligner function internals function as expected."""
    
    
    def test_mock_taps_pair_end(self):
        """Tests whether the aligner function will run with pair-end TAPS reads."""
        samtools = mock_open(read_data = current + '/test_data/samtools').return_value
        bwa = mock_open(read_data = current + '/test_data/bwa').return_value
        bwa = 'bwa_'
        samtools = 'samtools_'
        run_alignment(current + '/test_data/small_lambda_synth_taps_lambda_1.fq.gz', current + '/test_data/small_lambda_synth_taps_lambda_2.fq.gz', current + '/test_data/lambda_phage.fa', bwa, samtools, current + '/test_data/', 'mCtoT', 'BAM', 1, False, 1, 19, 100, 100, 1.5, 20, 500, 0.5, 0, 50, False, False, 1, 4, [6,6], [1,1], [5,5], 17, 'null', False)
        self.assertEqual(bwa, 'bwa_')
        self.assertEqual(samtools, 'samtools_')
        remove = 'rm {}'.format(current + '/test_data/small_lambda_synth_taps_lambda_mCtoT.bam')
        subprocess.Popen(remove, shell=True)
        
        
    def test_mock_taps_single_end(self):
        """Tests whether the aligner function will run with single-end TAPS reads."""
        package = version_testing_builtin()
        samtools = mock_open(read_data = current + '/test_data/samtools').return_value
        bwa = mock_open(read_data = current + '/test_data/bwa').return_value
        bwa = 'bwa_'
        samtools = 'samtools_'
        run_alignment(current + '/test_data/small_lambda_synth_taps_lambda_1.fq.gz', '', current + '/test_data/lambda_phage.fa', bwa, samtools, current + '/test_data/', 'mCtoT', 'BAM', 1, False, 1, 19, 100, 100, 1.5, 20, 500, 0.5, 0, 50, False, False, 1, 4, [6,6], [1,1], [5,5], 17, 'null', True)
        self.assertEqual(bwa, 'bwa_')
        self.assertEqual(samtools, 'samtools_')
        remove = 'rm {}'.format(current + '/test_data/small_lambda_synth_taps_lambda_mCtoT.bam')
        subprocess.Popen(remove, shell=True)
        
        
    def test_mock_wgbs_pair_end(self):
        """Tests whether the aligner function will run with pair-end WGBS reads."""
        package = version_testing_builtin()
        samtools = mock_open(read_data = current + '/test_data/samtools').return_value
        bwa = mock_open(read_data = current + '/test_data/bwameth.py').return_value
        bwa = 'bwameth.py_'
        samtools = 'samtools_'
        run_alignment(current + '/test_data/small_lambda_synth_wgbs_lambda_1.fq.gz', current + '/test_data/small_lambda_synth_wgbs_lambda_2.fq.gz', current + '/test_data/lambda_phage.fa', bwa, samtools, current + '/test_data/', 'CtoT', 'BAM', 1, False, 1, 19, 100, 100, 1.5, 20, 500, 0.5, 0, 50, False, False, 1, 4, [6,6], [1,1], [5,5], 17, 'null', False)
        self.assertEqual(bwa, 'bwameth.py_')
        self.assertEqual(samtools, 'samtools_')
        remove = 'rm {}'.format(current + '/test_data/small_lambda_synth_wgbs_lambda_CtoT.bam')
        subprocess.Popen(remove, shell=True)
        
        

if __name__ == '__main__':
    unittest.main()
