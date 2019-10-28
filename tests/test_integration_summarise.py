import pdb
import csv
import sys
import gzip
import unittest
import subprocess
from os import path

from astair.summary import read_summariser

current = path.abspath(path.dirname(__file__))


class FindOutputTest(unittest.TestCase):
    """Tests whether the same cytosine positions will be called as modified for both TAPS and WGBS."""

    def test_taps_default_CpG(self):
        """Looks for cytosine positions and compares them with the expected ones by context."""
        read_summariser(current + '/test_data/small_real_taps_chr10:20000-60000_pos.bam', current + '/test_data/hg38_chr10_20000-60000.fa',  current + '/test_data/GRCh38p7_common_snps_sample.vcf.gz', 'CpG', None, 'directional', 'mCtoT', (None, None, None), 10, 1, None, 1,  current + '/test_data/', False, False)
        data_generated = list()
        if sys.version[0] == '3':
            variants_file = gzip.open(current + '/test_data/small_real_taps_chr10:20000-60000_pos_mCtoT_CpG_read_summary.txt.gz','rt') 
        elif sys.version[0] == '2':
            file_ = subprocess.Popen('gunzip {}'.format(current + '/test_data/small_real_taps_chr10:20000-60000_pos_mCtoT_CpG_read_summary.txt.gz'), shell=True)
            exit_code = file_.wait()
            if exit_code == 0:
                variants_file = open(current + '/test_data/small_real_taps_chr10:20000-60000_pos_mCtoT_CpG_read_summary.txt', 'r')
        for row in variants_file.readlines():
            data_generated.append(row.splitlines()[0].split('\t'))
        self.assertEqual(data_generated[0:20], [['#CHROM', 'START', 'END', 'READ_NAME', 'POSSIBLE_MODIFICATION', 'FLAG', 'CONTEXT', 'SPECIFIC_CONTEXT', 'BQ', 'MAPQ', 'KNOWN_SNP'], ['chr10', '881', '882', 'K00198:385:H32VHBBXY:8:2210:16224:28833', '0', '163', 'CpG', 'CGA', '41', '40', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:2210:16224:28833', '0', '163', 'CpG', 'CGT', '41', '40', '*'], ['chr10', '881', '882', 'K00198:385:H32VHBBXY:8:1110:12124:16418', '0', '163', 'CpG', 'CGA', '37', '40', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:1110:12124:16418', '1', '163', 'CpG', 'CGT', '41', '40', '*'], ['chr10', '881', '882', 'K00198:385:H32VHBBXY:8:1218:4219:20586', '0', '163', 'CpG', 'CGA', '41', '40', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:1218:4219:20586', '1', '163', 'CpG', 'CGT', '41', '40', '*'], ['chr10', '881', '882', 'K00198:385:H32VHBBXY:8:2126:16021:40508', '0', '163', 'CpG', 'CGA', '37', '3', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:2126:16021:40508', '0', '163', 'CpG', 'CGT', '41', '3', '*'], ['chr10', '975', '976', 'K00198:385:H32VHBBXY:8:2126:16021:40508', '0', '163', 'CpG', 'CGG', '27', '3', '*'], ['chr10', '976', '977', 'K00198:385:H32VHBBXY:8:2126:16021:40508', '0', '163', 'CpG', 'CGA', '27', '3', '*'], ['chr10', '881', '882', 'K00198:385:H32VHBBXY:8:1108:24038:6888', '0', '163', 'CpG', 'CGA', '12', '40', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:1108:24038:6888', '0', '163', 'CpG', 'CGT', '27', '40', '*'], ['chr10', '975', '976', 'K00198:385:H32VHBBXY:8:1108:24038:6888', '1', '163', 'CpG', 'CGG', '32', '40', '*'], ['chr10', '976', '977', 'K00198:385:H32VHBBXY:8:1108:24038:6888', '1', '163', 'CpG', 'CGA', '12', '40', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:1211:18781:38785', '0', '99', 'CpG', 'CGT', '32', '40', '*'], ['chr10', '975', '976', 'K00198:385:H32VHBBXY:8:1211:18781:38785', '0', '99', 'CpG', 'CGG', '41', '40', '*'], ['chr10', '976', '977', 'K00198:385:H32VHBBXY:8:1211:18781:38785', '0', '99', 'CpG', 'CGA', '41', '40', '*'], ['chr10', '975', '976', 'K00198:385:H32VHBBXY:8:2223:2595:20093', '0', '99', 'CpG', 'CGG', '12', '27', '*'], ['chr10', '976', '977', 'K00198:385:H32VHBBXY:8:2223:2595:20093', '0', '99', 'CpG', 'CGA', '22', '27', '*']])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_chr10:20000-60000_pos_mCtoT_CpG_read_summary*')
        subprocess.Popen(remove, shell=True)
        

    
    def test_wgbs_default_CpG(self):
        """Looks for cytosine positions and compares them with the expected ones by context."""
        read_summariser(current + '/test_data/small_real_taps_chr10:20000-60000_pos.bam', current + '/test_data/hg38_chr10_20000-60000.fa',  current + '/test_data/GRCh38p7_common_snps_sample.vcf.gz', 'CpG', None, 'directional', 'CtoT', (None, None, None), 10, 1, None, 1,  current + '/test_data/', False, False)
        data_generated = list()
        if sys.version[0] == '3':
            variants_file = gzip.open(current + '/test_data/small_real_taps_chr10:20000-60000_pos_CtoT_CpG_read_summary.txt.gz','rt') 
        elif sys.version[0] == '2':
            file_ = subprocess.Popen('gunzip {}'.format(current + '/test_data/small_real_taps_chr10:20000-60000_pos_CtoT_CpG_read_summary.txt.gz'), shell=True)
            exit_code = file_.wait()
            if exit_code == 0:
                variants_file = open(current + '/test_data/small_real_taps_chr10:20000-60000_pos_CtoT_CpG_read_summary.txt', 'r')
        for row in variants_file.readlines():
            data_generated.append(row.splitlines()[0].split('\t'))
        self.assertEqual(data_generated[0:20],[['#CHROM', 'START', 'END', 'READ_NAME', 'POSSIBLE_MODIFICATION', 'FLAG', 'CONTEXT', 'SPECIFIC_CONTEXT', 'BQ', 'MAPQ', 'KNOWN_SNP'], ['chr10', '881', '882', 'K00198:385:H32VHBBXY:8:2210:16224:28833', '0', '163', 'CpG', 'CGA', '41', '40', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:2210:16224:28833', '1', '163', 'CpG', 'CGT', '41', '40', '*'], ['chr10', '881', '882', 'K00198:385:H32VHBBXY:8:1110:12124:16418', '0', '163', 'CpG', 'CGA', '37', '40', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:1110:12124:16418', '0', '163', 'CpG', 'CGT', '41', '40', '*'], ['chr10', '881', '882', 'K00198:385:H32VHBBXY:8:1218:4219:20586', '0', '163', 'CpG', 'CGA', '41', '40', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:1218:4219:20586', '0', '163', 'CpG', 'CGT', '41', '40', '*'], ['chr10', '881', '882', 'K00198:385:H32VHBBXY:8:2126:16021:40508', '0', '163', 'CpG', 'CGA', '37', '3', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:2126:16021:40508', '1', '163', 'CpG', 'CGT', '41', '3', '*'], ['chr10', '975', '976', 'K00198:385:H32VHBBXY:8:2126:16021:40508', '0', '163', 'CpG', 'CGG', '27', '3', '*'], ['chr10', '976', '977', 'K00198:385:H32VHBBXY:8:2126:16021:40508', '1', '163', 'CpG', 'CGA', '27', '3', '*'], ['chr10', '881', '882', 'K00198:385:H32VHBBXY:8:1108:24038:6888', '1', '163', 'CpG', 'CGA', '12', '40', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:1108:24038:6888', '1', '163', 'CpG', 'CGT', '27', '40', '*'], ['chr10', '975', '976', 'K00198:385:H32VHBBXY:8:1108:24038:6888', '0', '163', 'CpG', 'CGG', '32', '40', '*'], ['chr10', '976', '977', 'K00198:385:H32VHBBXY:8:1108:24038:6888', '0', '163', 'CpG', 'CGA', '12', '40', '*'], ['chr10', '882', '883', 'K00198:385:H32VHBBXY:8:1211:18781:38785', '0', '99', 'CpG', 'CGT', '32', '40', '*'], ['chr10', '975', '976', 'K00198:385:H32VHBBXY:8:1211:18781:38785', '1', '99', 'CpG', 'CGG', '41', '40', '*'], ['chr10', '976', '977', 'K00198:385:H32VHBBXY:8:1211:18781:38785', '0', '99', 'CpG', 'CGA', '41', '40', '*'], ['chr10', '975', '976', 'K00198:385:H32VHBBXY:8:2223:2595:20093', '1', '99', 'CpG', 'CGG', '12', '27', '*'], ['chr10', '976', '977', 'K00198:385:H32VHBBXY:8:2223:2595:20093', '0', '99', 'CpG', 'CGA', '22', '27', '*']])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_chr10:20000-60000_pos_CtoT_CpG_read_summary*')
        subprocess.Popen(remove, shell=True)

            
            
if __name__ == '__main__':
    unittest.main()
