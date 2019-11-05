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
        cytosine_modification_finder(current + '/test_data/small_lambda.bam', None, None, current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'directional', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False, False)
        data_generated = list()
        with open(current + '/test_data/small_lambda_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('#CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '63.338', '6225', '46', '2023', '1171'), ('*', 'CGA', '53.383', '1210', '12', '434', '379'), ('*', 'CGC', '96.688', '1730', '11', '613', '21'), ('*', 'CGG', '54.018', '1847', '13', '531', '452'), ('*', 'CGT', '58.246', '1438', '10', '445', '319'), ('CHG', '*', '0.125', '6451', '44', '5', '3998'), ('*', 'CAG', '0.07', '2302', '15', '1', '1433'), ('*', 'CCG', '0.097', '1847', '13', '1', '1034'), ('*', 'CTG', '0.196', '2302', '16', '3', '1531'), ('CHH', '*', '25.23', '11503', '81', '1701', '5041'), ('*', 'CTT', '18.36', '1349', '15', '215', '956'), ('*', 'CAT', '29.035', '1802', '14', '349', '853'), ('*', 'CCT', '29.564', '1182', '14', '332', '791'), ('*', 'CTA', '0.0', '501', '1', '0', '23'), ('*', 'CAA', '41.746', '1432', '7', '263', '367'), ('*', 'CCA', '31.18', '1610', '8', '251', '554'), ('*', 'CTC', '8.447', '1116', '10', '68', '737'), ('*', 'CAC', '40.594', '1474', '5', '205', '300'), ('*', 'CCC', '3.766', '1037', '7', '18', '460')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)
        
        
    def test_call_default_taps_known_snps(self):
        """Looks for TAPS modified cytosine positions and compares their modification level with
        the expected one by context and labels known SNPs."""
        cytosine_modification_finder(current + '/test_data/small_real_taps_chr10:20000-60000_pos.bam', current + '/test_data/GRCh38p7_common_snps_sample.vcf.gz', None, current + '/test_data/hg38_chr10_20000-60000.fa', 'all', False, False, 13, None, 'directional', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False, False)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_chr10:20000-60000_pos_mCtoT_all.mods','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[4500:4530], [('chr10', '25922', '25923', '0.0', '0', '13', 'C', 'T', 'CTC', 'CHH', 'WGS_known', '27'), ('chr10', '25924', '25925', '0.0', '0', '14', 'C', 'T', 'CAC', 'CHH', 'No', '28'), ('chr10', '25926', '25927', '0.0', '0', '13', 'C', 'T', 'CAC', 'CHH', 'No', '27'), ('chr10', '25928', '25929', '0.0', '0', '13', 'C', 'T', 'CAT', 'CHH', 'No', '27'), ('chr10', '25931', '25932', '0.0', '0', '11', 'C', 'T', 'CTC', 'CHH', 'No', '25'), ('chr10', '25933', '25934', '0.0', '0', '13', 'C', 'T', 'CTC', 'CHH', 'No', '27'), ('chr10', '25935', '25936', '0.0', '0', '13', 'C', 'T', 'CAT', 'CHH', 'No', '26'), ('chr10', '25940', '25941', '0.0', '0', '13', 'G', 'A', 'CAA', 'CHH', 'No', '25'), ('chr10', '25942', '25943', '0.0', '0', '14', 'C', 'T', 'CAG', 'CHG', 'No', '27'), ('chr10', '25944', '25945', '0.0', '0', '12', 'G', 'A', 'CTG', 'CHG', 'No', '24'), ('chr10', '25945', '25946', '0.0', '0', '13', 'C', 'T', 'CTT', 'CHH', 'WGS_known', '24'), ('chr10', '25948', '25949', '0.0', '0', '11', 'C', 'T', 'CCC', 'CHH', 'No', '22'), ('chr10', '25949', '25950', '0.0', '0', '13', 'C', 'T', 'CCC', 'CHH', 'No', '23'), ('chr10', '25950', '25951', '0.0', '0', '13', 'C', 'T', 'CCA', 'CHH', 'No', '23'), ('chr10', '25951', '25952', '0.0', '0', '12', 'C', 'T', 'CAA', 'CHH', 'No', '22'), ('chr10', '25955', '25956', '0.0', '0', '9', 'G', 'A', 'CTT', 'CHH', 'No', '21'), ('chr10', '25956', '25957', '0.0', '0', '10', 'C', 'T', 'CAG', 'CHG', 'No', '20'), ('chr10', '25958', '25959', '0.0', '0', '9', 'G', 'A', 'CTG', 'CHG', 'No', '22'), ('chr10', '25959', '25960', '0.0', '0', '9', 'G', 'A', 'CCT', 'CHH', 'No', '22'), ('chr10', '25960', '25961', '0.0', '0', '13', 'C', 'T', 'CTT', 'CHH', 'No', '22'), ('chr10', '25964', '25965', '0.0', '0', '9', 'G', 'A', 'CAA', 'CHH', 'No', '22'), ('chr10', '25965', '25966', '0.0', '0', '13', 'C', 'T', 'CAG', 'CHG', 'No', '22'), ('chr10', '25967', '25968', '0.0', '0', '9', 'G', 'A', 'CTG', 'CHG', 'No', '22'), ('chr10', '25969', '25970', '0.0', '0', '9', 'C', 'T', 'CCA', 'CHH', 'No', '18'), ('chr10', '25970', '25971', '0.0', '0', '11', 'C', 'T', 'CAG', 'CHG', 'No', '20'), ('chr10', '25972', '25973', '0.0', '0', '9', 'G', 'A', 'CTG', 'CHG', 'No', '20'), ('chr10', '25973', '25974', '0.0', '0', '9', 'G', 'A', 'CCT', 'CHH', 'No', '20'), ('chr10', '25974', '25975', '0.0', '0', '10', 'C', 'T', 'CTC', 'CHH', 'No', '19'), ('chr10', '25976', '25977', '0.0', '0', '11', 'C', 'T', 'CCC', 'CHH', 'No', '20'), ('chr10', '25977', '25978', '0.0', '0', '12', 'C', 'T', 'CCC', 'CHH', 'No', '20')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_chr10:20000-60000_pos_mCtoT_all.*')
        subprocess.Popen(remove, shell=True) 
        
          
        
    def test_call_taps_mark_ends(self):
        """Looks for TAPS modified cytosine positions and compares their modification level with
        the expected one by context for a sample with mark_ends issue."""
        cytosine_modification_finder(current + '/test_data/small_real_taps_chr10:20000-60000_pos.bam', None, None, current + '/test_data/hg38_chr10_20000-60000.fa', 'all', False, False, 13, None, 'directional', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False, False)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_chr10:20000-60000_pos_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:20], [('#CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '59.47', '804', '567', '1727', '1177'), ('*', 'CGA', '52.542', '170', '115', '341', '308'), ('*', 'CGC', '74.0', '191', '135', '444', '156'), ('*', 'CGG', '65.306', '233', '169', '512', '272'), ('*', 'CGT', '49.369', '210', '148', '430', '441'), ('CHG', '*', '0.174', '3583', '2464', '32', '18319'), ('*', 'CAG', '0.205', '1675', '1147', '18', '8760'), ('*', 'CCG', '0.38', '233', '170', '3', '787'), ('*', 'CTG', '0.125', '1675', '1147', '11', '8772'), ('CHH', '*', '0.174', '12229', '7772', '104', '59709'), ('*', 'CTT', '0.204', '1503', '900', '15', '7336'), ('*', 'CAT', '0.221', '1521', '878', '15', '6766'), ('*', 'CCT', '0.092', '1398', '966', '7', '7561'), ('*', 'CTA', '0.046', '952', '550', '2', '4330'), ('*', 'CAA', '0.027', '1626', '938', '2', '7312'), ('*', 'CCA', '0.237', '1606', '1073', '19', '8009'), ('*', 'CTC', '0.214', '1253', '843', '14', '6515'), ('*', 'CAC', '0.314', '1176', '748', '18', '5714'), ('*', 'CCC', '0.194', '1194', '876', '12', '6166')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_chr10:20000-60000_pos_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)     
        
    
    def test_call_default_taps_gzip(self):
        """Looks for TAPS modified cytosine positions and compares their modification level with
        the expected one by context with GZIP compressed fasta reference."""
        cytosine_modification_finder(current + '/test_data/small_lambda.bam', None, None, current + '/test_data/lambda_phage_.fa.gz', 'all', False, False, 13, None, 'directional', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False, False)
        data_generated = list()
        with open(current + '/test_data/small_lambda_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('#CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '63.338', '6225', '46', '2023', '1171'), ('*', 'CGA', '53.383', '1210', '12', '434', '379'), ('*', 'CGC', '96.688', '1730', '11', '613', '21'), ('*', 'CGG', '54.018', '1847', '13', '531', '452'), ('*', 'CGT', '58.246', '1438', '10', '445', '319'), ('CHG', '*', '0.125', '6451', '44', '5', '3998'), ('*', 'CAG', '0.07', '2302', '15', '1', '1433'), ('*', 'CCG', '0.097', '1847', '13', '1', '1034'), ('*', 'CTG', '0.196', '2302', '16', '3', '1531'), ('CHH', '*', '25.23', '11503', '81', '1701', '5041'), ('*', 'CTT', '18.36', '1349', '15', '215', '956'), ('*', 'CAT', '29.035', '1802', '14', '349', '853'), ('*', 'CCT', '29.564', '1182', '14', '332', '791'), ('*', 'CTA', '0.0', '501', '1', '0', '23'), ('*', 'CAA', '41.746', '1432', '7', '263', '367'), ('*', 'CCA', '31.18', '1610', '8', '251', '554'), ('*', 'CTC', '8.447', '1116', '10', '68', '737'), ('*', 'CAC', '40.594', '1474', '5', '205', '300'), ('*', 'CCC', '3.766', '1037', '7', '18', '460')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        
        
    def test_call_default_taps_SE(self):
        """Looks for TAPS modified cytosine positions and compares their modification level with
        the expected one by context."""
        cytosine_modification_finder(current + '/test_data/small_real_taps_lambda_mCtoT_SE.bam', None, None, current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'directional', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, True, False)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_SE_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:6],[('#CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '98.039', '6225', '103', '50', '1'), ('*', 'CGA', '100.0', '1210', '24', '12', '0'), ('*', 'CGC', '100.0', '1730', '21', '7', '0'), ('*', 'CGG', '94.118', '1847', '33', '16', '1'), ('*', 'CGT', '100.0', '1438', '25', '15', '0')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT_SE_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)
        

    def test_call_default_taps_reversed(self):
        """Looks for TAPS modified cytosine positions and compares their modification level with
        the expected one by context in TAPS data created with reversed directionality (i.e OT has G>A modifications and OB has C>T)."""
        cytosine_modification_finder(current + '/test_data/small_real_taps_lambda_mCtoT_reversed.bam', None, None, current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'reverse', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False, False)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_reversed_mCtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('#CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '90.877', '6225', '2121', '24793', '2489'), ('*', 'CGA', '93.626', '1210', '462', '7212', '491'), ('*', 'CGC', '86.154', '1730', '552', '5681', '913'), ('*', 'CGG', '91.339', '1847', '609', '6022', '571'), ('*', 'CGT', '91.959', '1438', '498', '5878', '514'), ('CHG', '*', '0.222', '6451', '2257', '71', '31942'), ('*', 'CAG', '0.164', '2302', '820', '22', '13370'), ('*', 'CCG', '0.475', '1847', '610', '30', '6286'), ('*', 'CTG', '0.154', '2302', '827', '19', '12286'), ('CHH', '*', '0.163', '11503', '4625', '160', '98077'), ('*', 'CTT', '0.179', '1349', '558', '21', '11739'), ('*', 'CAT', '0.101', '1802', '773', '18', '17763'), ('*', 'CCT', '0.14', '1182', '443', '13', '9249'), ('*', 'CTA', '0.091', '501', '260', '4', '4376'), ('*', 'CAA', '0.089', '1432', '630', '12', '13501'), ('*', 'CCA', '0.207', '1610', '606', '27', '12987'), ('*', 'CTC', '0.151', '1116', '444', '14', '9258'), ('*', 'CAC', '0.188', '1474', '562', '21', '11128'), ('*', 'CCC', '0.37', '1037', '349', '30', '8076')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT_reversed_mCtoT_all.*')
        subprocess.Popen(remove, shell=True)



    def test_call_default_wgbs(self):
        """Looks for WGBS modified cytosine positions and compares their modification level with
        the expected one by context."""
        cytosine_modification_finder(current + '/test_data/small_lambda.bam', None, None, current + '/test_data/lambda_phage.fa', 'all', False, False, 13, None, 'directional', 'CtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False, False)
        data_generated = list()
        with open(current + '/test_data/small_lambda_CtoT_all.stats','r') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('#CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '36.662', '6225', '46', '1171', '2023'), ('*', 'CGA', '46.617', '1210', '12', '379', '434'), ('*', 'CGC', '3.312', '1730', '11', '21', '613'), ('*', 'CGG', '45.982', '1847', '13', '452', '531'), ('*', 'CGT', '41.754', '1438', '10', '319', '445'), ('CHG', '*', '99.875', '6451', '44', '3998', '5'), ('*', 'CAG', '99.93', '2302', '15', '1433', '1'), ('*', 'CCG', '99.903', '1847', '13', '1034', '1'), ('*', 'CTG', '99.804', '2302', '16', '1531', '3'), ('CHH', '*', '74.77', '11503', '81', '5041', '1701'), ('*', 'CTT', '81.64', '1349', '15', '956', '215'), ('*', 'CAT', '70.965', '1802', '14', '853', '349'), ('*', 'CCT', '70.436', '1182', '14', '791', '332'), ('*', 'CTA', '100.0', '501', '1', '23', '0'), ('*', 'CAA', '58.254', '1432', '7', '367', '263'), ('*', 'CCA', '68.82', '1610', '8', '554', '251'), ('*', 'CTC', '91.553', '1116', '10', '737', '68'), ('*', 'CAC', '59.406', '1474', '5', '300', '205'), ('*', 'CCC', '96.234', '1037', '7', '460', '18')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_CtoT_all.*')
        subprocess.Popen(remove, shell=True)
        
   
        
if __name__ == '__main__':
    unittest.main()
