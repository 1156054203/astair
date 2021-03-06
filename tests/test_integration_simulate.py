import pdb
import csv
import gzip
import unittest
import subprocess
from os import path

from astair.simulator import modification_simulator
from astair.caller import cytosine_modification_finder


current = path.abspath(path.dirname(__file__))

class SimulateOutputTest(unittest.TestCase):
    """Tests whether the same positions will be modified for TAPS and WGBS set with the same seed.
    Also, ensures that the modified positions are correct as modification level and context."""


    def test_simulate_CHG_taps(self):
        """Tests whether CHG positions will be fully modified in TAPS."""
        modification_simulator(current + '/test_data/lambda_phage.fa', None, None, 75, current + '/test_data/small_lambda.bam', 'bam', 'mCtoT', 100, 'directional', None, 'CHG', None,  1,  (None, None, None), False, None, 0.3, 0.1, 1, False, current + '/test_data/', False, False, False, 20, 0, 0, True, False, True, True, 250)
        data_generated = list()
        with gzip.open(current + '/test_data/small_lambda_mCtoT_100_CHG_modified_positions_information.txt.gz','rt') as simulated_file:
            mod_reader = csv.reader(simulated_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated, [('__________________________________________________________________________________________________',),
                                          ('Absolute modified positions: 44   |   Percentage to all positions of the desired context: 0.682 %',), ('__________________________________________________________________________________________________',),
                                          ('lambda', '5', '6'), ('lambda', '16', '17'), ('lambda', '41', '42'), ('lambda', '44', '45'), ('lambda', '57', '58'),
                                          ('lambda', '102', '103'), ('lambda', '104', '105'), ('lambda', '121', '122'), ('lambda', '123', '124'), ('lambda', '127', '128'), ('lambda', '129', '130'), ('lambda', '150', '151'), ('lambda', '152', '153'), ('lambda', '166', '167'), ('lambda', '168', '169'), ('lambda', '176', '177'), ('lambda', '208', '209'), ('lambda', '210', '211'), ('lambda', '211', '212'), ('lambda', '213', '214'), ('lambda', '215', '216'), ('lambda', '217', '218'), ('lambda', '227', '228'), ('lambda', '237', '238'), ('lambda', '247', '248'), ('lambda', '249', '250'), ('lambda', '252', '253'), ('lambda', '254', '255'), ('lambda', '256', '257'), ('lambda', '258', '259'), ('lambda', '262', '263'), ('lambda', '264', '265'), ('lambda', '272', '273'), ('lambda', '277', '278'), ('lambda', '279', '280'), ('lambda', '287', '288'), ('lambda', '317', '318'), ('lambda', '319', '320'), ('lambda', '320', '321'), ('lambda', '323', '324'), ('lambda', '341', '342'), ('lambda', '353', '354'), ('lambda', '355', '356'), ('lambda', '373', '374')])
        cytosine_modification_finder(current + '/test_data/small_lambda_mCtoT_100_CHG.bam', None, 'none', current + '/test_data/lambda_phage.fa',
                                     'CHG', False, False, 13, None, 'directional','mCtoT', 0, 0, True, False, True, True, 250, None,
                                     1, current + '/test_data/', False, False, False,0, 0, 0)
        data_generated = list()
        with open(current + '/test_data/small_lambda_mCtoT_100_CHG_mCtoT_CHG.mods','rt') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('#CHROM', 'START', 'END', 'MOD_LEVEL', 'MOD', 'UNMOD', 'REF', 'ALT', 'SPECIFIC_CONTEXT', 'CONTEXT', 'SNV', 'TOTAL_DEPTH'), ('lambda', '5', '6', '1.0', '4', '0', 'G', 'A', 'CCG', 'CHG', 'No', '9'), ('lambda', '16', '17', '1.0', '18', '0', 'G', 'A', 'CCG', 'CHG', 'No', '33'), ('lambda', '41', '42', '1.0', '37', '0', 'C', 'T', 'CCG', 'CHG', 'No', '95'), ('lambda', '44', '45', '1.0', '59', '0', 'G', 'A', 'CCG', 'CHG', 'No', '100'), ('lambda', '57', '58', '1.0', '48', '0', 'C', 'T', 'CCG', 'CHG', 'No', '121'), ('lambda', '102', '103', '1.0', '70', '0', 'C', 'T', 'CTG', 'CHG', 'No', '150'), ('lambda', '104', '105', '1.0', '80', '0', 'G', 'A', 'CAG', 'CHG', 'No', '149'), ('lambda', '121', '122', '1.0', '75', '0', 'C', 'T', 'CAG', 'CHG', 'No', '150'), ('lambda', '123', '124', '1.0', '73', '0', 'G', 'A', 'CTG', 'CHG', 'No', '150'), ('lambda', '127', '128', '1.0', '80', '0', 'C', 'T', 'CTG', 'CHG', 'No', '155'), ('lambda', '129', '130', '1.0', '72', '0', 'G', 'A', 'CAG', 'CHG', 'No', '156'), ('lambda', '150', '151', '1.0', '95', '0', 'C', 'T', 'CTG', 'CHG', 'No', '164'), ('lambda', '152', '153', '1.0', '68', '0', 'G', 'A', 'CAG', 'CHG', 'No', '163'), ('lambda', '166', '167', '1.0', '100', '0', 'C', 'T', 'CTG', 'CHG', 'No', '165'), ('lambda', '168', '169', '1.0', '64', '0', 'G', 'A', 'CAG', 'CHG', 'No', '161'), ('lambda', '176', '177', '1.0', '102', '0', 'C', 'T', 'CCG', 'CHG', 'No', '183'), ('lambda', '208', '209', '1.0', '118', '0', 'C', 'T', 'CAG', 'CHG', 'No', '249'), ('lambda', '210', '211', '1.0', '130', '0', 'G', 'A', 'CTG', 'CHG', 'No', '248'), ('lambda', '211', '212', '1.0', '115', '0', 'C', 'T', 'CTG', 'CHG', 'No', '246'), ('lambda', '213', '214', '1.0', '131', '0', 'G', 'A', 'CAG', 'CHG', 'No', '248'), ('lambda', '215', '216', '1.0', '117', '0', 'C', 'T', 'CTG', 'CHG', 'No', '249'), ('lambda', '217', '218', '1.0', '129', '0', 'G', 'A', 'CAG', 'CHG', 'No', '246'), ('lambda', '227', '228', '1.0', '131', '0', 'G', 'A', 'CCG', 'CHG', 'No', '250'), ('lambda', '237', '238', '1.0', '114', '0', 'C', 'T', 'CCG', 'CHG', 'No', '246'), ('lambda', '247', '248', '1.0', '115', '0', 'C', 'T', 'CAG', 'CHG', 'No', '247'), ('lambda', '249', '250', '1.0', '126', '0', 'G', 'A', 'CTG', 'CHG', 'No', '244'), ('lambda', '252', '253', '1.0', '115', '0', 'C', 'T', 'CTG', 'CHG', 'No', '240'), ('lambda', '254', '255', '1.0', '122', '0', 'G', 'A', 'CAG', 'CHG', 'No', '239'), ('lambda', '256', '257', '1.0', '117', '0', 'C', 'T', 'CAG', 'CHG', 'No', '238'), ('lambda', '258', '259', '1.0', '119', '0', 'G', 'A', 'CTG', 'CHG', 'No', '240'), ('lambda', '262', '263', '1.0', '117', '0', 'C', 'T', 'CAG', 'CHG', 'No', '227'), ('lambda', '264', '265', '1.0', '105', '0', 'G', 'A', 'CTG', 'CHG', 'No', '218'), ('lambda', '272', '273', '1.0', '121', '0', 'C', 'T', 'CCG', 'CHG', 'No', '226'), ('lambda', '277', '278', '1.0', '119', '0', 'C', 'T', 'CTG', 'CHG', 'No', '216'), ('lambda', '279', '280', '1.0', '99', '0', 'G', 'A', 'CAG', 'CHG', 'No', '216'), ('lambda', '287', '288', '1.0', '104', '0', 'G', 'A', 'CCG', 'CHG', 'No', '227'), ('lambda', '317', '318', '1.0', '117', '0', 'C', 'T', 'CTG', 'CHG', 'No', '210'), ('lambda', '319', '320', '1.0', '92', '0', 'G', 'A', 'CAG', 'CHG', 'No', '207'), ('lambda', '320', '321', '1.0', '115', '0', 'C', 'T', 'CCG', 'CHG', 'No', '205'), ('lambda', '323', '324', '1.0', '107', '0', 'C', 'T', 'CCG', 'CHG', 'No', '194'), ('lambda', '341', '342', '1.0', '75', '0', 'C', 'T', 'CCG', 'CHG', 'No', '136'), ('lambda', '353', '354', '1.0', '50', '0', 'C', 'T', 'CTG', 'CHG', 'No', '90'), ('lambda', '355', '356', '1.0', '35', '0', 'G', 'A', 'CAG', 'CHG', 'No', '80')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_mCtoT_100_CHG*')
        subprocess.Popen(remove, shell=True)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
   
   
    def test_simulate_CHG_wgbs(self):
        """Tests whether CHG positions will be fully modified in WGBS."""
        modification_simulator(current + '/test_data/lambda_phage.fa', None, None, 75, current + '/test_data/small_lambda.bam', 'bam', 'CtoT', 100, 'directional', None, 'CHG', None,  1,  (None, None, None), False, None, 0.3, 0.1, 1, False, current + '/test_data/', False, False, False, 20, 0, 0, True, False, True, True, 250)
        data_generated = list()
        with gzip.open(current + '/test_data/small_lambda_CtoT_100_CHG_modified_positions_information.txt.gz','rt') as simulated_file:
            mod_reader = csv.reader(simulated_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated, [('__________________________________________________________________________________________________',),
                                          ('Absolute modified positions: 44   |   Percentage to all positions of the desired context: 0.682 %',), ('__________________________________________________________________________________________________',),
                                          ('lambda', '5', '6'), ('lambda', '16', '17'), ('lambda', '41', '42'), ('lambda', '44', '45'), ('lambda', '57', '58'),
                                          ('lambda', '102', '103'), ('lambda', '104', '105'), ('lambda', '121', '122'), ('lambda', '123', '124'), ('lambda', '127', '128'), ('lambda', '129', '130'), ('lambda', '150', '151'), ('lambda', '152', '153'), ('lambda', '166', '167'), ('lambda', '168', '169'), ('lambda', '176', '177'), ('lambda', '208', '209'), ('lambda', '210', '211'), ('lambda', '211', '212'), ('lambda', '213', '214'), ('lambda', '215', '216'), ('lambda', '217', '218'), ('lambda', '227', '228'), ('lambda', '237', '238'), ('lambda', '247', '248'), ('lambda', '249', '250'), ('lambda', '252', '253'), ('lambda', '254', '255'), ('lambda', '256', '257'), ('lambda', '258', '259'), ('lambda', '262', '263'), ('lambda', '264', '265'), ('lambda', '272', '273'), ('lambda', '277', '278'), ('lambda', '279', '280'), ('lambda', '287', '288'), ('lambda', '317', '318'), ('lambda', '319', '320'), ('lambda', '320', '321'), ('lambda', '323', '324'), ('lambda', '341', '342'), ('lambda', '353', '354'), ('lambda', '355', '356'), ('lambda', '373', '374')])
        cytosine_modification_finder(current + '/test_data/small_lambda_CtoT_100_CHG.bam', None, 'none', current + '/test_data/lambda_phage.fa', 'CHG', False, False, 13, None, 'directional',  'CtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False, False, 0, 0, 0)
        data_generated = list()
        with open(current + '/test_data/small_lambda_CtoT_100_CHG_CtoT_CHG.mods','rt') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated[:-1], [('#CHROM', 'START', 'END', 'MOD_LEVEL', 'MOD', 'UNMOD', 'REF', 'ALT', 'SPECIFIC_CONTEXT', 'CONTEXT', 'SNV', 'TOTAL_DEPTH'), ('lambda', '5', '6', '1.0', '4', '0', 'G', 'A', 'CCG', 'CHG', 'No', '9'), ('lambda', '16', '17', '1.0', '18', '0', 'G', 'A', 'CCG', 'CHG', 'No', '33'), ('lambda', '41', '42', '1.0', '37', '0', 'C', 'T', 'CCG', 'CHG', 'No', '95'), ('lambda', '44', '45', '1.0', '59', '0', 'G', 'A', 'CCG', 'CHG', 'No', '100'), ('lambda', '57', '58', '1.0', '48', '0', 'C', 'T', 'CCG', 'CHG', 'No', '121'), ('lambda', '102', '103', '1.0', '70', '0', 'C', 'T', 'CTG', 'CHG', 'No', '150'), ('lambda', '104', '105', '1.0', '80', '0', 'G', 'A', 'CAG', 'CHG', 'No', '149'), ('lambda', '121', '122', '1.0', '75', '0', 'C', 'T', 'CAG', 'CHG', 'No', '150'), ('lambda', '123', '124', '1.0', '73', '0', 'G', 'A', 'CTG', 'CHG', 'No', '150'), ('lambda', '127', '128', '1.0', '80', '0', 'C', 'T', 'CTG', 'CHG', 'No', '155'), ('lambda', '129', '130', '1.0', '72', '0', 'G', 'A', 'CAG', 'CHG', 'No', '156'), ('lambda', '150', '151', '1.0', '95', '0', 'C', 'T', 'CTG', 'CHG', 'No', '164'), ('lambda', '152', '153', '1.0', '68', '0', 'G', 'A', 'CAG', 'CHG', 'No', '163'), ('lambda', '166', '167', '1.0', '100', '0', 'C', 'T', 'CTG', 'CHG', 'No', '165'), ('lambda', '168', '169', '1.0', '64', '0', 'G', 'A', 'CAG', 'CHG', 'No', '161'), ('lambda', '176', '177', '1.0', '102', '0', 'C', 'T', 'CCG', 'CHG', 'No', '183'), ('lambda', '208', '209', '1.0', '118', '0', 'C', 'T', 'CAG', 'CHG', 'No', '249'), ('lambda', '210', '211', '1.0', '130', '0', 'G', 'A', 'CTG', 'CHG', 'No', '248'), ('lambda', '211', '212', '1.0', '115', '0', 'C', 'T', 'CTG', 'CHG', 'No', '246'), ('lambda', '213', '214', '1.0', '131', '0', 'G', 'A', 'CAG', 'CHG', 'No', '248'), ('lambda', '215', '216', '0.991', '116', '1', 'C', 'T', 'CTG', 'CHG', 'No', '249'), ('lambda', '217', '218', '1.0', '129', '0', 'G', 'A', 'CAG', 'CHG', 'No', '246'), ('lambda', '227', '228', '1.0', '131', '0', 'G', 'A', 'CCG', 'CHG', 'No', '250'), ('lambda', '237', '238', '1.0', '114', '0', 'C', 'T', 'CCG', 'CHG', 'No', '246'), ('lambda', '247', '248', '1.0', '115', '0', 'C', 'T', 'CAG', 'CHG', 'No', '247'), ('lambda', '249', '250', '0.992', '125', '1', 'G', 'A', 'CTG', 'CHG', 'No', '244'), ('lambda', '252', '253', '1.0', '115', '0', 'C', 'T', 'CTG', 'CHG', 'No', '240'), ('lambda', '254', '255', '0.992', '121', '1', 'G', 'A', 'CAG', 'CHG', 'No', '239'), ('lambda', '256', '257', '1.0', '117', '0', 'C', 'T', 'CAG', 'CHG', 'No', '238'), ('lambda', '258', '259', '0.992', '118', '1', 'G', 'A', 'CTG', 'CHG', 'No', '240'), ('lambda', '262', '263', '1.0', '117', '0', 'C', 'T', 'CAG', 'CHG', 'No', '227'), ('lambda', '264', '265', '1.0', '105', '0', 'G', 'A', 'CTG', 'CHG', 'No', '218'), ('lambda', '272', '273', '1.0', '121', '0', 'C', 'T', 'CCG', 'CHG', 'No', '226'), ('lambda', '277', '278', '1.0', '119', '0', 'C', 'T', 'CTG', 'CHG', 'No', '216'), ('lambda', '279', '280', '1.0', '99', '0', 'G', 'A', 'CAG', 'CHG', 'No', '216'), ('lambda', '287', '288', '1.0', '104', '0', 'G', 'A', 'CCG', 'CHG', 'No', '227'), ('lambda', '317', '318', '1.0', '117', '0', 'C', 'T', 'CTG', 'CHG', 'No', '210'), ('lambda', '319', '320', '1.0', '92', '0', 'G', 'A', 'CAG', 'CHG', 'No', '207'), ('lambda', '320', '321', '0.991', '114', '1', 'C', 'T', 'CCG', 'CHG', 'No', '205'), ('lambda', '323', '324', '1.0', '107', '0', 'C', 'T', 'CCG', 'CHG', 'No', '194'), ('lambda', '341', '342', '1.0', '75', '0', 'C', 'T', 'CCG', 'CHG', 'No', '136'), ('lambda', '353', '354', '1.0', '50', '0', 'C', 'T', 'CTG', 'CHG', 'No', '90'), ('lambda', '355', '356', '1.0', '35', '0', 'G', 'A', 'CAG', 'CHG', 'No', '80')])
        remove = 'rm {}'.format(current + '/test_data/small_lambda_CtoT_100_CHG*')
        subprocess.Popen(remove, shell=True)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        


    def test_reverse_CpG_taps(self):
        """Tests whether truely modified CpG positions will be will be reverted in TAPS."""
        modification_simulator(current + '/test_data/lambda_phage.fa', None, None, 80, current + '/test_data/small_real_taps_lambda_mCtoT.bam', 'bam', 'mCtoT', 100, 'directional', None, 'all', None,  1,  (None, None, None), False, None, 0.3, 0.1, 1, False, current + '/test_data/', False, False, False, 20, 0, 0, True, False, True, True, 250)
        modification_simulator(current + '/test_data/lambda_phage.fa', None, None, 80, current + '/test_data/small_real_taps_lambda_mCtoT_mCtoT_100_all.bam', 'bam', 'mCtoT', 100, 'directional', None, 'CpG', None,  1,  (None, None, None), False, None, 0.3, 0.1, 1, True, current + '/test_data/', False, False, False, 20, 0, 0, True, False, True, True, 250)
        data_generated = list()
        cytosine_modification_finder(current + '/test_data/small_real_taps_lambda_mCtoT_mCtoT_100_all_mCtoT_100_CpG_reversed.bam', None, 'none', current + '/test_data/lambda_phage.fa', 'CpG', False, False, 13, None, 'directional', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False, False, 0, 0, 0)
        data_generated = list()
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_mCtoT_100_all_mCtoT_100_CpG_reversed_mCtoT_CpG.mods','rt') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row[0:4])))
        self.assertEqual(data_generated[0:15],[('#CHROM', 'START', 'END', 'MOD_LEVEL'), ('lambda', '3', '4', '0.0'), ('lambda', '4', '5', '0.0'), ('lambda', '6', '7', '0.0'), ('lambda', '7', '8', '0.0'), ('lambda', '12', '13', '0.0'), ('lambda', '13', '14', '0.0'), ('lambda', '14', '15', '0.0'), ('lambda', '15', '16', '0.0'), ('lambda', '22', '23', '0.0'), ('lambda', '23', '24', '0.0'), ('lambda', '42', '43', '0.0'), ('lambda', '43', '44', '0.0'), ('lambda', '52', '53', '0.0'), ('lambda', '53', '54', '0.0')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT_mCtoT_100_all*')
        subprocess.Popen(remove, shell=True)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        


    def test_by_list_CpG_taps(self):
        """Tests whether truely modified CpG positions will be modified in TAPS if a list of positions is given."""
        cytosine_modification_finder(current + '/test_data/small_real_taps_lambda_mCtoT.bam', None, 'none', current + '/test_data/lambda_phage.fa', 'CpG', False, False, 13, None, 'directional', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False, False,0, 0, 0)
        modification_simulator(current + '/test_data/lambda_phage.fa', None, None, 80, current + '/test_data/small_real_taps_lambda_mCtoT.bam', 'bam', 'mCtoT', 100, 'directional', None, 'all', None,  1,  (None, None, None), False, None, 0.3, 0.1, 1, False, current + '/test_data/', False, False, False, 20, 0, 0, True, False, True, True, 250)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        modification_simulator(current + '/test_data/lambda_phage.fa', None, None, 80, current + '/test_data/small_real_taps_lambda_mCtoT.bam', 'bam', 'mCtoT', 100, 'directional',
              current + '/test_data/small_real_taps_lambda_mCtoT_mCtoT_CpG.mods', 'all',None,  1,  (None, None, None), False, None, 0.3, 0.1, 1, False, current + '/test_data/', False, False, False, 20, 0, 0, True, False, True, True, 250)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        cytosine_modification_finder(current + '/test_data/small_real_taps_lambda_mCtoT_mCtoT_user_provided_list_all.bam', None, 'none', current + '/test_data/lambda_phage.fa', 'CpG', False, False, 13, None, 'directional', 'mCtoT', 0, 0, True, False, True, True, 250, None, 1, current + '/test_data/', False, False, False,0, 0, 0)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)
        data_generated = list()
        modification_simulator(current + '/test_data/lambda_phage.fa', None, None, 80, current + '/test_data/small_real_taps_lambda_mCtoT.bam', 'bam', 'mCtoT', 100, 'directional', None, 'all', None,  1,  (None, None, None), False, None, 0.3, 0.1, 1, False, current + '/test_data/', False, False, False, 20, 0, 0, True, False, True, True, 250)
        with open(current + '/test_data/small_real_taps_lambda_mCtoT_mCtoT_user_provided_list_all_mCtoT_CpG.stats','rt') as call_file:
            mod_reader = csv.reader(call_file, delimiter='\t', lineterminator='\n')
            for row in mod_reader:
                data_generated.append(tuple((row)))
        self.assertEqual(data_generated,[('#CONTEXT', 'SPECIFIC_CONTEXT', 'MEAN_MODIFICATION_RATE_PERCENT', 'TOTAL_POSITIONS', 'COVERED_POSITIONS', 'MODIFIED', 'UNMODIFIED'), ('CpG', '*', '88.946', '6225', '48', '1722', '214'), ('*', 'CGA', '97.593', '1210', '12', '527', '13'), ('*', 'CGC', '66.667', '1730', '13', '362', '181'), ('*', 'CGG', '97.065', '1847', '13', '496', '15'), ('*', 'CGT', '98.538', '1438', '10', '337', '5')])
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT_mCtoT_user*')
        subprocess.Popen(remove, shell=True)
        remove = 'rm {}'.format(current + '/test_data/small_real_taps_lambda_mCtoT_mCtoT_CpG*')
        subprocess.Popen(remove, shell=True)
        if path.isfile(current + '/test_data/lambda_phage_.fa'):
            subprocess.Popen('gzip {}'.format(current + '/test_data/lambda_phage_.fa'), shell=True)



if __name__ == '__main__':
    unittest.main()

