import sys
import unittest

from astair import mbias as mbias

class MbiasContextSearchTest(unittest.TestCase):
    """Tests whether given the read flag and the method, the correct bases for modification will be found.""" 
        
    def test_mbias_calculator_taps_read1(self):
        read1_mods_CHH, read1_mods_CHG, read1_mods_CpG = mbias.initialise_data_counters(30)
        read1_umod_CHH, read1_umod_CHG, read1_umod_CpG = mbias.initialise_data_counters(30)
        read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH = mbias.mbias_calculator(99, 'test', [(3, 3), (8, 8), (10, 10), (15, 15), (18, 18), (21, 21), (22, 22), (23, 23), (25, 25), (27, 27), (28, 28)], [], [3, 8, 10, 15, 18, 21, 22, 23, 25, 27, 28], 'AAACAATACGCGgggCTGCATCCCGCGCCG', 'AAACAATACGCGgggCTGCATCCCGCGCCG', 30, read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH, 'mCtoT', False)
        self.assertEqual(read1_mods_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_umod_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 1, 9: 0, 10: 1, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 1, 24: 0, 25: 1, 26: 0, 27: 0, 28: 1, 29: 0})
        self.assertEqual(read1_mods_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_umod_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 1, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 1, 23: 0, 24: 0, 25: 0, 26: 0, 27: 1, 28: 0, 29: 0})
        self.assertEqual(read1_mods_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_umod_CHH, {0: 0, 1: 0, 2: 0, 3: 1, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 1, 19: 0, 20: 0, 21: 1, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        
        
    def test_mbias_calculator_taps_read2(self):
        read2_mods_CHH, read2_mods_CHG, read2_mods_CpG = mbias.initialise_data_counters(30)
        read2_umod_CHH, read2_umod_CHG, read2_umod_CpG = mbias.initialise_data_counters(30) 
        read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH = mbias.mbias_calculator(83, 'test', [(9, 9), (11, 11), (12, 12), (13, 13), (14, 14), (17, 17), (24, 24), (26, 26), (29, 29)], [], [9, 11, 12, 13, 14, 17, 24, 26, 29], 'AAACAATACGCGgggCTGCATCCCGCGCCG', 'AAACAATACGCGgggCTGCATCCCGCGCCG', 30, read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH, 'mCtoT', False)
        self.assertEqual(read2_mods_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_umod_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 1, 10: 0, 11: 1, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 1, 25: 0, 26: 1, 27: 0, 28: 0, 29: 1})
        self.assertEqual(read2_mods_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_umod_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 1, 13: 0, 14: 0, 15: 0, 16: 0, 17: 1, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_mods_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_umod_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 1, 14: 1, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})

    def test_mbias_calculator_bs_read1(self):
        read1_mods_CHH, read1_mods_CHG, read1_mods_CpG = mbias.initialise_data_counters(30)
        read1_umod_CHH, read1_umod_CHG, read1_umod_CpG = mbias.initialise_data_counters(30)
        read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH = mbias.mbias_calculator(147, 'test', [(3, 3), (8, 8), (10, 10), (15, 15), (18, 18), (21, 21), (22, 22), (23, 23), (25, 25), (27, 27), (28, 28)], [3, 8, 10, 15, 18, 21, 22, 23, 25, 27, 28], [], 'AAACAATACGCGgggCTGCATCCCGCGCCG', 'AAACAATACGCGgggCTGCATCCCGCGCCG', 30, read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH, 'CtoT', False)
        self.assertEqual(read1_mods_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 1, 9: 0, 10: 1, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 1, 24: 0, 25: 1, 26: 0, 27: 0, 28: 1, 29: 0})
        self.assertEqual(read1_umod_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_mods_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 1, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 1, 23: 0, 24: 0, 25: 0, 26: 0, 27: 1, 28: 0, 29: 0})
        self.assertEqual(read1_umod_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_mods_CHH, {0: 0, 1: 0, 2: 0, 3: 1, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 1, 19: 0, 20: 0, 21: 1, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_umod_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        
        
    def test_mbias_calculator_bs_read2(self):
        read2_mods_CHH, read2_mods_CHG, read2_mods_CpG = mbias.initialise_data_counters(30)
        read2_umod_CHH, read2_umod_CHG, read2_umod_CpG = mbias.initialise_data_counters(30)
        read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH = mbias.mbias_calculator(163, 'test', [(9, 9), (11, 11), (12, 12), (13, 13), (14, 14), (17, 17), (24, 24), (26, 26), (29, 29)], [9, 11, 12, 13, 14, 17, 24, 26, 29], [], 'AAACAATACGCGgggCTGCATCCCGCGCCG', 'AAACAATACGCGgggCTGCATCCCGCGCCG', 30, read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH, 'CtoT', False)
        self.assertEqual(read2_mods_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 1, 10: 0, 11: 1, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 1, 25: 0, 26: 1, 27: 0, 28: 0, 29: 1})
        self.assertEqual(read2_umod_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_mods_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 1, 13: 0, 14: 0, 15: 0, 16: 0, 17: 1, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_umod_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_mods_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 1, 14: 1, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_umod_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        
    
if __name__ == '__main__':
    unittest.main()
