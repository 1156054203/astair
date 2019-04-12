import sys
import unittest

from astair import mbias as mbias

class MbiasContextSearchTest(unittest.TestCase):
    """Tests whether given the read flag and the method, the correct bases for modification will be found."""

    def test_flags_modified_taps_flag_99(self):
        """Tests whether cytosines converted into thymines will be recognised as modified given that the method 
        is mCtoT and the read flag is 99."""
        self.assertEqual(mbias.strand_and_method(99, 'ACTGCTCCCTGGaaaTCG', 'ATTGCTCCTTGGaaaTCG', 'mCtoT', False), [8, 1])
        
        
    def test_flags_modified_taps_flag_147(self):
        """Tests whether cytosines converted into thymines will be recognised as modified given that the method 
        is mCtoT and the read flag is 147."""
        self.assertEqual(mbias.strand_and_method(147, 'ACTGCTCCCTGGaaaTCG', 'ATTGCTCCTTGGaaaTCG', 'mCtoT', False), [8, 1])


    def test_flags_modified_taps_flag_83(self):
        """Tests whether guanines converted into adeninies will be recognised as modified given that the method 
        is mCtoT and the read flag is 83."""
        self.assertEqual(mbias.strand_and_method(83, 'ACTGCTCCCTGGaaaTCG', 'ATTACTCCTTGAaaaTCG', 'mCtoT', False), [3, 11])
        
        
    def test_flags_modified_taps_flag_163(self):
        """Tests whether guanines converted into adeninies will be recognised as modified given that the method 
        is mCtoT and the read flag is 163."""
        self.assertEqual(mbias.strand_and_method(163, 'ACTGCTCCCTGGaaaTCG', 'ATTGCTCCTTAAaaaTCG', 'mCtoT', False), [10, 11])
        

    def test_flags_modified_bs_flag_99(self):
        """Tests whether cytosines that were not converted into thymines will be recognised as modified given that
        the method is CtoT and the read flag is 99."""
        self.assertEqual(mbias.strand_and_method(99, 'ACTGCTCCCTGGaaaTCG', 'ATTGCTCCTTGGaaaTCG', 'CtoT', False), [16, 4, 6, 7])
        
        
    def test_flags_modified_bs_flag_147(self):
        """Tests whether cytosines that were not converted into thymines will be recognised as modified given that
        the method is CtoT and the read flag is 147."""
        self.assertEqual(mbias.strand_and_method(147, 'ACTGCTCCCTGGaaaTCG', 'ATTGCTCCTTGGaaaTCG', 'CtoT', False), [16, 4, 6, 7])


    def test_flags_modified_bs_flag_83(self):
        """Tests whether guanines that were not converted into adenines will be recognised as modified given that
        the method is CtoT and the read flag is 83."""
        self.assertEqual(mbias.strand_and_method(83, 'ACTGCTCCCTGGaaaTCG', 'ATTACTCCTTGAaaaTCG', 'CtoT', False), [17, 10])
        
        
    def test_flags_modified_bs_flag_163(self):
        """Tests whether guanines that were not converted into adenines will be recognised as modified given that
        the method is CtoT and the read flag is 83."""
        self.assertEqual(mbias.strand_and_method(163, 'ACTGCTCCCTGGaaaTCG', 'ATTGCTCCTTAAaaaTCG', 'CtoT', False), [17, 3])
        
        
    def test_mbias_calculator_taps_read1(self):
        read1_mods_CHH, read1_mods_CHG, read1_mods_CpG = mbias.initialise_data_counters(30)
        read1_umod_CHH, read1_umod_CHG, read1_umod_CpG = mbias.initialise_data_counters(30)
        read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH = mbias.mbias_calculator(99, 'AAACAATACGCGgggCTGCATCCCGCGCCG', 'AAACAATACGCGgggCTGCATCCCGCGCCG', 30, read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH, 'mCtoT', False)
        self.assertEqual(read1_mods_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_umod_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 1, 9: 0, 10: 1, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 1, 24: 0, 25: 1, 26: 0, 27: 0, 28: 1, 29: 0})
        self.assertEqual(read1_mods_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_umod_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 1, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 1, 23: 0, 24: 0, 25: 0, 26: 0, 27: 1, 28: 0, 29: 0})
        self.assertEqual(read1_mods_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_umod_CHH, {0: 0, 1: 0, 2: 0, 3: 1, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 1, 19: 0, 20: 0, 21: 1, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        
        
    def test_mbias_calculator_taps_read2(self):
        read2_mods_CHH, read2_mods_CHG, read2_mods_CpG = mbias.initialise_data_counters(30)
        read2_umod_CHH, read2_umod_CHG, read2_umod_CpG = mbias.initialise_data_counters(30)
        read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH = mbias.mbias_calculator(163, 'AAACAATACGCGgggCTGCATCCCGCGCCG', 'AAACAATACGCGgggCTGCATCCCGCGCCG', 30, read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH, 'mCtoT', False)
        self.assertEqual(read2_mods_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_umod_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 1, 10: 0, 11: 1, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 1, 25: 0, 26: 1, 27: 0, 28: 0, 29: 1})
        self.assertEqual(read2_mods_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_umod_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 1, 13: 0, 14: 0, 15: 0, 16: 0, 17: 1, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_mods_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_umod_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 1, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})

    def test_mbias_calculator_bs_read1(self):
        read1_mods_CHH, read1_mods_CHG, read1_mods_CpG = mbias.initialise_data_counters(30)
        read1_umod_CHH, read1_umod_CHG, read1_umod_CpG = mbias.initialise_data_counters(30)
        read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH = mbias.mbias_calculator(99, 'AAACAATACGCGgggCTGCATCCCGCGCCG', 'AAACAATACGCGgggCTGCATCCCGCGCCG', 30, read1_mods_CpG, read1_mods_CHG, read1_mods_CHH, read1_umod_CpG, read1_umod_CHG, read1_umod_CHH, 'CtoT', False)
        self.assertEqual(read1_mods_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 1, 9: 0, 10: 1, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 1, 24: 0, 25: 1, 26: 0, 27: 0, 28: 1, 29: 0})
        self.assertEqual(read1_umod_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_mods_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 1, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 1, 23: 0, 24: 0, 25: 0, 26: 0, 27: 1, 28: 0, 29: 0})
        self.assertEqual(read1_umod_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_mods_CHH, {0: 0, 1: 0, 2: 0, 3: 1, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 1, 19: 0, 20: 0, 21: 1, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read1_umod_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        
        
    def test_mbias_calculator_bs_read2(self):
        read2_mods_CHH, read2_mods_CHG, read2_mods_CpG = mbias.initialise_data_counters(30)
        read2_umod_CHH, read2_umod_CHG, read2_umod_CpG = mbias.initialise_data_counters(30)
        read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH = mbias.mbias_calculator(163, 'AAACAATACGCGgggCTGCATCCCGCGCCG', 'AAACAATACGCGgggCTGCATCCCGCGCCG', 30, read2_mods_CpG, read2_mods_CHG, read2_mods_CHH, read2_umod_CpG, read2_umod_CHG, read2_umod_CHH, 'CtoT', False)
        self.assertEqual(read2_mods_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 1, 10: 0, 11: 1, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 1, 25: 0, 26: 1, 27: 0, 28: 0, 29: 1})
        self.assertEqual(read2_umod_CpG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_mods_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 1, 13: 0, 14: 0, 15: 0, 16: 0, 17: 1, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_umod_CHG, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_mods_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 1, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        self.assertEqual(read2_umod_CHH, {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0})
        
    
if __name__ == '__main__':
    unittest.main()
