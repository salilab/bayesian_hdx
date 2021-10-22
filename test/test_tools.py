'''
Test features revolving around data handling
'''
from __future__ import print_function
import unittest
import os
import utils

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

import tools


class TestDataset(unittest.TestCase):

    def test_get_neighbor_effects(self):
        # Test the intrinsic neighbor effects
        # Make sure that the pH and temp dependence formulas work
        # All values are compared to those calculated from the
        # Englander EXCEL worksheet: FBMME_HD.XLS (20-Jan-10 update)

        ne = tools.get_residue_neighbor_effects("D", 8.5, 283)
        std_ne = (0.89996047, 0.57996785, 0.10011608, -0.17979819)

        for i in range(len(ne)):
            self.assertAlmostEqual(ne[i], std_ne[i], places=3)

        ne = tools.get_residue_neighbor_effects("CT", 8.5, 283)
        std_ne = (0.95990099, "", -1.8, "")
        for i in range(len(ne)):
            self.assertAlmostEqual(ne[i], std_ne[i], places=3)

        ne = tools.get_residue_neighbor_effects("E", 8.5, 283)
        std_ne = (-0.89988767, 0.30991680, -0.10986019, -0.14972157)
        for i in range(len(ne)):
            self.assertAlmostEqual(ne[i], std_ne[i], places=3)

        ne = tools.get_residue_neighbor_effects("H", 8.5, 283)
        std_ne = (-0.02304667, -0.01883294, 0.05425316, 0.23320880)
        for i in range(len(ne)):
            self.assertAlmostEqual(ne[i], std_ne[i], places=3)

        ne = tools.get_residue_neighbor_effects("H", 6.5, 283)
        std_ne = (-0.56856927, -0.39726092, 0.74651445, 0.78158214)
        for i in range(len(ne)):
            self.assertAlmostEqual(ne[i], std_ne[i], places=3)

        ne = tools.get_residue_neighbor_effects("H", 7.4, 310)
        std_ne = (-0.08222552, -0.06632590, 0.28438767, 0.39501777)
        for i in range(len(ne)):
            self.assertAlmostEqual(ne[i], std_ne[i], places=3)

    def test_intrinsic_calculator(self):
        # All intrinsic rate values are compared to values transcribed from the
        # Englander EXCEL worksheet: FBMME_HD.XLS (20-Jan-10 update)

        seq = "LALALAND"
        rates = tools.get_sequence_intrinsic_rates(seq, pH=7.0, T=293)
        std_rates = [0., 107.6226, 1.1108, 2.5817, 1.1108, 2.5817,
                     12.9391, 0.1751]

        for i in range(len(rates)):
            self.assertAlmostEqual(rates[i], std_rates[i], places=3)

        rates = tools.get_sequence_intrinsic_rates(seq, pH=7.0, T=274)
        std_rates = [0., 14.2077, 0.1466, 0.3408, 0.1466, 0.3408, 1.7081,
                     0.0231]

        for i in range(len(rates)):
            self.assertAlmostEqual(rates[i], std_rates[i], places=3)

        seq = "PINKPANTHER"
        rates = tools.get_sequence_intrinsic_rates(seq, pH=6.0, T=283)
        std_rates = [0., 0.6672, 0.2718, 0.2846, 0, 0.0859, 0.4616, 0.2679,
                     1.3522, 0.7444, 0.0022]
        for i in range(len(rates)):
            self.assertAlmostEqual(rates[i], std_rates[i], places=3)

    def test_calculate_isotopic_distribution(self):
        seq = "AA"
        seq = "DAREDEVIL"
        seq = "CAPTAINAMERICA"

    def test_calculate_deut(self):
        # Probably the most used function.  Better make sure it's right!!
        rate = -1
        time = 10
        self.assertAlmostEqual(
            tools.calculate_deut(rate, time), 0.63212, places=4)

    def test_subsequence_consistency(self):
        sequence = "SANDMAN"
        self.assertTrue(tools.subsequence_consistency(
            sequence, sequence[2:5], 3))

        self.assertFalse(tools.subsequence_consistency(
            sequence, sequence[2:5], 2))


if __name__ == '__main__':
    unittest.main()
