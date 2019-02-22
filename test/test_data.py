'''
Test features revolving around data handling
'''
from __future__ import print_function
import unittest
import os
import numpy
import utils

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

import system
import data
import tools
import hxio


class TestDataset(unittest.TestCase):

    def initialize_system(self):
        sys = system.System()
        mol = sys.add_macromolecule(name="Temp", sequence = "LALALAND")
        return mol

    def initialize_dataset(self):
        c = data.Conditions()
        d = data.Dataset("Test", c, sequence="LALALAND")
        return d

    def test_conditions(self):
        c = data.Conditions()
        self.assertEqual(c.saturation, 1.0)

    def test_dataset(self):
        d = self.initialize_dataset()
        d.create_peptide("LALAND", start_residue=3)

        mol = self.initialize_system()
        s = mol.get_state(name="Apo")

        s.add_dataset(d)
        self.assertEqual(s, d.get_state())

        rates = d.calculate_intrinsic_rates()

        std_rates = [numpy.nan, 2.03190342, 0.04565069, 0.41190342, 0.04565069, 0.41190342, 1.11190342, -0.75672096]

        for i in range(len(rates)):
            if numpy.isfinite(rates[i]):
                self.assertAlmostEqual(rates[i], std_rates[i], places=4)

    def test_peptides(self):

        d = self.initialize_dataset()

        p1 = d.create_peptide("LALA", start_residue=1)
        self.assertEqual(len(d.get_peptides()), 1)
        self.assertEqual(p1.get_dataset(), d)

        p2 = d.create_peptide("LALA", start_residue=3)
        self.assertEqual(len(d.get_peptides()), 2)

        p3 = data.Peptide(d, "LALA", 3, None)
        d.add_peptide(p3)
        self.assertEqual(len(d.get_peptides()), 3)

    def test_calculate_protection_factors(self):
        d = self.initialize_dataset()
        p = d.create_peptide("LALA", start_residue=1)
        p.add_timepoints([10, 100, 1000])

        threshold = 0.01
        bounds = d.calculate_observable_rate_bounds(threshold=threshold)
        set_bounds = (numpy.log10(-numpy.log(1-threshold)/1000), numpy.log10(-numpy.log(threshold)/10))
        self.assertAlmostEqual(bounds, set_bounds)

        threshold = 0.05
        bounds = d.calculate_observable_rate_bounds(threshold=threshold)
        set_bounds = (numpy.log10(-numpy.log(1-threshold)/1000), numpy.log10(-numpy.log(threshold)/10))
        self.assertAlmostEqual(bounds, set_bounds)

        obs_pfs = d.calculate_observable_protection_factors()

        self.assertEqual(obs_pfs[0][0], numpy.inf)
        self.assertAlmostEqual(obs_pfs[1][0], 2.5554004223717413)
        self.assertAlmostEqual(obs_pfs[3][1], 4.7018428308264344)

    def test_write_file(self):
        pass




class TestPeptides(unittest.TestCase):


    def initialize(self):
        c = data.Conditions()
        d = data.Dataset("Test", c, sequence="EWESEEESSEFF")
        return d


    def test_adding_timepoints(self):

        d = self.initialize()
        p = d.create_peptide("ESEEE", start_residue=3)
        p.add_timepoints([10, 100, 1000])

        self.assertEqual(len(p.get_timepoints()), 3)


        # Test setting all tp in this peptide to the same sigma
        p.set_timepoint_sigmas(1.5)
        for t in p.get_timepoints():
            self.assertEqual(1.5, t.get_sigma())









if __name__ == '__main__':
    unittest.main()
