'''
Test features of the HDX system representation
'''
from __future__ import print_function
import system
import data
import tools
import hxio
import model
import unittest
import os
import numpy


def initialize_system():
    sys = system.System()       
    return sys.add_macromolecule("MEGAMAN", "test_molecule")

def initialize_dataset():
    d = data.Dataset("Test", data.Conditions(), sequence="MEGAMAN")
    d.create_peptide("MEGA", start_residue=1)
    d.create_peptide("MEGAMA", start_residue=1)
    d.create_peptide("GAMAN", start_residue=3)
    d.create_peptide("AMAN", start_residue=4)
    return d

class TestSystem(unittest.TestCase):

    def test_initialize(self):
        sys = system.System()
        
        self.assertEqual(len(sys.get_macromolecules()), 0)

        sys.add_macromolecule("MEGAMAN", "test_molecule")
        self.assertEqual(len(sys.get_macromolecules()), 1)


class TestMacromolecule(unittest.TestCase):

    def test_macromolecule_functions(self):

        mol = initialize_system()

        self.assertEqual(mol.get_sequence(), "MEGAMAN")

        self.assertEqual(len(mol.get_states()), 1)
        mol.add_state("test", None)


class TestState(unittest.TestCase):

    def test_state_functions(self):

        mol = initialize_system()

        d = initialize_dataset()
        state = mol.get_state(0)

        self.assertFalse(state.has_data())

        state.add_dataset(d)

        self.assertTrue(state.has_data())

        # test calculating sectors
        sectors = state.calculate_sectors()
        self.assertEqual(len(sectors), 4)
        self.assertEqual(state.sector_dictionary[0], set([3,4]))
        self.assertEqual(len(state.sector_dictionary), 4)

        pep = d.create_peptide("EGAM", start_residue=2)
        sectors = state.calculate_sectors()
        self.assertEqual(len(sectors), 5)       
        self.assertEqual(state.sector_dictionary[0], set([3]))

        pep.add_timepoints([10,100])

        state.set_output_model(model.ResidueGridModel(state, 10))
        state.initialize()









if __name__ == '__main__':
    unittest.main()