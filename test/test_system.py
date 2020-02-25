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
    times=[10,30,90,900,3600]
    d = data.Dataset("Test", data.Conditions(), sequence="MEGAMAN")
    d.create_peptide("MEGA", start_residue=1)
    d.create_peptide("MEGAMA", start_residue=1)
    d.create_peptide("GAMAN", start_residue=3)
    d.create_peptide("AMAN", start_residue=4)
    for p in d.get_peptides():
        for t in times:
            p.add_timepoint(t)
    return d

class TestSystem(unittest.TestCase):

    def setup_system_with_data(self):
        mol = initialize_system()
        state = mol.get_apo_state()
        d = initialize_dataset()
        state.add_dataset(d)
        state.set_output_model(model.ResidueGridModel(state, 50, protection_factors=True))

        return mol, state

    def test_initialize(self):
        sys = system.System()

        self.assertEqual(len(sys.get_macromolecules()), 0)

        sys.add_macromolecule("MEGAMAN", "test_molecule")
        self.assertEqual(len(sys.get_macromolecules()), 1)

    def test_residue_information_content(self):
        mol, state = self.setup_system_with_data()
        ri = state.calculate_residue_information_content()
        self.assertEqual(len(state.sequence), len(ri))
        #print("RI:", ri)



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
