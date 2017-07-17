'''
Test the input of various file types
'''
from __future__ import print_function
import model
import system
import unittest
import data
import numpy
from copy import deepcopy
#import tools


def initialize_system_and_dataset():
    sequence = "MEGAMAN"
    sys = system.System()       
    mol = sys.add_macromolecule(sequence, "test_molecule")
    d = data.Dataset("Test", data.Conditions(), sequence=sequence)
    d.create_peptide("MEGA", start_residue=1)
    d.create_peptide("MEGAMA", start_residue=1)
    d.create_peptide("GAMAN", start_residue=3)
    d.create_peptide("AMAN", start_residue=4)

    for pep in d.get_peptides():
        pep.add_timepoints([10, 100, 1000])

    mol.get_state(0).add_dataset(d)
    return mol.get_state(0)


class TestResidueGridModel(unittest.TestCase):

    def test_change_residue(self):

        state = initialize_system_and_dataset()
        rgm = model.ResidueGridModel(state, grid_size=10)

        mod = rgm.generate_model(False, value=1, initialize=True)
        pfs = rgm.model_protection_factors

        dc_model = deepcopy(rgm.model)
        dc_model_pf = deepcopy(rgm.model_protection_factors)

        self.assertEqual(mod[1], rgm.model[1])
        self.assertEqual(pfs[1], rgm.model_protection_factors[1])

        rgm.change_residue(5, 4)

        print("MODEL", rgm.get_model_residue(5), rgm.model[4:6], rgm.model_protection_factors[4:6])
        print("REG_COPY", mod[4:6], pfs[4:6])
        print("DEEP_COPY", dc_model[4:6], dc_model_pf[4:6])


    def test_residue_grid_model(self):

        state = initialize_system_and_dataset()
        rgm = model.ResidueGridModel(state, grid_size=10)

        pf_grids = rgm.calculate_protection_factor_grids(threshold = 0.01)

        mod = rgm.generate_model(random=False, value=2, initialize=True)
        self.assertEqual(len(state.get_sequence()), len(mod))
        for i in mod:
            self.assertEqual(i, 2)

        mod2 = rgm.generate_model(random=True, initialize=True)
        self.assertEqual(len(mod2), len(mod))
        for i in range(len(mod2)):
            self.assertEqual(mod2[i], rgm.model[i])

        #print(rgm.grid_size, rgm.get_sampling_grid(), pf_grids[0])





if __name__ == '__main__':
    unittest.main()