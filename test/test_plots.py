'''
Test the scoring functions
'''
from __future__ import print_function
import unittest
import data
import numpy
import scipy

class TestSequenceBasedPlots(unittest.TestCase):

    def create_empty_datasets(self, n_sets, sequence="SEAQWEST"):
        datasets = []
        for i in range(nsets):
            datasets.append(data.Dataset("Data"+str(i), data.Conditions(), sequence))
        return datasets


    def test_basic_residue_plot(self):
        datasets = create_empty_datasets(2)

        fig, ax = plt.figure()

        ax = plots.plot()



if __name__ == '__main__':
    unittest.main()