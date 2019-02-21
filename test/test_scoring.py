'''
Test the scoring functions
'''
from __future__ import print_function
import unittest
import numpy
import scipy
import os
import utils

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

import model
import scoring
import data
#import tools



def initialize_system_and_dataset():
    sequence = "ILIKEPLINYTHEELDER"
    sys = system.System()       
    mol = sys.add_molecule(sequence, "test_molecule")
    d = data.Dataset("Test", data.Conditions(), sequence=sequence)
    d.create_peptide("MEGA", start_residue=1)
    d.create_peptide("MEGAMA", start_residue=1)
    d.create_peptide("GAMAN", start_residue=3)
    d.create_peptide("AMAN", start_residue=4)

    for pep in d.get_peptides():
        pep.add_timepoints([10, 100, 1000])

    mol.get_state(0).add_dataset(d)
    return mol.get_state(0)


class TestGaussianNoiseModel(unittest.TestCase):
    def test_functions(self):
        gnm = scoring.GaussianNoiseModel()
        model = 0.45
        exp = 0.45
        sigma = 0.05
        sigma0 = 0.1

        # Test the replicate scores
        for exp in numpy.linspace(0,1.0,10):
            score = numpy.exp(-((model-exp)**2)/(2*sigma**2))/(sigma*numpy.sqrt(2*numpy.pi))
            #print(model, exp, -numpy.log10(score))
            self.assertAlmostEqual(gnm.replicate_score(model, exp, sigma), score)

        exp = 0.65
        for sigma in numpy.linspace(0.02,0.5,25):
            score = numpy.exp(-((model-exp)**2)/(2*sigma**2))/(sigma*numpy.sqrt(2*numpy.pi))
            #print(model, exp, sigma,  -numpy.log10(score))
            self.assertAlmostEqual(gnm.replicate_score(model, exp, sigma), score)

        for sigma in numpy.linspace(0.02,0.5,25):
            likelihood = (20.0 / sigma**2) * numpy.exp(-sigma0**2 / sigma**2)
            print(sigma0, sigma, -1*numpy.log10(likelihood))
            self.assertEqual(likelihood, gnm.experimental_sigma_prior(sigma, sigma0))


class TestTruncatedGaussianNoiseModel(unittest.TestCase):
    def test_functions(self):
        gnm = scoring.GaussianNoiseModel()
        model = 0.45
        exp = 0.45
        sigma = 0.05
        sigma0 = 0.1

        # Test the replicate scores
        for exp in numpy.linspace(0,1.0,10):
            score = numpy.exp(-((model-exp)**2)/(2*sigma**2))/(sigma*numpy.sqrt(2*numpy.pi))
            #print(model, exp, -numpy.log10(score))
            self.assertAlmostEqual(gnm.replicate_score(model, exp, sigma), score)

        exp = 0.65
        for sigma in numpy.linspace(0.02,0.5,25):
            score = numpy.exp(-((model-exp)**2)/(2*sigma**2))/(sigma*numpy.sqrt(2*numpy.pi))
            #print(model, exp, sigma,  -numpy.log10(score))
            self.assertAlmostEqual(gnm.replicate_score(model, exp, sigma), score)

        for sigma in numpy.linspace(0.02,0.5,25):
            likelihood = (20.0 / sigma**2) * numpy.exp(-sigma0**2 / sigma**2)
            print(sigma0, sigma, -1*numpy.log10(likelihood))
            self.assertEqual(likelihood, gnm.experimental_sigma_prior(sigma, sigma0))

    def test_truncation(self):
        model = 0.45
        exp = 0.45
        sigma = 0.05
        sigma0 = 0.1
        ub = 1
        lb = 0

        tfact = 1/ ( 0.5 * ( scipy.special.erf( (ub-exp)/sigma * numpy.sqrt(3.1415) ) - scipy.special.erf( (lb-exp)/sigma * numpy.sqrt(3.1415) ) )  )


        gnm = scoring.GaussianNoiseModel()
        tgnmxx = scoring.GaussianNoiseModel(True)
        tgnm01 = scoring.GaussianNoiseModel(True, (0,1))
        tgnm0 = scoring.GaussianNoiseModel(True, (0,None))
        # Test the replicate scores
        for exp in numpy.linspace(0,1.0,101):
            score = gnm.replicate_score(model, exp, sigma)
            tscore = tgnmxx.replicate_score(model, exp, sigma)
            print(model, exp, tfact, score, tscore)
            self.assertAlmostEqual(score, tscore)       





if __name__ == '__main__':
    unittest.main()
