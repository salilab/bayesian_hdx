"""
   Scoring functions for HDX models
"""
from __future__ import print_function
import numpy
import scipy
import math


class GaussianNoiseModel(object):
    '''
    Takes a data.Dataset object, a set of protection factors and priors
    and calculates a score.
    '''
    def __init__(self, truncated=False, bounds=(None, None)):
        self.truncated = truncated

        if truncated:
            # 10 and -10 are numerically equivalent to no truncation
            # factor when data is in the range of 0-->1
            if bounds[1] is None:
                self.upper_bound = 10
            else:
                self.upper_bound = bounds[1]

            if bounds[0] is None:
                self.lower_bound = -10
            else:
                self.lower_bound = bounds[0]

    def replicate_score(self, model, exp, sigma):
        # Forward model
        raw_likelihood = math.exp(-((model-exp)**2)/(2*sigma**2)) \
            / (sigma*math.sqrt(2*numpy.pi))
        if self.truncated:
            raw_likelihood *= 1 / (0.5 * (
                scipy.special.erf((self.upper_bound-exp)/sigma
                                  * math.sqrt(3.1415))
                - scipy.special.erf((self.lower_bound-exp)/sigma
                                    * math.sqrt(3.1415))))
        return raw_likelihood

    def deuterium_incorporation_prior(self, exp):
        # Prior on the experimental observation
        # For now, choosing a uniform prior over all values of 2D
        # incorporation.
        return 1.0

    def protection_factor_prior(self, protection_factors):
        # Prior on the protection factor at each residue in the peptide
        # Currently, just using a uniform prior
        prior = 1.0
        for p in protection_factors:
            prior *= 1.0
        self.pf_prior = prior
        return prior

    def peptide_confidence_score(self, peptide):
        # User-definable function for converting peptide confidence into
        # a likelihood.
        # The function must be evaluatable between 0 and 1.
        pass

    def experimental_sigma_prior(self, sigma, sigma0):
        # Prior on the sigma value. Long tailed to allow for outlier values.
        return (20.0 / sigma**2) * math.exp(-sigma0**2 / sigma**2)
