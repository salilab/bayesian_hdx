"""
   Classes that store the representation and the sampled parameters for each state
   Models are defined for each system state.
   They contain the sampled parameters for each state.
"""
from __future__ import print_function
import system
import sampling
import numpy
import scipy
from numpy import linalg
import sys
from copy import deepcopy
import os
import os.path
import time
from random import randint

class ResidueGridModel(object):
    ''' 
    Models the system as individual residues with protection factors along a grid 
    Defined by a grid_size and the parameters of the datasets included in the state

    Calculates the grid of observable protection factors
    and converts between grid values and protection factors.

    @param state - the state that this model applies to
    @param grid_size - the size of the sampling grid
    @param protection_factors - Boolean. Set to true to calculate protection factors. False to calculate rates only.
    '''
    def __init__(self, state, grid_size, protection_factors=False, 
                sample_back_exchange=True, sample_only_observed_residues=True):
        self.state = deepcopy(state)
        self.length = len(state.get_sequence())

        self.grid_size = grid_size
        self.protection_factors = protection_factors
        self.model = numpy.zeros(self.length)   # The model values that are used in the sampler
        self.model_protection_factors = numpy.ones(self.length)*-1  # The protection factors
        self.sampler_type = "int"
        self.sampler_size = range(1, grid_size+1)
        self.state.set_output_model(self)
        self.calculate_protection_factor_grids()
        self.sample_back_exchange = sample_back_exchange
        self.sample_only_observed_residues = sample_only_observed_residues

    def generate_model(self, random=True, value=1, initialize=False):
        # Using the sequence from the state, generate a random model
        # from the sampling grid
        sequence = self.state.get_sequence()
        model = numpy.zeros(self.length)

        if not random and (type(value) is not int or value > self.grid_size or value < 1):
            raise Exception("ResidueGridModel.generate_model: Value error. Either allow random assignment or set an integer value from 1 to grid_size")

        for i in range(self.length):

            if sequence[i] == "P" or (i+1 not in self.state.get_observable_residue_numbers() and self.sample_only_observed_residues):
                model[i] = 0
            else:
                if random:
                    model[i] = numpy.random.randint(1, high=self.grid_size, size=1)[0]
                else:
                    model[i] = int(value)
        if initialize:
            self.model = model

        #if self.protection_factors:
        self.convert_model_to_protection_factors(model)
        #else:
        #    self.convert_model_to_rates(model)
        return model

    def get_sampling_grid(self):
        return range(1, self.grid_size+1)

    def sampler_range(self):
        # Function that returns that available values to each residue
        return self.get_sampling_grid()

    def convert_model_to_rates(self, model):
        return self.convert_model_to_protection_factors(model)

    def get_pf(self, res, i):
        # Given an integer and residue number, return the corresponding grid value that
        # Corresponds to the protection factor.
        #print("PFG", res, i)
        return self.pf_grids[res-1][i-1]

    def convert_model_to_protection_factors(self, model):
        # For a vector of grid values (the model), return a vector of protection values
        pf_grids = self.pf_grids

        for i in range(self.length):
            if int(model[i]-1) == -1:
                self.model_protection_factors[i] = -1
            else:
                self.model_protection_factors[i] = self.get_pf(i, int(model[i]-1))
        
        #self.model_protection_factors = [pf_grids[i][int(model[i]-1)] for i in range(self.length)]
        return self.model_protection_factors

    def calculate_protection_factor_grids(self, threshold = 0.01):
        '''
        We theorize that the protection factor will be constant over multiple experiments.
        Therefore, the grid of observable protection factors for a given residue is dependent
        on all datasets.

        @param threshold - Means that the grid will start at a protection factor value where (1-threshold) * 100%
                    will be observed at the first timepoint and end at a protection factor value where threshold * 100%
                    of exchange will be observed for this site at the longest timepoint.
        '''
        bounds = []

        # For each dataset, find the observable thresholds
        if not self.state.has_data:
            raise Exception("Cannot calculate protection factor grid because state " + self.state.name + " has no associated datasets")
        orb_slow = []
        orb_fast = []
        for d in self.state.get_datasets():
            orb_slow.append(d.calculate_observable_rate_bounds(threshold)[0])
            orb_fast.append(d.calculate_observable_rate_bounds(threshold)[0])

        # Get the fastest and slowest log(rate)
        #min_rate = min(orb_slow)
        #max_rate = max(orb_fast)

        pf_grids = []
        for n in range(self.length):
            #pf_ranges = [pf[n] for pf in observable_pfs]
            if self.protection_factors:
                # Protection factors are set up between 0 and 10
                pf_grid = numpy.linspace( 0, 14, self.grid_size )
            else:
                pf_grid = numpy.linspace( 0, 14, self.grid_size )
            pf_grids.append(pf_grid)

        self.pf_grids = pf_grids

        return self.pf_grids

    def get_conversion_grid(self):
        '''
        Returns a list (of length = # of residues) of lists (of length = grid_size)
        with the protection factor values for each residue.
        '''
        return self.pf_grids
    '''
    def output_conversion_model(self, outfile):
        model_string = ""
        for i in self.model_protection_factors:
            the_string+=" "+str(i)
        the_string = str(self.model)

        # Conversion_string is a list of
        pf_grid_value_string = str(pf_grid)
        return pf_grid_value_string, model_string[0:-1]
    '''
    def change_residue(self, residue, value):
        # Changes the grid value of a single residue in the model in both the grid and protection factor. 
        if type(value) is not int:
            raise Exception("ResidueGridModel.change_residue : value must be of type int. Currently ", type(value))
        if value >=1 and value < self.grid_size +1:
            self.model[residue-1] = value
            self.model_protection_factors[residue-1] = self.pf_grids[residue-1][int(value-1)]
        else:
            raise Exception("ResidueGridModel.change_residue : value must be between 1 and ", grid_size +1, ". Currently ",  value)           
    
    def get_model_residue(self, residue):
        # Returns the integer value of a single residue in the current model
        return self.model[residue-1]

    def get_model(self):
        # Returns the model, a list of integers.  
        return self.model

    def get_masked_model(self, observed_residues):
        # get model only returning those in the observed_residues list.
        # All other residues are = 0
        mod = []
        for i in range(len(self.model)):
            if i+1 in observed_residues:
                mod.append(int(self.model[i]))
            else:
                mod.append(0)
        return mod

    def get_current_model(self):
        return self.convert_model_to_protection_factors(self.model)

