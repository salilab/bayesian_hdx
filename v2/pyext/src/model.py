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
    def __init__(self, state, grid_size, protection_factors=False, sample_back_exchange=True):
        self.state = deepcopy(state)
        self.length = len(state.get_sequence())

        self.grid_size = grid_size
        self.protection_factors = protection_factors
        self.model_protection_factors = numpy.zeros(self.length)
        self.sampler_type = "int"
        self.sampler_size = range(1, grid_size+1)
        self.state.set_output_model(self)
        self.calculate_protection_factor_grids()
        self.sample_back_exchange = sample_back_exchange

    def generate_model(self, random=True, value=1, initialize=False):
        # Using the sequence from the state, generate a random model
        # from the sampling grid
        sequence = self.state.get_sequence()
        model = numpy.zeros(self.length)

        if not random and (type(value) is not int or value > self.grid_size or value < 1):
            raise Exception("ResidueGridModel.generate_model: Value error. Either allow random assignment or set an integer value from 1 to grid_size")

        for i in range(self.length):
            if sequence[i] == "P" or i+1 not in self.state.get_observable_residue_numbers():
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

    def convert_model_to_protection_factors(self, model):
        # For a vector of grid values (the model), return a vector of protection values
        pf_grids = self.pf_grids

        for i in range(self.length):
            if int(model[i]-1) == -1:
                self.model_protection_factors[i] = -1
            else:
                self.model_protection_factors[i] = pf_grids[i][int(model[i]-1)]
        
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
        min_rate = min(orb_slow)
        max_rate = max(orb_fast)

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

    def output_conversion_model(self, outfile):
        model_string = ""
        for i in self.model:
            the_string+=" "+str(i)
        the_string = str(self.model)

        # Conversion_string is a list of
        pf_grid_value_string = str(pf_grid)
        return pf_grid_value_string, model_string[0:-1]

    def change_residue(self, residue, value):
        # Changes the grid value of a single residue in the model
        if type(value) is not int:
            raise Exception("ResidueGridModel.change_residue : value must be of type int. Currently ", + type(int))
        if value >=1 and value < self.grid_size +1:
            self.model[residue-1] = value
            self.model_protection_factors[residue-1] = self.pf_grids[residue-1][int(value-1)]
        else:
            raise Exception("ResidueGridModel.change_residue : value must be between 1 and ", grid_size +1, ". Currently ", + int)           
    
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





class MultiExponentialModel(object):
    """ A multi-exponential model for exchange that approximates each independent sector
        as a sum of exponentials, corresponding to the number of exchanging amides in that sector.
        The exchange rate constant, kex,  for each residue is sampled along a grid of log(kex) values.

        The model can be initialized to a certain bin in the exponential grid, or to random values.
        The Bayesian likelihood of the model is calculated against a set of hdx.system_setup.Fragment objects
        which contain experimental data.
    """

    def __init__(self, model, state, frags=None, sigma=10.0, num_exp_bins=10, init="enumerate", 
                 error_model="gaussian", output_directory=None, noclobber=True, first_tp=10, last_tp=3600, marginalize=False):
        if frags is None:
            self.frags=state.frags
        else:
            self.frags=frags
            #print "using frags"
        tpoints=[]
        times=[]

        # get first and last time points
        for f in state.frags:
            tpoints += [tp for tp in f.timepoints]
            times+=[t.time for t in tpoints]
        first_tp = numpy.min(times)
        last_tp = numpy.max(times)

        self.first_tp = first_tp
        self.last_tp = last_tp

        self.state=state
        self.state.exp_model=self
        self.exp_grid=self.calc_exp_grid(num_exp_bins, first_tp, last_tp)
        self.exp_seq=[0]*len(model.seq)
        self.name=state.state_name
        state.add_exp_model(self)
        self.offset=model.offset
        self.sigma=sigma
        self.oldsigma=-1
        self.scores=[]
        self.seq=model.seq
        self.flag_outliers=True
        self.output_directory=output_directory
        self.noclobber=noclobber
        self.marginalize=marginalize

        # These grids contain the values at which the integral of sigma is 
        # computed 
        self.sigma0_grid=numpy.linspace(0.1,10,100) #stdev between 0.01 and 0.5
        self.omega_grid=numpy.linspace(1,100,100) #stdev between 0.01 and 0.5
        self.fmod_grid=numpy.linspace(0,100,101) #fmod between 0.00 and 1.00
        self.fexp_grid=self.fmod_grid #fmod between 0.00 and 1.00

        #Initialize each log(k) value to one of the grid points.
        #If there is no coverage, or there is a proline, make it zero
        #so we do not use this position in calculations.  Otherwise, use
        #a random integer to refer to the exp_grid value to use for log(k).
        #Or use a user defined value or just the first value.

        if init=="enumerate":
            self.guess_init_exp_sequence(5)
        else:
            for n in range(0,model.num_res-1):
                if model.seq[n]=='P' or model.seq[n]=='p' or state.coverage[n]==0:
                    self.exp_seq[n]=int(0)
                elif init=="random":
                    self.exp_seq[n]=numpy.random.randint(len(self.exp_grid)-1)+1
                elif isinstance(init,int) and init > 0 and init < len(self.exp_grid)-1:
                    self.exp_seq[n]=init
                elif isinstance(init,int) and len(init)==len(seq):
                    self.exp_seq=init
                else:
                    self.exp_seq[n]=1

        state.exp_model=self
        self.exp_models=[]

        if marginalize:
            self.calculate_marginal_probabilities()

        self.score=self.calculate_bayesian_score(self.frags, sigma, error_model=error_model)

    def set_sigma(self,sigma):
        self.sigma=sigma

    def calculate_bayesian_score(self, frags, sig=None, error_model="gaussian", report_stats = False):
        total_score=0
        joint_likelihood=0
        if sig is None:
            sig = self.sigma
        # Get the sigma for this timepoint

        j = self.get_nearest_value(sig, self.sigma0_grid)
        #print("SIG:: CBS:", sig, j, self.sigma0_grid[j])

        for f in frags:

            frag_likelihood = 0
            frag_replicates = 0

            for tp in f.timepoints:

                tp.clear_model_deuteration_values()
                tp_avg, tp_sd=tp.get_avg_sd()
                tp_likelihood=0

                model_deut=0

                sigma=tp.sigma

                # Calculate model deuteration values
                # for r in range(f.start_res-self.offset,f.end_res-self.offset-1):
                for r in range(f.start_res+self.offset+2,f.end_res+self.offset+1):
                    if self.seq[r] != 'P' and self.seq[r] != 'p':
                        kr=10**self.exp_grid[self.exp_seq[r]-1]
                        model_deut+=(1-numpy.exp(-kr*tp.time))

                model_deut=100.*model_deut/f.num_observable_amides
                tp.add_model_deuteration_value(model_deut)

                #print(f.seq, tp.time, sigma, tp_avg, tp_sd, kr, model_deut, self.exp_grid)

                #fm=self.get_nearest_value(model_deut, self.fmod_grid)

                # Calculate (or lookup) replicate likelihood
                for n in range(len(tp.replicates)):
                    frag_replicates += 1
                    if self.marginalize:

                        #print("Marginalized step")

                        #delfmod=tp.replicates[n].deut-model_deut
                        #sig=100.*tp.get_sigma()/f.num_observable_amides

                        #fe=self.get_nearest_value(tp.replicates[n].deut, self.fexp_grid)   
                        fe = tp.replicates[n].get_nearest_gridpoint(self.fexp_grid) 

                        replicate_score=self.prob_grid[j,fm,fe]
                        
                    else:
                        if error_model=="gaussian":
                            replicate_score = self.gaussian_model(tp.replicates[n].deut, model_deut, sigma)
                        elif error_model=="truncated_gaussian":
                            replicate_score = self.gaussian_model(tp.replicates[n].deut, model_deut, sigma)
                            replicate_score += self.truncated_gaussian_factor(tp.replicates[n].deut,-10,120,sigma)
                        elif error_model=="lognormal":
                            replicate_score = self.lognormal_model(tp.replicates[n].deut, model_deut, sigma)
                        else:
                            print("Unrecognized Error Model:", error_model)
                            sys.exit(1)

                    # Hack to deal with log-based overflows from likelihood=0
                    if replicate_score != 0:
                        tp_likelihood = tp_likelihood + -1*numpy.log(replicate_score) #numpy.log(replicate_score) #+ -1.0*numpy.log(self.get_sigma_score(self.sigma,tp.sigma))
                    else:
                        tp_likelihood = tp_likelihood + 10000




                #print(f.seq, tp.time, model_deut, [r.deut for r in tp.replicates], sigma, replicate_score, tp_likelihood)

                #avg_diff=(model_deut-numpy.average([r.deut for r in tp.replicates]))
            frag_likelihood=frag_likelihood+tp_likelihood
            f.likelihood = frag_likelihood

            joint_likelihood = joint_likelihood + frag_likelihood

            #print(f.seq, f.start_res, frag_likelihood, joint_likelihood)
        return joint_likelihood

    def guess_init_exp_sequence(self, size_enum):
        ''' Guess initial state by finding most likely exp values for each fragment
            Returns a sequence of exp values
        '''
        frags = deepcopy(self.frags)

        exp_seq=[0]*len(self.exp_seq)
        exp_grid=self.calc_exp_grid(size_enum, self.first_tp, self.last_tp)
        print("\n##### Initial Guess: Enumerate #####")
        print("# Enumerating all fragments with gridsize of", size_enum)
        print("# State:", self.state.state_name, "\n")
        print("> Fragment start_res end_res exp_grid")
        for f in frags:
            f.exp_vals = sampling.enumerate_fragment(f, exp_grid, 2.0, num_models = 1)
            print(f.seq, f.exp_vals)
        print("Calculating overlaps and casting to gridsize", len(self.exp_grid))  
        sectors = self.state.get_sectors(frags)
        # sort sectors by size
        sectors.sort(key=lambda x: x.num_amides, reverse=True)

        for s in sectors:
            sector_exp_vals = numpy.zeros(size_enum)
            sector_frags = deepcopy(s.fragments)
            # if there is only one fragment, just place the residues.  We are done.
            if len(sector_frags) == 1 and sum(sector_frags[0].exp_vals) == s.num_amides:
                method = "Single fragment"
                sector_exp_vals = sector_frags[0].exp_vals
                sector_frags[0].exp_vals = deepcopy(sector_frags[0].exp_vals) - sector_exp_vals
            else:
                # Populate sector with overlapping exp values
                # Is total overlap greater or less than s.num_amides?
                # Make an array of all sector fragment exp_grid vals
                # Overlap is the minimum value for each field
                # How to deal with the rest?  High to low or low to high?
                # Three possibilities:   
                #    * Equal - assign overlap to each amide.  Location doesn't matter 
                #    * sum(overlap) > len(seq) - take (random?) subset of overlap and assign to sequence
                #    * sum(overlap) < len(seq) - Assign overlap to sequence. Find closest non-overlapping values

                all_exp_vals = [] #numpy.zeros((len(frags),len(self.exp_grid)))              
                for f in sector_frags:
                    all_exp_vals.append(f.exp_vals)

                overlap = numpy.amin(numpy.array(all_exp_vals), axis=0)
                #------------------
                # For equal number of overlap values and amides:
                #------------------
                if numpy.sum(overlap)==s.num_amides:
                    method = "Equal Overlap"
                    sector_exp_vals = overlap
                    #print("Overlap = num_amides", s.seq, overlap, s.start_res)
                    for f in sector_frags:
                        f.exp_vals = deepcopy(f.exp_vals) - overlap
            
                    #print(exp_seq[s.start_res:s.end_res+1])
                #------------------
                # For more overlap values than amides:
                #   - randomly pick amides to assign
                #------------------
                elif numpy.sum(overlap) > s.num_amides:
                    method = "Overlap > num_amides"
                    overlap_sub = numpy.zeros(len(overlap)).astype(int)
                    i=1
                    while i <= s.num_amides:
                        picked_val = randint(0,len(overlap)-1)
                        #print(i, picked_val, overlap, overlap[picked_val], overlap_sub)
                        if overlap[picked_val] != 0:
                            overlap_sub[picked_val] += 1
                            overlap[picked_val] -= 1
                            i+=1
                    #print("Overlap > num_amides", s.seq, overlap, overlap_sub, s.start_res)
                    
                    sector_exp_vals = overlap_sub

                    for f in sector_frags:
                        f.exp_vals = deepcopy(f.exp_vals) - overlap_sub   

                #------------------
                # For fewer overlap values than amides:
                #   - assign overlap
                #   - start from smallest fragment and assign values until you hit n_remaining
                #------------------ 
                else:  
                    method = "Overlap < num_amides"
                    overlap_sub = numpy.zeros(len(overlap)).astype(int)
                    n_remaining = s.num_amides - numpy.sum(overlap) 
                    for f in sector_frags:
                        f.exp_vals = deepcopy(f.exp_vals) - overlap

                    #sort frags from smallest to largest               
                    sector_frags.sort(key=lambda x: x.num_observable_amides) 
                    #print("***Overlap < num_amides", s.seq, s.num_amides - sum(overlap), overlap, sector_frags[0].exp_vals, s.start_res)
                    i=1
                    #randomly pick value from smallest fragment
                    while i <= n_remaining:               
                        picked_val = randint(0,len(overlap)-1)
                        if sector_frags[0].exp_vals[picked_val] != 0:
                            overlap_sub[picked_val] += 1
                            sector_frags[0].exp_vals[picked_val] = deepcopy(sector_frags[0].exp_vals[picked_val])-1
                            #print(i, n_remaining, picked_val, overlap_sub, sector_frags[0].exp_vals)
                            i = i+1
                    sector_exp_vals = overlap + overlap_sub
                    
                    # For each other frag, delete exp_vals closest to overlap_sub
                    for f in sector_frags[1:-1]:
                        nonzero_elements=f.exp_vals.nonzero()[0]
                        for i in range(len(exp_grid)):
                            if overlap_sub[i] != 0:
                                if i==0:  # if we are at i=0, delete first nonzero
                                    first_nze = nonzero_elements[0]
                                    f.exp_vals[first_nze] = deepcopy(f.exp_vals[first_nze]) - 1

                                elif i==len(exp_grid)-1:  # if we are at i=N, delete last nonzero
                                    last_nze = nonzero_elements[-1]
                                    f.exp_vals[last_nze] = deepcopy(f.exp_vals[last_nze]) - 1

                                else:  #Otherwise a bit more complicated
                                    direction = int(2*(randint(0,1)) - 1) #decide whether to go + or -
                                    #print(i, direction, f.exp_vals)
                                    if f.exp_vals[i+direction] != 0: # try 1st direction
                                        f.exp_vals[i+direction] = deepcopy(f.exp_vals[i+direction]) - 1
                                    elif f.exp_vals[i-direction] != 0: #if not, try second direction
                                        f.exp_vals[i-direction] = deepcopy(f.exp_vals[i-direction]) - 1
                                    else: #cop out for this...just use the first one.
                                        f.exp_vals[nonzero_elements[0]] = deepcopy(f.exp_vals[nonzero_elements[0]]) - 1

            if numpy.sum(sector_exp_vals) != s.num_amides:
                print("ERROR: sector ", s.seq, "does not have right number of amides,", s.num_amides)
                print("Overlap not correctly calculated. Sector exp_vals:", sector_exp_vals)
                print("Method:", method)
                sys.exit(1)

                              
            # cast exp_grid to self.exp_grid for each sector
            long_sector_exp_vals=numpy.zeros(len(self.exp_grid)).astype(int)
            # first and last are identical for each.
            long_sector_exp_vals[0]=sector_exp_vals[0]
            long_sector_exp_vals[-1]=sector_exp_vals[-1]

            #Cast middle values    

            #something is wrong with these indicies...need to fix. long_sector -1 is getting overwritten
            factor = float(len(sector_exp_vals))/(float(len(long_sector_exp_vals))) 

            for i in range(1,len(sector_exp_vals)-1):
                for j in range(1,len(long_sector_exp_vals)-1):
                    if j * factor <= i and (j+1)*factor >= i:
                        ran_int = randint(0,1)
                        long_sector_exp_vals[j+ran_int] += sector_exp_vals[i]
                        break
                    #print(j*factor, i,long_sector_exp_vals[j], sector_exp_vals[i])

            #print(s.seq, sector_exp_vals, method)
            #print(s.seq, long_sector_exp_vals)

            if numpy.sum(long_sector_exp_vals) != s.num_amides:
                print("ERROR: sector ", s.seq, "does not have right number of amides,", s.num_amides)
                print(sector_exp_vals, "not correctly cast to", long_sector_exp_vals)
                print("Method:", method)
                sys.exit(1)

            for n in range(s.start_res, s.end_res+1):
                if self.seq[n]=='P' or self.seq=='p':
                    exp_seq[n]=0
                else:
                    for i in range(len(self.exp_grid)):
                        if long_sector_exp_vals[i] != 0:
                            exp_seq[n] = i+1
                            long_sector_exp_vals[i] = long_sector_exp_vals[i] - 1
                            break
            #print("SSSS", s.seq, sector_exp_vals, exp_seq[s.start_res:s.end_res+1])

        #for f in frags:
            #print(f.seq, f.exp_vals)

        #for n in range(len(self.exp_seq)):
            #print(n,self.state.seq[n], exp_seq[n])
        #print(exp_seq)

        self.exp_seq = exp_seq
        return exp_seq
                    

                #print("Fragment Amides:",s.num_amides,"Overlap:",overlap, "Diff:" , s.num_amides - numpy.sum(overlap), s.seq)

    '''

    def calc_overlap_score(self, seq):
        coverage = self.state.get_coverage(self.frags)
        seclist = self.state.get_sectors().sort
        sorted_sec = 

            if coverage == 0 or self.state.seq[n] == 'P':
                seq[n] = 0
            else
                for f in self.frags:
                    if n > f.start_res + 2 and n < f.end_res:
                    
    '''   
        

    def get_sigma(self, inval, sig):
        return inval*sig+0.05

    def gaussian_model(self,exp,model,sig):
        #print("hello")
        return numpy.exp(-((model-exp)**2)/(2*sig**2))/(sig*numpy.sqrt(2*numpy.pi))

    def lognormal_model(self,exp,model,sig):
        return numpy.exp(-( (numpy.log(model)-numpy.log(exp) )**2)/(2*sig**2))/(sig*numpy.sqrt(numpy.pi)*model)

    def truncated_gaussian_factor(self,exp,a,b,sig):
        return 1/ ( 0.5 * ( scipy.special.erf( (b-exp)/sig * numpy.sqrt(3.1415) ) - scipy.special.erf( (a-exp)/sig * numpy.sqrt(3.1415) ) )  )

    def set_exp_grid(self, grid):
        self.exp_grid=grid

    def calc_exp_grid(self, num_bins, first_time_point=10, last_time_point=3600):

        if num_bins < 2:
            raise Exception("get_exp_bins: num_bins must be greater than 1" )
        # slow_exp_val should satisfy 0.99 = exp(-10**slowest_rate_bin*first_time_point) 
        # fast_exp_val should satisfy 0.01 = exp(-10**fastest_rate_bin*last_time_point) 

        if first_time_point == 0.0:
            first_time_point = 1.0
        slowest_rate_bin=numpy.log10(-numpy.log(0.01)/first_time_point)
        fastest_rate_bin=numpy.log10(-numpy.log(0.99)/last_time_point)

        print("RATE BINS", slowest_rate_bin, fastest_rate_bin, first_time_point, last_time_point)

        return numpy.linspace(fastest_rate_bin, slowest_rate_bin, num_bins)

    def get_sector_averaged_protection_values(self, seq, sectors):
        sector_averaged_protection_values=[0.]*len(seq)
        for s in sectors:
            sum_exp=0
            for n in range(s.start_res, s.end_res+1):
                #print n, seq[n], sum_exp, self.exp_seq[n], self.exp_grid[int(seq[n])-1]
                sum_exp=sum_exp+self.exp_grid[int(seq[n])-1]
            frac_exp=sum_exp/s.num_amides
            for n in range(s.start_res, s.end_res+1):
                if self.seq[n]!='p' and self.seq[n]!='P':
                    sector_averaged_protection_values[n]=frac_exp
        return sector_averaged_protection_values

    def get_indiv_sector_averaged_protection_values(self,seq, sectors):
        sector_averaged_protection_values=[0.]*len(sectors)
        for snum in range(len(sectors)):
            s = sectors[snum]
            sum_exp=0
            for n in range(s.start_res, s.end_res+1):
                #print n, seq[n], sum_exp, self.exp_seq[n], self.exp_grid[int(seq[n])-1]
                sum_exp=sum_exp+self.exp_grid[int(seq[n])-1]
            frac_exp=sum_exp/s.num_amides
            sector_averaged_protection_values[snum]=frac_exp
        return sector_averaged_protection_values

    def get_protection_values(self,seq):
        protection_values=[0.]*len(seq)
        for n in range(len(seq)):
            protection_values[n]=self.exp_grid[self.exp_seq[n]-1]
        return protection_values

    def get_protection_bins(self,seq):
        protection_values=[0.]*len(seq)
        for n in range(len(seq)):
            protection_values[n]=self.exp_seq[n]
        return protection_values

    def get_frag_deut_comparison(self, frags, exp_seq="none"):
        if exp_seq=="none":
            exp_seq=self.exp_seq
        for f in frags:
            for t in f.timepoints:
                model_deut=0
                grid=[]
                #print len(grid), f.start_res, f.end_res
                for r in range(f.start_res+2-self.offset, f.end_res-self.offset):
                    kr=10**self.exp_grid[exp_seq[r]-1]
                    model_deut=model_deut+(1-numpy.exp(-kr*t.time))*self.sat
                    grid.append(exp_seq[r]-1)
                #print(f.seq, t.time, grid, model_deut)
        return model_deut

    def calc_exp_chiscores_from_file(self, frags, infile, outdir):
        self.import_model_deuteration_from_file(frags, infile, outdir)
        for f in frags:
            chi=0
            totd=0
            for t in f.timepoints:
                t.avg_sd_model()
                for d in t.deut:
                    chi=chi+(d-t.model_avg)**2/self.sigma
                    totd=totd+1
                #print f.seq, chi, t.time, t.model_avg, t.deut
            f.chi=chi/totd
            #print f.seq, chi, [t.model_avg for t in f.timepoints], [t.model_avg for t in f.timepoints]

    def calc_exp_chiscores_from_gridvals(self, frags, infile, outdir):
        self.import_model_deuteration_from_gridvals(frags, infile, outdir)
        for f in frags:
            chi=0
            totd=0
            for t in f.timepoints:
                t.avg_sd_model()
                for d in t.deut:
                    chi=chi+(d-t.model_avg)**2/self.sigma
                    totd=totd+1
                #print f.seq, chi, t.time, t.model_avg, t.deut
            f.chi=chi/totd
            print(f.seq, chi)

    def get_model_distributions(self, frags):
        for f in frags:
            for t in f.timepoints:
                print(t.model)

    def get_model_deuteration(self, time, frag, exp_seq):
        model_deut=0
        for r in range(frag.start_res+1, frag.end_res):
            #print f.seq, t.time, r, self.seq[r]
            if self.seq[r] !='P' and self.seq[r] != 'p':
                kr=10**self.exp_grid[int(float(exp_seq[r])-1)]
                model_deut=model_deut+(1-numpy.exp(-kr*time))
        return model_deut/frag.num_observable_amides*100


    def import_model_deuteration_from_file(self, frags, infile, firstline=1, lastline=-1, append=False):
        data=open(infile,"r")
        if append==False:
            self.exp_models=[]
        for line in data:
            #print line
            exp=map(int, line.split(' '))
            self.exp_models.append(exp)
            for f in frags:
                #print f.seq
                for t in f.timepoints:
                    model_deut=0
                    if append==False:
                        t.clear_model_deuteration_values()
                        append=True
                    for r in range(f.start_res+1, f.end_res):
                        #print f.seq, t.time, r, self.seq[r]
                        if self.seq[r] !='P' and self.seq[r] != 'p':
                            kr=10**self.exp_grid[int(float(exp[r])-1)]
                            model_deut=model_deut+(1-numpy.exp(-kr*t.time))
                            #print f.seq, self.seq[r], r, f.start_res, r-f.start_res-1, (1-numpy.exp(-kr*t.time))*self.sat
                    model_deut=model_deut/f.num_observable_amides*100
                    t.add_model_deuteration_value(model_deut)
                    #print(f.seq, t.time, len(t.models), model_deut, t.get_model_avg(), t.get_model_sd(), [r.deut for r in t.replicates])
                    #print f.seq, t.time, t.get_model_avg(), t.get_model_sd(), t.get_avg_sd()

    def import_model_deuteration_from_gridvals(self, frags, models):
        for f in frags:
            timepoints=f.timepoints
            for t in timepoints:
                for m in models:
                    model_deut=0
                    for r in range(f.start_res+self.offset+2, f.end_res+self.offset+1):
                        #print r, f.seq, f.start_res+1, f.end_res, self.seq[r]
                        if self.seq[r] !='P' and self.seq[r] != 'p':
                            #print(r, self.seq[r], m[r], "|", m)
                            kr=10**self.exp_grid[int(float(m[r])-1)]
                            model_deut+=(1-numpy.exp(-kr*t.time))*t.replicates[0].sat
                            #print(">>", r, m[r], int(float(m[r])), self.exp_grid[int(float(m[r])-1)])
                            #print f.seq, self.seq[r], r, f.start_res, r-f.start_res-1, (1-numpy.exp(-kr*t.time))*self.sat
                    #print(model_deut, t.time, t.sigma, f.seq, f.start_res, f.end_res, m[f.start_res+2+self.offset:f.end_res+1+self.offset], m, self.exp_grid)
                    model_deut=model_deut/f.num_observable_amides*100
                    t.add_model_deuteration_value(model_deut)
                #t.avg_sd_model()
                #print t.time, t.model_avg, t.model_sd, t.deut

    def import_models_from_gridvals(self, models, append=False):
        if append=='False':
            self.exp_models=[]
        for i in models:
            self.exp_models.append(i)


    def import_models_from_file(self, modelfile, append='False'):
        data=open(modelfile,"r")
        if append=='False':
            self.exp_models=[]
        for line in data:
            self.exp_models.append(map(int, line.split(' ')))

    def calc_model_scores(self, frags=None, sig=1.0, error_model="gaussian"):
        if self.exp_models==[]:
            print("No models to calculate score")
            return
        if frags is None:
            frags=self.frags
            
        self.model_scores=[]
        for m in self.exp_models:
            self.exp_seq=deepcopy(m)
            score=self.calculate_bayesian_score(frags,sig,error_model)
            self.model_scores.append(score)
        return self.model_scores

    def import_scores(self, scorefile):
        try:
            data=open(scorefile,"r")
            for line in data:
                self.scores.append(line.split().strip()) 
        except:
            for i in scorefile:
                self.scores.append(i)  

    def get_sigma_score(self,sigma0,sig):
        return (2*sigma0/ (numpy.sqrt(numpy.pi) * sig**2 ))*numpy.exp(-sigma0**2 /sig**2)

    def get_model_average(self):
        #print(len(self.exp_models))
        if len(self.exp_models) > 0:
            #print(len(self.exp_models), len(self.exp_models[0]), numpy.array(self.exp_models).shape, self.exp_models[0])
            modelavg=numpy.average(numpy.array(self.exp_models), axis=0)
            return modelavg
        else:
            print("No models to get the average of =(")

    def unimodal_prior(self, omj, sigma0=0.1):
        return (2.0 * omj / ( numpy.sqrt(numpy.pi) * omj**2 )) * numpy.exp(-sigma0**2 / omj**2)
        
    def get_nearest_value(self, value, array):
        '''
        Given a value and 1D numpy array, return the index of the closest value
        '''
        return (numpy.abs(array-value)).argmin(), 

    def calculate_marginal_probabilities(self, error_model="truncated_gaussian"):
        # Calculates marginal probabilities given grids of omega, sigma0, fmod, and fexp

        self.prob_grid=numpy.zeros( (len(self.sigma0_grid), len(self.fmod_grid), len(self.fexp_grid)) )

        for fm in range(len(self.fmod_grid)):
            fmod=self.fmod_grid[fm]

            for f in range(len(self.fexp_grid)):
                fexp=self.fexp_grid[f]

                for s in range(len(self.sigma0_grid)):
                    sigma0=self.sigma0_grid[s]
                    cumul=0

                    for j in range(1,len(self.omega_grid)):
                    # We're going to calculate the likelihoods at each point
                        omj=self.omega_grid[j]
                        omjm1=self.omega_grid[j-1]
                        priorj = self.unimodal_prior(omj, sigma0)
                        priorjm1 = self.unimodal_prior(omjm1, sigma0)
                        dom=omj-omjm1

                        
                        if error_model=="truncated_gaussian":
                            # If we are using the truncated gaussian, calculate those pre-factors
                            # use the truncated gaussian factors at -0.1 and 1.2
                            factor=self.truncated_gaussian_factor(fexp,-10,120,omj)
                            factorm1=self.truncated_gaussian_factor(fexp,-10,120,omjm1)

                        else:
                            factor=1.0
                            factorm1=1.0


                        # Calculate the likelihood at each point
                        pj=self.gaussian_model(fexp,fmod,omj)*factor*priorj
                        pjm1=self.gaussian_model(fexp,fmod,omjm1)*factorm1*priorjm1

                        cumul = cumul+(pj+pjm1)/2.0/dom;
                        #print(pj, pjm1, omj, omjm1, sigma0, "| ", factor, factorm1, "| ", priorj, priorjm1, cumul)
                    #print(fexp, fmod, fexp-fmod, sigma0, "| ", cumul, -1.0 * numpy.log(cumul))
                        #print j,fexp,fmod,omj,(pj+pjm1)/2.0/dom,cumul, cumul2, cumul2/cumul
                    self.prob_grid[s,fm,f]=-1.0 * numpy.log(cumul)
        return self.prob_grid

    def get_marginalized_sigma_probabilities(self,sigma0,error_model="gaussian"):
        '''
        Precomputes marginal probabilities for sigma levels using a cauchy distribution
        @param delfmod - difference between forward model and noise model in delta pct-D
                        units to normalize for fragment length (num amides)
        @param sigma0 - Estimated standard deviation of the instrument in pct-D units
        
        '''
        self.prob_grid=numpy.zeros( (len(self.fmod_grid), len(self.fexp_grid)) )


        for fm in range(len(self.fmod_grid)):
            fmod=self.fmod_grid[fm]
            # Unimodal Distributions for sigma
            #
            for f in range(len(self.fexp_grid)):
                fexp=self.fexp_grid[f]
                if error_model=="truncated_gaussian":
                    # If we are using the truncated gaussian, calculate those pre-factors
                    # use the truncated gaussian factors at -0.1 and 1.2
                    factor=self.truncated_gaussian_factor(fexp,-10,120,omj)
                    factorm1=self.truncated_gaussian_factor(fexp,-10,120,omjm1)
                    o=1
                else:
                    factor=1.0
                    factorm1=1.0

                for j in range(len(self.omega_grid)):
                    # We're going to calculate the likelihoods at each point
                    omj=self.omega_grid[j]
                    omjm1=self.omega_grid[j-1]
                    priorj = self.unimodal_prior(omj, sigma0)
                    priorjm1 = self.unimodal_prior(omjm1, sigma0)
                    dom=omj-omjm1

                    cumul=0

                    # Calculate the likelihood at each point
                    pj=self.gaussian_model(fexp,fmod,omj)*factor*priorj
                    pjm1=self.gaussian_model(fexp,fmod,omjm1)*factorm1*priorjm1

                    cumul = cumul+(pj+pjm1)/2.0/dom;

                    #print(j, fexp, fmod, sigma0, omj,(pj+pjm1)/2.0/dom,cumul)

                self.prob_grid[fm,f]=-1.0 * numpy.log(cumul)
        #print "OMEGA: ",fexp,fmod,omega0,cumul,cumul2,cumul2/cumul
        return self.prob_grid




