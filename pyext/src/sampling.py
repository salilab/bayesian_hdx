"""
   Sampling functions for HDX models
"""
from __future__ import print_function
import scoring
import numpy
import numpy.random
import math
import scipy
import system
from random import shuffle
from numpy import linalg
from copy import deepcopy
import os
from itertools import combinations_with_replacement
from itertools import permutations

# if exp(-diff/temp) > rand: accept

def metropolis_criteria(old_score, new_score, temp, proposal_ratio=0.4):
    # Returns True if move is accepted, False if move is not accepted
    #print(old_score, new_score)
    if new_score <= old_score:
        return True
    rand_num = numpy.random.rand()
    boltzmann_prob = numpy.exp(-1*(new_score - old_score)/temp)
    #print("  MCMC", boltzmann_prob * proposal_ratio > rand_num, boltzmann_prob, rand_num, old_score-new_score, temp)
    if boltzmann_prob * proposal_ratio > rand_num:
        # Upward step taken
        return True
    else:
        return False

# Connect this object to a Residue
class SampledInt(object):
    def __init__(self, allowed_range, random=True, adjacency=3, is_sampled=True, uniform=True):
        self.range=allowed_range # range of values this integer can be
        self.random=random # Are we doing random assignment?
        self.adjacency=adjacency # What is the range of values that delta can be?
        self.moved=False
        self.uniform=uniform # Choose step magnitude from uniform distribution or Normal distribution (False)

    def initialize():
        self.propose_move()

    def set_adjacency(self, is_adjacency, adjacency, uniform=True):
        if is_adjacency:
            self.random = False
            self.adjacency = adjacency
        else:
            self.random = True
        self.uniform = uniform

    def set_range(self, allowed_range):
        self.range=allowed_range

    def set_moved(self):
        self.moved=True

    def propose_move(self, previous):
        if self.moved:
            print("This mover has already been moved")

        if self.random:
            # Remove the previous value from the list of potential moves
            this_range = [s for s in self.range if s != previous]
            new_index = numpy.random.randint(0, len(this_range))

            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # Instead of this, change the value in the residue object
            # self.object.set_value(new_index)
            return this_range[new_index]

        else:
            self.old_index = self.range.index(previous)
            new_index = -1
            while new_index < 0 or new_index >= len(self.range):
                sign = numpy.random.randint(0,2) * 2 - 1
                magnitude = self.get_magnitude()
                new_index = int(self.old_index + magnitude * sign)

            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # Instead of this, change the value in the residue object
            # self.object.set_value(new_index)
            return self.range[new_index]

        self.moved=True

    def get_magnitude(self):
        if self.uniform:
            return numpy.random.randint(1, self.adjacency+1)
        else:
            return int(numpy.ceil(abs(numpy.random.normal(0,self.adjacency))))

    def accept(self):
        # Keep the residue object the same
        self.moved=True

    def reject(self):
        # Change the residue object back to its original value
        # NOT IMPLEMENTED!!!
        # self.object.set_value(self.range[self.old_index])
        self.moved=False


class SampledFloat(object):
    def __init__(self, lower_bound, upper_bound, maxdel, distribution="uniform", random=False, is_sampled=True):
        self.upper_bound = upper_bound
        self.lower_bound = lower_bound
        self.maxdel = maxdel
        self.random = random

    def propose_move(self, previous):
        if self.random:
            new_value = numpy.random.rand() * (self.upper_bound - self.lower_bound) + self.lower_bound
            return new_value

        else:
            new_value = self.lower_bound - 1
            while new_value >= self.upper_bound or new_value <= self.lower_bound:
                sign = numpy.random.randint(0,2) * 2 - 1
                magnitude = numpy.random.rand() * self.maxdel 
                new_value = previous + sign * magnitude
            #print("PROPOSE_MOVE:", new_value, previous, sign, magnitude, self.upper_bound, self.lower_bound)
            return new_value
            '''
            if new_value <= self.upper_bound and new_value >= self.lower_bound:
                print("PROPOSE_MOVE:", new_value, previous, sign, magnitude, self.upper_bound, self.lower_bound)
                return new_value
            else:
                print("NO_MOVE")
                self.propose_move(previous)
            '''

class EnumerationSampler(object):
    def __init__(self, sys, sigma_sample_level=None):
        if type(sys) is system.State:
            states = [sys]
        elif type(sys) is system.Macromolecule:
            states = sys.get_states()
        elif type(sys) is system.System:
            states = sum([mol.get_states() for mol in sys.get_macromolecules()], [])
        else:
            raise Exception("You must pass a System, State or Macromolecule object")

        self.states = states 
        
        m = self.states[0].output_model 

    def uun(self, write=False):
        """Enumerates and scores all possible models for the given an exp_grid
        returns the top num_models scoring exp grids"""
        output = self.states[0].macromolecule.system.get_output()
        acceptance = 0.0
        # Get the residue types
        for s in self.states:
            s.initialize(1)

            if write:
                output_file = open(output.get_output_file(s), "a")

            resis = s.observed_residues

            n = len(resis) 
            exp_grid = s.output_model.get_sampling_grid()
            nbin = len(exp_grid)
            num = n
            possible_number_combinations = list(combinations_with_replacement(range(1,nbin+1), n))
            #all_possible_combinations = list(product(range(n), repeat=nbin))
            print("Assessing", possible_number_combinations, "combinations")
            for model in possible_number_combinations:

                for n in range(len(resis)):
                    resid = resis[n]
                    mod_val = model[n]
                    s.output_model.change_residue(resid, mod_val)
                    s.change_single_residue_incorporation(resid, int(mod_val))

                score = s.calculate_score(s.output_model.model)

                if write:
                    output.write_model_to_file(output_file, s, s.output_model.get_masked_model(s.get_observed_residues()), score, acceptance, sigmas=True)
            if write:
                output_file.close()

class MCSampler(object):
    '''
    # 
    
    '''
    def __init__(self, sys, initialize=True, 
                sigma_sample_level=None, 
                pct_moves=25, 
                accept_range=(0.3, 0.8),
                num_independent_sims=1):
        # Ensure that all states in system has a dataset and a model and a scoring function
        '''
        @param sigma_sample_level - None: Don't sample sigmas. "dataset": sample one sigma per dataset. 
                            "peptide": 1 sigma per peptide, "timepoint": 1 sigma per timepoint
                            "replicate": 1 sigma per replicate_id
        '''
        if type(sys) is system.State:
            states = [sys]
        elif type(sys) is system.Macromolecule:
            states = sys.get_states()
        elif type(sys) is system.System:
            states = sum([mol.get_states() for mol in sys.get_macromolecules()], [])
        else:
            raise Exception("You must pass a System, State or Macromolecule object")

        self.states = states

        m = self.states[0].output_model

        if m.sampler_type == "int":
            self.residue_sampler = SampledInt(m.sampler_range())
        elif m.sampler_type == "float":
            raise Exception("Sampler.__init__: Floating point representation for residues is not implemented yet")

        # Should put Sigma in the scoring function
        self.sigma_sample_level = sigma_sample_level
        self.sigma_sampler = SampledFloat(0.1, 20, 0.05)

        self.pct_moves = pct_moves
        
        self.acceptance_range = accept_range
        
        self.num_sims = num_independent_sims
        # Recalculates all sectors and timepoint data.
        if initialize:
            for s in self.states:
                s.initialize()


    def run_exponential_temperature_decay(self, tmax=100, tmin=2.0, 
                                annealing_steps=200, 
                                steps_per_anneal=1, 
                                write=False, 
                                adjacency=10,
                                uniform_adjacency=True,
                                quiet=False):
        '''
        A simulated annealing that starts from high temperature and gradually cools
        to a low temperature
        '''
        
        # Set how far each residue can move.
        self.residue_sampler.set_adjacency(True, adjacency, uniform=uniform_adjacency)

        print("********")
        print("Starting Exponential Temperature Decay from")
        print("T =", tmax, "to T =", tmin, "over", annealing_steps, "steps.")
        print("********")
        deltaT = math.exp(math.log(tmin/tmax)/annealing_steps)
        for s in range(annealing_steps):

            Tm = tmax * (deltaT ** s)
            for i in range(steps_per_anneal):
                score, model_avg, accept = self.run_one_step(Tm, write)

            if not quiet:
                print('Temp: %2.2f | Score: %2.1f' % (Tm,score))


    def run_benchmark(self):
        import time
        print("Starting Benchmark")
        '''Run 100 MC steps and report the time to run 1000 steps
        '''       
        times = []
        for i in range(10):
            start = time.time()
            self.run(NSTEPS=10)
            end = time.time()
            times.append((end-start)*100)
        time = numpy.average(times)
        sd = numpy.std(times)
        print("This system will take about ", int(time), "+/-", int(sd*2) , " seconds per 1000 steps")

    def get_acceptable_temperature(self, init_temp, acceptance_range, nsteps=100):
        accept = False
        low = False
        high = False

        orig_pct_moves = self.pct_moves
        self.pct_moves = 100

        scale = 1.0

        temp = init_temp

        while not accept:
            acceptance_total=1.0
            init_temp = temp
            for i in range(nsteps):
                #print(i,)
                score, model_avg, acceptance = self.run_one_step(temp)
                #print(score, model_avg, acceptance)
                acceptance_total += acceptance

            acceptance_ratio = acceptance_total/nsteps

            if acceptance_ratio < acceptance_range[0]:
                # We have an acceptance ratio that is too low. Raise the temperature
                low = True
                if high:
                    # If we were previously too high, raise the divisor by a factor of two    
                    high = False
                    scale = scale * 2.0
                    
                temp = temp + temp / scale

            elif acceptance_ratio > acceptance_range[1]:
                # We have an acceptance ratio that is too high.  Lower the temperature.
                high = True
                if low:
                    # If we were previously too low, raise the divisor by a factor of two                
                    low = False
                    scale = scale * 2.0

                temp = temp - temp / scale / 2.0

            else:
                accept = True

            print(high, low, accept, scale, init_temp, acceptance_ratio)
        self.pct_moves = orig_pct_moves
        return temp

    def run(self, NSTEPS=100, init_temp=10, write=False, write_all=False, acceptance_range=None, find_temperature=False):
        print("Starting run")
        if acceptance_range is None:
            acceptance_range = self.acceptance_range
        acceptance_total = 1.0       
        init_score = 0
        output_files = []
        output = self.states[0].macromolecule.system.get_output()
        for s in self.states:
            state_score = s.calculate_score(s.output_model.model)
            s.set_score(state_score)
            init_score += state_score
            if write:
                output_files.append(open(output.get_output_file(s), "a"))
            for d in s.data:
                d.collect_times()
        if find_temperature:
            print("Finding initial running temperature")
            temperature = self.get_acceptable_temperature(init_temp=init_temp, acceptance_range=acceptance_range)
        else:
            temperature = init_temp

        print("Simulation temperature =", temperature)
        print("Step score | states_avg_protection_factor | mc_acceptance_ratio")
        for i in range(NSTEPS):
            #print("Step:", i)
            score, model_avg_str, acceptance = self.run_one_step(temperature, write_all)
            acceptance_total += acceptance
            if i%100 == 0:
                print("Step %i, %0.1f | %s, %0.2f" % (i, score, model_avg_str, acceptance))

            if write:
                for s in range(len(self.states)):
                    st = self.states[s]
                    if st.output_model.sample_only_observed_residues == True:
                        model = st.output_model.get_masked_model(st.get_observed_residues())
                    else:
                        model = st.output_model.get_model()

                    output.write_model_to_file(output_files[s], st, model, st.score, acceptance, sigmas=True)

        acceptance_ratio = acceptance_total/NSTEPS
        print("Average acceptance ratio for this run =", round(acceptance_ratio,3), " |  Temp = ", temperature)

        for of in output_files:
            of.close()  

    def run_one_step(self, temperature, write=False):
        # Running one MC step over all model states and sigma parameters
        # Loop over all states in self.states
        total_score = 0
        acceptance_ratio = 0

        for state in self.states:

            init_model = deepcopy(state.output_model.model)
            init_score = state.get_score()

            #(state.name, state.output_model.model[5:15])
            ###########################
            # This should be movable particles
            ###########################

            if state.output_model.sample_only_observed_residues==True:
                resis = state.observed_residues
            else:
                resis = state.get_exchanging_residues()

            shuffle(resis)
            flips = int(max(math.ceil((self.pct_moves * len(resis))/100.), 1))
            # Flip a number of residues


            for r in resis[:flips]:
                # Get the sector that holds this residue           
                #r_sector = state.residue_sector_dictionary[r]
                # Get the current value for this residue
                oldval = int(state.output_model.get_model_residue(r))
                # Propose a new value given the current state
                oldscore = state.get_score()
                #print(r, oldval, oldscore, state.output_model.get_model())
                newval = self.residue_sampler.propose_move(oldval) 

                # Change the residue incorporation values in each sector and calculate the new score:
                state.change_single_residue_incorporation(r, newval)
                newscore = state.calculate_peptides_score(state.get_all_peptides(), state.output_model.get_current_model())
                state.set_score(newscore)

                accept = metropolis_criteria(oldscore, newscore, temperature)
                #print(" ", oldscore, newscore, r, oldval, newval, accept, flips)
                if not accept:
                    flips -= 1
                    state.output_model.change_residue(r, oldval)
                    state.change_single_residue_incorporation(r, oldval)
                    state.set_score(oldscore)

            # Determine the total number of moves performed
            flips2=0
            for i in range(len(init_model)):
                if init_model[i] != state.output_model.model[i]:
                    flips2 += 1

                # Measures the total change in Pf value over all residues. Might be useful.
                #m_sq_change = sum((state.output_model.model - init_model)**2)


            # Now sample the sigma values if we are doing that
            if self.sigma_sample_level is not None:
                for d in state.data:
                    self.sample_sigma(state, temperature)

            # Recalculate the state score and add it to the total score
            #print(state.get_score())
            #final_state_score = state.calculate_peptides_score(state.get_all_peptides(), state.output_model.get_current_model())
            #print("-",final_state_score)
            #state.set_score(final_state_score)
            total_score += state.get_score() #final_state_score

            '''
            mpf = state.output_model.get_current_model()
            tot_mpf = 0
            num_mpf = 0
            for pep in state.data[0].get_peptides():
                for i in pep.get_observable_residue_numbers():
                    if not math.isnan(mpf[i-1]) and mpf[i-1] != numpy.inf:
                        tot_mpf += mpf[i-1]
                        num_mpf += 1 
            '''

            acceptance_ratio += float(flips2)/int(max(math.ceil((self.pct_moves * len(resis))/100.), 1))
            #print(state.name, state.output_model.get_model()[5:15])
            #print("NUMFLIPS", state.name, flips, flips2, float(flips2)/int(max(math.ceil((self.pct_moves * len(resis))/100.), 1)))
        model_avg = [numpy.average(s.output_model.get_current_model()) for s in self.states]
        model_avg_str = ""
        for m in model_avg:
            model_avg_str+=str(round(m,3))+" "

        return total_score, model_avg_str[0:-1], acceptance_ratio/len(self.states)#, m_sq_change #acceptance_ratio/len(self.states)


    def sample_sigma(self, state, temperature):
        # Sample the sigma values in this dataset.
        # Returns the acceptance boolean
        for dataset in state.data:

            init_score = state.get_score()

            if self.sigma_sample_level == "dataset":

                init_score = state.calculate_peptides_score(dataset.get_peptides(), state.output_model.get_current_model())
                
                init_sigma = deepcopy(dataset.get_peptides()[0].get_timepoints()[0].get_sigma())
                new_sigma = self.sigma_sampler.propose_move(init_sigma)

                dataset.set_sigma(new_sigma)
                new_score = state.calculate_peptides_score(dataset.get_peptides(), state.output_model.get_current_model())

                if not metropolis_criteria(init_score, new_score, temperature):
                    # Reset the sigma back to the original one
                    dataset.set_sigma(init_sigma)
                    # do we need to reset the scores then????
                    state.set_score(init_score)

            elif self.sigma_sample_level == "peptide":

                # For each peptide, only need to calculate the score for each peptide
                for pep in dataset.get_peptides():
        
                    init_score = state.calculate_peptides_score([pep], state.output_model.get_current_model())
                    init_sigma = deepcopy(pep.sigma)
                    new_sigma = self.sigma_sampler.propose_move(init_sigma)
                    pep.set_sigma(new_sigma)
                    new_score = state.calculate_peptides_score([pep], state.output_model.get_current_model())
                    if not metropolis_criteria(init_score, new_score, temperature):
                        # Reset the sigma back to the original one
                        pep.set_sigma(init_sigma)
                        # do we need to reset the scores then????

            elif self.sigma_sample_level == "timepoint":
                for pep in dataset.get_peptides():
                    for tp in pep.get_timepoints():

                        init_sigma = deepcopy(tp.get_sigma())
                        tp_model_deut = tp.model_deuteration * 100

                        init_score = 0
                        for rep in tp.get_replicates():
                            init_score += state.scoring_function.forward_model.replicate_score(tp_model_deut, rep.deut, init_sigma)

                        new_sigma = self.sigma_sampler.propose_move(init_sigma)
                        new_score = 0
                        for rep in tp.get_replicates():
                            new_score += state.scoring_function.forward_model.replicate_score(tp_model_deut, rep.deut, new_sigma)

                        if metropolis_criteria(init_score, new_score, temperature):
                            tp.set_sigma(new_sigma)
                            tp.set_score(new_score)

    def perturb_specific_residues(self, state, residues, adjacency="random", metropolis=False):
        # not implemented yet

        # Certain residues get stuck.  We can see this though convergence PSRF analysis.
        # To help this, we should implement a method for perturbing these specific residues.

        # Each state is different.  

        # Perhaps include the metropolis criteria
        for r in residues:
            oldval = int(state.output_model.model[r])
            if adjacency == "random":
                newval = self.residue_sampler.propose_move(oldval, random=True)
            else:
                newval = self.residue_sampler.propose_move(oldval, random=False, adjacency=adjacency)
            state.output_model.model[r] = newval


def benchmark(model, sample_sigma):
    import time
    '''Run 100 MC steps and report the time to run 1000 steps
    '''
    times=[]
    for i in range(10):
        start=time.time()
        do_mc_sampling(model, NSTEPS=10, print_t=1000)
        end=time.time()
        times.append((end-start)*100)
    time=numpy.average(times)
    sd=numpy.std(times)
    print("This system will take about ", int(time), "+/-", int(sd*2) , " seconds per 1000 steps")


def simulated_annealing(model, sigma, sample_sig=True, equil_steps=10000, annealing_steps=100, save_all=False, 
                        outdir="./sampling_output/", outfile_prefix="", print_t=10, sample_sigma=False, noclobber=False):

    for temp in [20, 10, 5, 1]:
        do_mc_sampling(model, temp, sigma, NSTEPS=annealing_steps, sample_sigma=sample_sig, print_t=print_t,
                        save_results=save_all, outdir=outdir, outfile_prefix=outfile_prefix, noclobber=noclobber)

    # Equilibrium run  
    do_mc_sampling(model, 1.0, sigma, NSTEPS=equil_steps, sample_sigma=sample_sigma, print_t=print_t,
                        save_results=True, outdir=outdir, outfile_prefix=outfile_prefix, noclobber=noclobber)


def enumerate_fragment(frag, exp_grid, sig, num_models = 1):
    """Enumerates and scores all possible models for the given an exp_grid
    returns the top num_models scoring exp grids"""
    n = frag.num_observable_amides  
    nbin = len(exp_grid)
    num = n
    possible_number_combinations = list(combinations_with_replacement(range(n), nbin))
    #all_possible_combinations = list(product(range(n), repeat=nbin))
    score = 0
    minscore = 1.0*pow(10,34)

    for i in possible_number_combinations:
        if sum(i) == n:
            numbers = list(i)
            all_possible_permutations = permutations(numbers)
            seen = set()
            for grid in all_possible_permutations:
                if grid in seen:
                    continue
                seen.add(grid)
                score = frag.calculate_frag_score(grid, exp_grid, sig, force=True)
                #print(frag.seq, score, grid)
                sum_exp = 0
                for x in range(len(exp_grid)):
                    sum_exp += exp_grid[x]*grid[x]
                if score < minscore:
                    minmodel = grid
                    minscore = score

    #print(numpy.array(minmodel), minscore)
    return numpy.array(minmodel)

