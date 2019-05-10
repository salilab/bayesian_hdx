"""
   Sampling functions for HDX models
"""
from __future__ import print_function
import hdx_models
import scoring
import numpy
import numpy.random
import math
import scipy
import system
from random import shuffle
from numpy import linalg
from copy import deepcopy
import os.path
from itertools import combinations_with_replacement
from itertools import permutations

def metropolis_criteria(old_score, new_score, temp):
    # Returns True if move is accepted, False if move is not accepted
    if new_score < old_score:
        return True
    rand_num = numpy.random.rand()
    boltzmann_prob = numpy.exp(-1*(old_score - new_score)/temp) - 1
    #print("MCMC", boltzmann_prob, rand_num, old_score, new_score)
    if boltzmann_prob < rand_num:
        return True
    else:
        return False


class SampledInt(object):
    def __init__(self, allowed_range, random=True, adjacency=3, is_sampled=True):
        self.range = allowed_range
        self.random = random
        self.adjacency = adjacency

    def set_adjacency(self, is_adjacency, adjacency):
        if is_adjacency:
            self.random = False
            self.adjacency = adjacency
        else:
            self.random = True

    def set_range(self, allowed_range):
        self.allowed_range = allowed_range

    def propose_move(self, previous):
        if self.random:
            # Remove the previous value from the list of potential moves
            this_range = [s for s in self.range if s != previous]
            randint = numpy.random.randint(0, len(this_range))
            return this_range[randint]

        else:
            old_index = self.range.index(previous)
            new_index = -1
            while new_index < 0 or new_index >= len(self.range):
                sign = numpy.random.randint(0,2) * 2 - 1
                magnitude = numpy.random.randint(1, self.adjacency+1)
                new_index = int(old_index + magnitude * sign)
            return self.range[new_index]

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
            #print(self.upper_bound, self.lower_bound)
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

class MCSampler(object):
    '''
    Base class that handles importing the system and writing output.

        1) Import system / model / data
        2) Gather all sectors / movers
        3) **define sampling scheme**
        4) Communicate results to hxio.Output
    '''
    def __init__(self, sys, initialize=True, sigma_sample_level=None):
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

        self.sigma_sample_level = sigma_sample_level
        self.sigma_sampler = SampledFloat(0.1, 20, 0.5)
        # Recalculates all sectors and timepoint data.
        if initialize:
            for s in self.states:
                s.initialize(9)

    def run_benchmark(self):
        import time
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

    def run(self, NSTEPS=100, temperature=10, write=False, write_all=False):
        acceptance_total = 1.0
        init_score = 0
        output_files = []
        output = self.states[0].macromolecule.system.get_output()
        for s in self.states:
            init_score += s.calculate_score(s.output_model.model)
            output_files.append(open(output.get_output_file(s), "a"))
            for d in s.data:
                d.collect_times()

        print("Step, score, model_avg_protection_factor, mc_acceptance_ratio")
        for i in range(NSTEPS):
            score, model_avg, acceptance = self.run_one_step(temperature, write_all)
            acceptance_total += acceptance
            print("Step %i, %0.1f | %0.2f, %0.3f" % (i, score, model_avg, acceptance))
            if write:
                for s in range(len(self.states)):
                    st = self.states[s]
                    output.write_model_to_file(output_files[s], st, st.output_model.get_masked_model(st.get_observed_residues()), score, acceptance, sigmas=True)

        acceptance_ratio = acceptance_total/NSTEPS
        print("Average acceptance ratio for this run = ", acceptance_ratio)

        for of in output_files:
            of.close()

    def get_system_model_parameters(self, system):
        # The model parameters are the set of sectors and the system.output_model
        # representation for the exchange factors of these sectors.
        #
        output_model = system.output_model
        system.get_sectors()

    def run_one_step(self, temperature, write=False):
        # For one monte carlo step:
        # Loop over all states in self.states
        total_score = 0
        acceptance_ratio = 0
        for state in self.states:
            # This model is invariant with changes
            init_model = deepcopy(state.output_model.model)
            init_score = state.calculate_score(init_model)
            # Precalculate all the sector scores

            #print(init_score)
            for s in state.sectors:
                #print(s, s.get_residues(), state.observed_residues)
                oldscore = state.calculate_peptides_score(s.get_peptides(), state.output_model.model_protection_factors)
                s.set_score(oldscore)
                #print(s.get_residues(), s.get_score())
                #randint = numpy.random.randint(0,len(s.get_residues()))
                #res = s.get_residues()[randint]
                #state.output_model.change_residue(r, newval)
                # Try changing all residues in sector and recalculating
                '''
                for r in s.get_residues():
                    oldval = state.output_model.get_model_residue(r)
                    newval = self.residue_sampler.propose_move(oldval)
                    # Change the residue in the output_model (changes the PF as well)
                    state.output_model.change_residue(r, newval)
                    # Collect all
                newscore = s.calculate_sector_score(state.output_model.model_protection_factors)
                if oldscore < newscore:
                    #print("Whole Flip Not Accepted", oldscore, newscore)
                    state.output_model.model = init_model
                    s.set_score(oldscore)
                    #print("Whole Flip Not Accepted", oldscore, newscore, state.calculate_score(state.output_model.model), init_score)
                else:

                    #print("Whole Flip Accepted", oldscore, newscore, state.calculate_score(state.output_model.model), init_score)
                    pass
                '''

            resis = state.observed_residues
            shuffle(resis)
            flips = len(resis)
            init_model = deepcopy(state.output_model.model)
            for r in resis:
                # Get the sector that holds this residue
                r_sector = state.residue_sector_dictionary[r]
                oldscore = r_sector.get_score()
                oldval = int(state.output_model.get_model_residue(r))
                newval = self.residue_sampler.propose_move(oldval)

                # First, change the residue in the model
                state.output_model.change_residue(r, newval)
                # Then, change the deuteration state of each timepoint
                state.change_single_residue_incorporation(r, newval)
                # Now, calculate the score again.
                newscore = state.calculate_peptides_score(r_sector.get_peptides(), state.output_model.model_protection_factors)

                accept = metropolis_criteria(oldscore, newscore, temperature)
                #print("SCORE:", oldscore, newscore, accept)
                if not accept:
                    flips -= 1
                    state.output_model.change_residue(r, oldval)
                    state.change_single_residue_incorporation(r, oldval)
                    r_sector.set_score(oldscore)
            flips2 = 0
            for i in range(len(init_model)):
                if init_model[i] != state.output_model.model[i]:
                    flips2 += 1
                cumchange = (state.output_model.model[i] - init_model[i])**2
            #print(math.sqrt(cumchange), init_model)
            #print(math.sqrt(cumchange), state.output_model.model)
            if self.sigma_sample_level is not None:
                for d in state.data:
                    self.sample_sigma(state, temperature)
                    #print(d.get_peptides()[10].get_timepoints()[2].get_sigma(), d.get_peptides()[5].get_timepoints()[2].get_sigma())

            state_score = state.calculate_peptides_score(None, state.output_model.model_protection_factors)

            total_score += state_score

            acceptance_ratio += float(flips)/len(resis)
            #print("XX", oldscore, newscore, state.calculate_peptides_score(r_sector.get_peptides(), state.output_model.model_protection_factors))
            mpf = state.output_model.model_protection_factors
            tot_mpf = 0
            num_mpf = 0
            for pep in state.data[0].get_peptides():
                for i in pep.get_observable_residue_numbers():
                    if not math.isnan(mpf[i-1]) and mpf[i-1] != numpy.inf:
                        tot_mpf += mpf[i-1]
                        num_mpf += 1

                #print(mpf)

        return total_score, state.data[0].get_peptides()[0].get_timepoints()[0].get_model_deuteration(), acceptance_ratio/len(self.states)


    def sample_sigma(self, state, temperature):
        # Sample the sigma values in this dataset.
        # Returns the acceptance boolean
        for dataset in state.data:
            if self.sigma_sample_level == "dataset":
                init_score = dataset.get_score()
                init_sigma = deepcopy(dataset.get_peptides()[0].get_timepoints()[0].get_sigma())
                new_sigma = self.sigma_sampler.propose_move(init_sigma)
                dataset.set_sigma(new_sigma)
                new_score = -1*math.log(dataset.state.scoring_function.experimental_sigma_prior(new_sigma, dataset.sigma_estimate))
                new_score += state.calculate_peptides_score(dataset.get_peptides(), state.output_model.model_protection_factors)

                if not metropolis_criteria(init_score, new_score, temperature):
                    # Reset the sigma back to the original one
                    dataset.set_sigma(init_sigma)
                    # do we need to reset the scores then????


            elif self.sigma_sample_level == "peptide":

                for pep in dataset.get_peptides():
                    init_score = pep.get_score()
                    init_sigma = deepcopy(pep.sigma)
                    new_sigma = self.sigma_sampler.propose_move(init_sigma)
                    dataset.set_sigma(new_sigma)
                    new_score = -1*math.log(dataset.state.scoring_function.experimental_sigma_prior(new_sigma, dataset.sigma_estimate))
                    new_score += dataset.system.calculate_peptides_score([pep], state.output_model.model_protection_factors)
                    if not metropolis_criteria(init_score, new_score, temperature):
                        # Reset the sigma back to the original one
                        pep.set_sigma(init_sigma)
                        # do we need to reset the scores then????

            elif self.sigma_sample_level == "timepoint":
                for pep in dataset.get_peptides():
                    for tp in pep.get_timepoints():

                        init_sigma = deepcopy(tp.get_sigma())
                        tp_model_deut = tp.model_deuteration
                        #init_score = tp.get_score()
                        init_score = -1*math.log(state.scoring_function.experimental_sigma_prior(init_sigma, dataset.sigma_estimate))
                        for rep in tp.get_replicates():
                            init_score += state.scoring_function.replicate_score(tp_model_deut, rep.deut, init_sigma)
                        #print(pep.sequence, tp.time)
                        new_sigma = self.sigma_sampler.propose_move(init_sigma)
                        #if new_sigma is not None:
                        #print(new_sigma, tp.get_sigma())
                        new_score = -1*math.log(state.scoring_function.experimental_sigma_prior(new_sigma, dataset.sigma_estimate))

                        for rep in tp.get_replicates():
                            new_score += state.scoring_function.replicate_score(tp_model_deut, rep.deut, new_sigma)
                        #print(init_score, new_score)

                        if metropolis_criteria(init_score, new_score, temperature):
                            tp.set_sigma(new_sigma)
                            tp.set_score(new_score)
                            #print("NEW SIGMA", pep.sequence, tp.time, "||", new_sigma, init_sigma, new_score, init_score)
                            #print("NO CHANGE IN SIGMA")

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




def do_mc_sampling(model, temp=1, sigma=5.0, error_model="gaussian", NSTEPS=100, save_results=False,
                    outdir="./", outfile_prefix="", sample_sigma=False, noclobber=True, print_t=10):


    # First create output directory
    if save_results:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        elif noclobber:
            raise Exception("Output directory ", outdir, " already exists")
        else:
            print("WARNING: Running this sampling protocol will potentially overwrite results in directory ", outdir)
        for s in model.states:
            s.sigmafile=outdir + outfile_prefix + "_" + s.state_name + "_s" +str(sigma)+ "_T" +str(temp) + "_sigmas.dat"
            outfile3= open(s.sigmafile, "w")

            string="#"
            for f in s.exp_model.frags:
                for tp in f.timepoints:
                    string+=f.seq+"_"+str(tp.time)+", "

            string+="\n"
            outfile3.write(string)
            outfile3.close()


    model_exp_sequences=numpy.zeros((NSTEPS, len(model.seq), len(model.states)))
    model_scores=numpy.zeros((NSTEPS, len(model.states)))
    for t in range(NSTEPS):
        timepoint_seq=[]
        ns = 0
        step_scores=[]
        for s in model.states:

            exp_model=s.exp_model
            score=exp_model_sector_MC_step(s, exp_model, exp_model.frags, temp, sigma, error_model, sample_sigma)
            model_exp_sequences[t,:,ns]=exp_model.exp_seq
            model_scores[t,ns]=score
            step_scores.append(score)
            ns+=1
            sigs=[]
            for f in exp_model.frags:
                for tp in f.timepoints:
                    sigs.append(tp.sigma)
                    if save_results:
                        deut=exp_model.get_model_deuteration(tp.time,f,exp_model.exp_seq)
                        tp.add_model_deuteration_value(deut)
                #print(f.seq, f.start_res, f.end_res, score, sum(exp_model.exp_seq))

            #print(['%0.2f' % val for val in exp_model.get_indiv_sector_averaged_protection_values(exp_model.exp_seq,s.sectors)])
            #print("sig", numpy.average(sigs), numpy.std(sigs), numpy.max(sigs), numpy.min(sigs))
            if save_results:
                outfile3 = open(s.sigmafile, "ab")
                numpy.savetxt( outfile3 , numpy.atleast_2d(sigs), "%f "*len(sigs))
                outfile3.close()
        if t%print_t == 0:
            string=""
            for i in range(len(step_scores)):
                #print(i, model.states[i].state_name)
                string+=model.states[i].state_name + "_score= " + str(step_scores[i]) + " "
            print(t, string)

    if save_results:
        state_num = 0
        ns = 0
        for s in model.states:
            s.modelfile=outdir + outfile_prefix + "_" + s.state_name + "_s" +str(sigma)+ "_T" +str(temp) + "_exp-seq.dat"
            s.scorefile=outdir + outfile_prefix + "_" + s.state_name + "_s" +str(sigma)+ "_T" +str(temp) + "_score.dat"
            outfile1= open(s.modelfile, "w")
            outfile2= open(s.scorefile, "w")

            numpy.savetxt( outfile1 , model_exp_sequences[:,:,ns], "%i")
            numpy.savetxt( outfile2 , model_scores[:,ns], "%f")

            outfile1.close()
            outfile2.close()

            ns += 1



def exp_model_sector_MC_step(state,exp_model,frags,temp,sigma,error_model="truncated_gaussian", sample_sigma=False):

    # randomly select a position in the chain.  If it is not a proline or have zero coverage, flip

    if sample_sigma:
        sigma=exp_model.sigma
    state_score=0
    flips=0
    nonflips=0
    seq0=exp_model.exp_seq
    #score0=exp_model.calculate_bayesian_score(frags,1.0,sigma,error_model="truncated_gaussian")
    state.get_all_sectors()
    i = 0
    flipped = 0
    for s in state.sectors:
            #print s.seq, [f.seq for f in s.fragments]
        for n in range(len(s.seq)):
            i = i+1
            sector_score=exp_model.calculate_bayesian_score(s.fragments, sigma, error_model)
            ipick = numpy.random.random_integers(0,len(s.seq)-1)
            exp_model.exp_seq=deepcopy(seq0) #sequence from last round
            resnum = s.start_res+ipick
            #ensure it is not proline...just skip it for now
            if s.seq[ipick]=='P' or s.seq[ipick]=='p':
                delta=0
                #print ipick, "is a proline"
                continue
            else:
                delta=numpy.random.random_integers(0,len(exp_model.exp_grid)-1)

                exp_model.exp_seq[resnum] = delta
                new_sector_score = exp_model.calculate_bayesian_score(s.fragments,sigma,error_model)
                randn = numpy.random.random_sample()
                delscore = new_sector_score - sector_score
                if delscore < 0.0:
                    # If change is negative, accept the move
                    sector_score = new_sector_score
                    #exp_model.exp_seq=deepcopy(seq0)
                    seq0 = deepcopy(exp_model.exp_seq)
                    flipped = flipped + 1

                elif randn < numpy.exp(-(delscore)/1):#temp):
                    # Negative values of delscore give energies > 1; so this will always be true
                    # Some positive values of delscore will return a value greater than than randn
                    sector_score = new_sector_score
                    flipped = flipped + 1
                    # Change the seq0 to the modified rate
                    #exp_model.exp_seq=deepcopy(seq0)
                    seq0 = deepcopy(exp_model.exp_seq)

                    #print(s.seq, resnum, "||", old, delta, "|FLIP|", delscore, flipped, "|", randn, numpy.exp(-(delscore)/0.01))
                else:
                    # If randn is higher than

                    exp_model.exp_seq=deepcopy(seq0)

    if sample_sigma==True:

        sample_individual_sigmas(exp_model, frags, temp)
        '''
        oldsig=exp_model.sigma
        oldscore=exp_model.calculate_bayesian_score(frags, oldsig)
        factor=numpy.random.random_sample()*2-1
        newsig = oldsig+factor*0.5  #continuous sampling
        if newsig <= 0:
            newsig=0.1
        #t.sigma=newsig
        newscore=exp_model.calculate_bayesian_score(frags, newsig)
        #print("XXX", f.seq, t.time, t.sigma, factor, oldscore, newscore)
        #Just minimization for now
        if oldscore >= newscore:
            exp_model.set_sigma(newsig)
            state_score=newscore
        else:
            state_score=oldscore
        print("SIG::", oldsig, newsig, oldscore, newscore, exp_model.state.state_name)
        '''
    state_score = exp_model.calculate_bayesian_score(frags,sigma,error_model)

    return state_score

def sample_individual_sigmas(exp_model, frags, temp):
    # For each timepoint for each fragment, propose a step to the sigma and recalculate
    # the timepoint score.

    # Sigma must be greater than 0.1

    for f in frags:
        noa = f.get_num_observable_amides()
        exp_seq = exp_model.exp_seq[f.start_res-1+2:f.end_res]

        for t in f.timepoints:
            sigma=t.sigma
            factor=numpy.random.random_sample()*2-1
            newsig = sigma+factor*0.5  #continuous sampling
            if newsig <= 0.1:
                newsig=0.1

            randn=numpy.random.random_sample()*2

            # Sum up the integers in exp_model.exp_grid

            freq_grid = numpy.zeros(len(exp_model.exp_grid))

            for i in range(len(exp_model.exp_grid)):
                freq_grid[i]=exp_seq.count(i+1)

            #print(f.seq, f.start_res, exp_seq, exp_model.exp_seq, freq_grid, len(exp_model.exp_grid))

            oldscore=t.calculate_tp_score(freq_grid, exp_model.exp_grid, sigma, noa)-1.0*numpy.log(exp_model.unimodal_prior(sigma))
            newscore=t.calculate_tp_score(freq_grid, exp_model.exp_grid, newsig, noa)-1.0*numpy.log(exp_model.unimodal_prior(newsig))

            if randn < numpy.exp(-(newscore-oldscore)/temp):
                #print("XX", f.seq, t.time, sigma, factor*0.5, " || ", oldscore, oldscore-newscore, numpy.exp(-(newscore-oldscore)/temp), " | ", randn, " | ", t.model_deut, t.get_avg_sd(), [(r.deut, -1*numpy.log(r.calculate_replicate_score(t.model_deut, sigma))) for r in t.replicates])
                t.sigma=newsig
            #else:
            #    print("--", f.seq, t.time, sigma, factor*0.5, " || ", oldscore, oldscore-newscore, numpy.exp(-(newscore-oldscore)/temp), " | ", randn, " | ", t.model_deut, t.get_avg_sd(), [(r.deut, -1*numpy.log(r.calculate_replicate_score(t.model_deut, sigma))) for r in t.replicates])


def enumerate_fragment(frag, exp_grid, sig, num_models = 1):
    """Enumerates and scores all possible models for the given an exp_grid
    returns the top num_models scoring exp grids"""
    n = frag.num_observable_amides
    nbin = len(exp_grid)
    num = n
    possible_number_combinations = list(combinations_with_replacement(range(n), nbin))
    #all_possible_combinations = list(product(range(n), repeat=nbin))
    score = 0
    minscore = 1.0*10^34
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
