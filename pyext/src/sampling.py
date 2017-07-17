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
                if new_sector_score-sector_score < 0.0:
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
    minscore = 10000
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

