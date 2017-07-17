"""@namespace IMP.hdx.analysis
   Analysis functions for HDX simulations
"""
from __future__ import print_function
import hdx_models
import input_data
import plots
#from scipy.stats import cumfreq
#from scipy.stats import chi2_contingency
import numpy
import math
import time
import os.path
from pylab import *
from matplotlib import *


class ParseOutputFile(object):
    def __init__(self, output_file, state):
        self.output_file = output_file
        self.state = state
        self.datafiles = []
        self.sectors = []
        self.pf_grids = {}
        self.observed_residues = []
        self.parse_header()

    def parse_header(self):
        '''
        Function that parses the header of output files.
        Stores the logk grid, along with other experimental information
        '''
        f = open(self.output_file)
        for line in f.readlines():
            
            # > means model data (so header is over.)
            if line[0]==">":
                break
            
            # #-symbol meas datasets
            elif line[0:2]=="# ":
                self.datafiles.append( (line[2:].split("|")[0].strip(), float(line[2:].split("|")[1].strip())) )
            
            # @-symbol means sectors
            elif line[0:2]=="@ ":
                for s_string in line[2:].strip().split("|"):
                    sector = []
                    for r in s_string.strip().split(" "):
                        if r != "":
                            sector.append(int(r))
                            self.observed_residues.append(int(r))
                    self.sectors.append(sector)

            elif line[0:2]=="$ ":
                if line[2:].split("|")[0].strip() != "Residue_number":
                    #print(line[2:].split("|")[0].strip())
                    res = int(line[2:].split("|")[0].strip())
                    grid = []
                    for pf in line[2:].split("|")[1].strip().split(" "):
                        #print(line[2:].split("|")[1].strip().split(" "))
                        grid.append(float(pf))
                    self.pf_grids[res] = grid

            elif line.split(":")[0].strip() == "grid_size":
                self.grid_size = int(line.split(":")[1].strip())

            elif line.split(":")[0].strip() == "State":
                self.state_name = line.split(":")[1].strip()

            elif line.split(":")[0].strip() == "Molecule_Name":
                self.molecule_name = line.split(":")[1].strip()

    def get_best_scoring_models(self, N, sigmas=False, return_pf=False):
        ''' Get the N best scoring models from the output file
        Returns a list of tuples of best_scoring_models 
            [(score, [model])]
        and (if sigmas=True)
        a grid of the timepoint sigma values.
        '''
        # Model entries are marked with a > as the first character

        f = open(self.output_file, "r")
        best_scoring_models = []
        # Cycle over all lines
        for line in f.readlines():
            if line[0]==">":
                score = float(line.split("|")[1].strip())

                # if the score is better than the last best score
                if len(best_scoring_models) < N or score < best_scoring_models[-1][0]:
                    if len(best_scoring_models) >= N:
                        del best_scoring_models[-1]
                    model_string = line[1:].split("|")[0].strip()
                    model_list = []
                    for m in model_string.split(" "):
                        model_list.append(int(m))
                    if return_pf:
                        model_list = self.models_to_protection_factors(model_list)
                    best_scoring_models.append((score, model_list))
                    best_scoring_models = sorted(best_scoring_models, key=lambda x: x[0])

        self.best_scoring_models = best_scoring_models
        return best_scoring_models


    def models_to_protection_factors(self, models):
        # Input a list of list of integers.  
        # CHECK THAT THE MODEL SIZE IS CORRECT!
        if type(models[0]) != list:
            models = [models]

        output = []
        for m in models:
            pf_model = []
            for res in range(len(m)):
                if res + 1 in self.observed_residues:
                    #print(res+1, m[res-1], self.pf_grids[res+1])
                    pf_model.append(float(self.pf_grids[res+1][m[res-1]-1]))

                else:
                    pf_model.append(numpy.nan)
            output.append(pf_model)

        return output





def get_best_scoring_models(modelfile, scorefile, num_best_models=100, prefix=None, write_file=True):
    #takes a model and score file and writes a new model file with the best X scoring models.
    #This new file can then be imported into an HDXModel class for analysis
    i=0
    if prefix is None:
        outfile="./best_models.dat"
    else:
        outfile="./" + prefix + "_best_models.dat"

    # You have one chance to not overwrite your best models file
    if os.path.isfile(outfile):
        print("WARNING: ", outfile, " exists, renamed to ", outfile, ".old")
        os.rename(outfile, outfile +".old")

    scores=[]
    models=[]
    top_models=[]
    top_score_indices=[]
    top_scores=[]
    infile=open(scorefile, "r")
    for line in infile:
        scores.append(float(line.split()[0].strip()))
    infile=open(modelfile, "r")

    for line in infile:
        models.append(line)

    for i in range(num_best_models):
        top_score_indices.append(scores.index(min(scores)))
        top_scores.append(min(scores))
        #print(scores, min(scores), scores.index(min(scores)), top_score_indices)
        scores[scores.index(min(scores))]=max(scores)+1
    if write_file:
        output_file=open(outfile, "w")
        return top_models, top_scores
    else:
        for i in top_score_indices:
            top_models.append(map(int, models[int(i)].split()) )
        return top_models, top_scores

def sector_sort(sectors):
    # Given a list of sectors, sort them by residue number
    last_first_res=0
    sorted_sectors=[]

    # get last residue of a sector start
    for s in sectors:
        if s.start_res > last_first_res:
            last_first_res = s.start_res

    for n in range(last_first_res+1):
        for s in sectors:
            if s.start_res==n:
                sorted_sectors.append(s)
                sectors.remove(s)

    return sorted_sectors

def array_frequency(a):
    #input is numpy array
    #output is list of 
    unique, inverse = numpy.unique(a, return_inverse=True)
    count=numpy.zeros(len(unique), numpy.int)
    numpy.add.at(count, inverse, 1)
    return numpy.vstack((unique, count)).T

def get_residue_rate_probabilities(modelfile, scorefile, sectors, seq, grid, num_models=5, outfile="rate_probabilities.dat", offset=0):
    # Given a set of models (from a best_models.dat file)
    # returns the probability of observing each rate
    # If the grid is given, it is outputted in the first line
    # sectors can be a list of Sector objects, or tuples (first_res, last_res)

    of=open(outfile, "w")

    if hasattr(grid, '__iter__'):
        grid=len(grid)
    
    best_models, best_scores=get_best_scoring_models(modelfile, scorefile, num_models, write_file=False)

    best_models=numpy.array(best_models)

    sorted_sectors=sector_sort(sectors)

    # Loop over all sectors and add up instances of each rate bin
    for s in sorted_sectors:
        # get all instances of each rate
        freq=array_frequency(best_models[:, s.start_res:s.end_res+1])
        model=numpy.zeros(grid)
        for i in freq:
            model[i[0]]=1.0*i[1]/s.num_amides/num_models

        for n in range(s.start_res, s.end_res+1):
            if seq[n+offset]=="P":
                of.write(n, "P", numpy.zeros(len(grid), numpy.int))
            else:
                of.write(str(n+1)+", "+seq[n+offset]+", " +str([m for m in model])+"\n")


def get_convergence(state, num_points=None):
    """ Takes all models in the exp_model of the given state and 
    divides them into two halves.  The average and SD for each residue is computed
    for each ensemble and compared via a Two-sample Kolmogorov-Smirnov Test"""
    these_states=model.states
    for state in these_states:
        es=state.exp_model
        sectors=state.sectors
        if num_points is None or num_points > len(es.exp_models)/2:
            num_points=len(es.exp_models)/2
        exp_model1=es.exp_models[len(es.exp_models)/2-num_points:len(es.exp_models)/2]
        exp_model2=es.exp_models[len(es.exp_models)-num_points:-1]
        zscore=calculate_zscore(exp_model1, exp_model2, state, state)
        #print(state.state_name)
        #print zscore
        #time.sleep(1)

def get_average_sector_values(exp_models, state):
    sector_model = get_sector_averaged_models(exp_models, state)
    avg = numpy.average(sector_model,0)
    std = numpy.std(sector_model,0)
    return (avg, std)

def get_cdf(exp_models):
    """ Takes a list of 1D exp_models and returns a sorted numpy array
    equivalent to the empirical density function for each residue. """
    exp_model_edf=numpy.empty((len(exp_models),len(exp_models[0])))
    A=numpy.array(exp_models)
    y=numpy.linspace(1./len(exp_models),1,len(exp_models))
    print(len(exp_models[0]))
    for i in range(len(exp_models[0])): 
        counts, edges = numpy.histogram(A[:,i], len(A), range=(-6,0), density=False) 
        #print i,A[:,i],counts,numpy.cumsum(counts*1.0/len(A)) 
        exp_model_edf[:,i]=numpy.cumsum(counts*1.0/len(A))
    return exp_model_edf

def get_chisq(exp_models1, exp_models2, nbins):
    """ Takes two lists of exp_models and returns the chi2 value along the second axis """
    A=numpy.array(exp_models1)
    B=numpy.array(exp_models2)
    #y=numpy.linspace(1./len(exp_models1),1,len(exp_models1))
    print(len(exp_models1[0]))
    for i in range(269,len(exp_models1[0])):
        meanA = numpy.mean(A[:,i])
        ssd = numpy.std(A[:,i])**2 + numpy.std(B[:,i])**2 
        sstdev = numpy.sqrt( ssd / 5000 )
        meanB = numpy.mean(B[:,i])
        t = 1.96
        ci = t * sstdev
        dm = meanA - meanB
        print(i, dm, ci, dm/ci)
        #fig=plt.figure()
        #ax1 = fig.add_subplot(111)
        #ax1.plot(range(nbins), countsA)
        #ax1.plot(range(nbins), countsB)
        #plt.show()

        #data = [countsA, countsB]
        #print(i, chi2_contingency(data))
    return exp_model_edf


def calculate_ks_statistic(edf1, edf2):
    """ Takes a two edfs and returns a vector of the Kolmogorov-Smirnov 
    statistic for each residue"""
    maxdiff=numpy.zeros(len(edf1[0]))
    threshold=1.98*numpy.sqrt(1.0*(len(edf1)+len(edf2))/(1.0*len(edf1)*len(edf2)))
    if len(edf1[0]) != len(edf2[0]):
        print("Different Number of Residues for EDFs in KS calculation: Exiting")
        exit()
    for r in range(len(edf1[0])):
        maxdiff[r]=0
        for m in range(len(edf1[:,0])):
            diff=abs(edf1[m,r]-edf2[m,r])
            if diff > maxdiff[r]:
                maxdiff[r]=diff
    return maxdiff, threshold

def get_sector_averaged_models(exp_models, state):
    sector_avg_models=[]
    for n in range(len(exp_models)):
        sector_avg_models.append(state.exp_model.get_sector_averaged_protection_values(exp_models[n], state.sectors))
    return sector_avg_models

def calculate_zscore(exp_models1, exp_models2, state1, state2):
    avg1, sd1=get_average_sector_values(exp_models1, state1)
    avg2, sd2=get_average_sector_values(exp_models2, state2)
    zscore=numpy.subtract(avg1,avg2)/numpy.sqrt(numpy.add(numpy.add(numpy.square(sd1),numpy.square(sd2)),0.00001))
    return zscore

def calculate_convergence(exp_models1, exp_models2):
    return 0
