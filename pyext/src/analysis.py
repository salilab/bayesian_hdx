"""@namespace IMP.hdx.analysis
   Analysis functions for HDX simulations
"""
from __future__ import print_function
import hdx_models
import input_data
import plots
from scipy.stats import cumfreq
from scipy.stats import chi2_contingency
import numpy
import time
import os.path
from pylab import *
from matplotlib import *


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
        return outfile
    else:
        for i in top_score_indices:
            top_models.append(map(int, models[int(i)].split()) )
        return top_models, top_scores

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
        print(state.state_name)
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
