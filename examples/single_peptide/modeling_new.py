## \example hdx/ERa/ERa_088074_modeling.py
import sys
sys.path.append( "../../../bayesian_hdx/pyext/src" )
import system_setup
import hdx_models
import sampling
from numpy import linspace
import numpy
import scipy

##### File/Directory Setup ######
outputdir="./output"
datadir="./"
workbench_file=datadir + "single_frag.csv"

####### Input Parameters #######
offset=0        # offset between fragment start/end values and FASTA sequence
sigma0=50.0        # Estimate for slope of experimental error (slope of sd vs. %D plot)
num_exp_bins=7   # Number of bins

exp_grid=linspace(-6,0,num_exp_bins)

inseq="YSMKSKNVVPLYDLL"

data  = {10   : [21.33928571428864, 15.702380952384537,  6.877976190478665],
         30   : [21.101190476192265,12.139880952382162, 23.738095238098282],
         60   : [29.78869047619391, 27.208333333334846, 17.913690476193878],
         300  : [32.130952380956586,20.264880952382892, 31.663690476193555],
         900  : [38.83035714286129, 30.360119047621147, 49.33035714286134],
         3600 : [51.97321428571825, 37.898809523810655, 41.9226190476218], }


data  = {10   : [6.31378],
         30   : [10.83015],
         60   : [13.55446],
         300  : [22.86113],
         900  : [31.88971],
         3600 : [42.51679] }


def truncated_gaussian_factor(exp,a,b,sig):
    print scipy.special.erf( (b-exp)/sig), scipy.special.erf( (a-exp)/sig ), b-exp, a-exp
    return 1/ ( 0.5 * ( scipy.special.erf( (b-exp)/sig * numpy.sqrt(3.1415) ) - scipy.special.erf( (a-exp)/sig * numpy.sqrt(3.1415) ) )  )

###############################
###   System Setup:
###############################


# Initialize model  (name, FASTA sequence, offset)
model=system_setup.HDXModel("simple",
                                    inseq,
                                    offset=offset)

# Create a fragment 
frag = system_setup.Fragment(inseq, 1, 15, 0)

# Add data to the fragment
for t in data:
    tp=frag.add_timepoint(t)
    for d in data[t]:
        tp.add_replicate(d)

# enumerate and score all solutions

best_grid = sampling.enumerate_fragment(frag, exp_grid, sigma0)

###############################
###   Analysis:
###############################

grid=[0, 3, 2, 3, 1, 1,2]
score=0
for t in frag.timepoints:
    model = t.calc_model_deut(grid, exp_grid, frag.num_observable_amides)
    diff = model-t.replicates[0].deut
    #likelihood = 
    this_tp_score=-1.0*numpy.log(t.replicates[0].gaussian_model(t.replicates[0].deut, model, sigma0))
    tg_factor = truncated_gaussian_factor(t.replicates[0].deut, -10, 120, sigma0)
    score=score + this_tp_score * tg_factor
    stuff = numpy.exp(-(numpy.log(model)-numpy.log(t.replicates[0].deut))**2/(2*sigma0**2))
    print t.time, model, t.replicates[0].deut, diff, stuff, tg_factor, this_tp_score, score
    #print score

'''
score=0
for t in frag.timepoints:
    model=t.calc_model_deut(best_grid, exp_grid, frag.num_observable_amides)
    score=score-1.0*numpy.log(t.replicates[0].lognormal_model(t.replicates[0].deut, model, sigma0))
    print t.time, model-t.replicates[0].deut, numpy.exp(-(numpy.log(model)-numpy.log(t.replicates[0].deut))**2/(2*sigma0**2)), 1/(model*numpy.sqrt(2*3.14)*sigma0),-1.0*numpy.log(t.replicates[0].lognormal_model(t.replicates[0].deut, model, sigma0))
    print score
'''
