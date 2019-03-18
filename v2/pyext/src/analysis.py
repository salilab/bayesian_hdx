"""@namespace IMP.hdx.analysis
   Analysis functions for HDX simulations
"""
from __future__ import print_function
import hdx_models
import input_data
import hxio
import plots
#from scipy.stats import cumfreq
#from scipy.stats import chi2_contingency
import numpy
import math
import time
import os.path
from pylab import *
from matplotlib import *

class Precision(object):
    '''
    Given a set of ParseOutputFile classes that have already undergone convergence testing, 
    calculate the precision of each cluster (and accuracy, if exact info available)
    '''
    def __init__(self, parse_output1, parse_output2, gsm="all"):
        self.sample1 = parse_output1
        self.sample2 = parse_output2


class Convergence(object):
    '''
    Given two ParseOutputFile classes, 
    allow testing of various convergence and clustering metrics
    '''
    def __init__(self, parse_output1, parse_output2, gsm="all"):
        self.sample1 = parse_output1
        self.sample2 = parse_output2

        self.models1 = self.sample1.get_all_models()
        self.models2 = self.sample2.get_all_models()

        self.num_gsm = gsm

    def calc_num_of_gsm(self, models):
        # For each set of models (models1, models2), calculate the optimal number
        # of good scoring models.

        for i in range(3,len(models)):
            j=0

    def get_scores(self, num_gsm="all"):
        #gsm1, gsm2 = self.get_gsm(num_gsm)
        return [g[0] for g in self.models1], [g[0] for g in self.models2]

    def get_models(self, num_gsm="all"):
        m1 = sorted(self.models1, key=lambda x: x[0])
        m2 = sorted(self.models1, key=lambda x: x[0])
        if num_gsm=="all":
            mod1 = [g[1] for g in m1]
            mod2 = [g[1] for g in m2]
        else:
            mod1 = [g[1] for g in m1[0:num_gsm]]
            mod2 = [g[1] for g in m2[0:num_gsm]]            
        return mod1, mod2         

    def total_score_pvalue_and_cohensd(self, num_gsm="all"):
        try:
            from scipy import stats
        except:
            print("Convergence:: scipy not installed. Skipping p-value test")
            return -1

        ts1, ts2 = self.get_scores(num_gsm)

        st, pvalue = stats.ttest_ind(ts1, ts2, equal_var=False)

        mean_diff = numpy.mean(ts1) - numpy.mean(ts2)

        vari1 = numpy.var(ts1)
        vari2 = numpy.var(ts2)

        #print(numpy.array(ts1).mean(), numpy.array(ts2).mean(), vari1, vari2, st)
        cohens_d = mean_diff / numpy.sqrt(len(ts1)*vari1 + len(ts2)*vari2) / (len(ts1)+len(ts2)-2)

        return pvalue, cohens_d

    def residue_pvalue_and_cohensd(self, num_gsm="all"):

        try:
            from scipy import stats
        except:
            print("Convergence:: scipy not installed. Skipping p-value test")
            return -1
        
        bsm1, bsm2 = self.get_models(num_gsm)
        output=[]

        # Loop over all residues
        for i in range(len(bsm1[0])):
            if bsm1[0][i] != 0:
                # now average over all
                b1 = []
                b2 = []
                for s in range(len(bsm1)):
                    b1.append(bsm1[s][i])
                for s in range(len(bsm2)):
                    b2.append(bsm2[s][i])

                #print(b1)
                #print(b2)

                st, pvalue = stats.ttest_ind(b1, b2, equal_var=False)

                mean_diff = numpy.mean(b1)- numpy.mean(b2)
                vari1 = numpy.var(b1)
                vari2 = numpy.var(b2)
                cohens_d = mean_diff / numpy.sqrt(len(b1)*vari1 + len(b2)*vari2) / (len(b1)+len(b2)-2)
                #print(numpy.array(b1).mean(), numpy.array(b2).mean(), vari1, vari2)
                output.append((pvalue, cohens_d))
            else:
                output.append((0,0))

        return output

    def get_distance_matrix(self, num_models='all'):
        from sklearn.metrics.pairwise import pairwise_distances
        if num_models == 'all':
            num_models = self.num_gsm

        bsm1, bsm2 = self.get_models(self.num_gsm)
        num_models = len(bsm1) + len(bsm2)

        print(len(bsm1+bsm2), self.num_gsm)
        #print(len(bsm1.join(bsm2)))

        self.full_distance_matrix = (pairwise_distances(numpy.array(bsm1+bsm2), metric='cosine', n_jobs=-1)* 10)**2
        return self.full_distance_matrix

    def precision_cluster(self, threshold):

        distmat = self.get_distance_matrix(num_models='all')

        # 1) get neighbors
        neighbors=[]
        for count in range(len(distamat[0])):
            neighbors.append([count]) # model is a neighbor of itself

        for i in range(num_models-1):
            for j in range(i+1,num_models):    
                if distmat[i][j]<=rmsd_cutoff:
                    neighbors[i].append(j)
                    neighbors[j].append(i)

        # 2). Get the weightiest cluster, and iterate
        unclustered=[]
        boolUnclustered=[]
        for i in range(num_models):
            unclustered.append(i)
            boolUnclustered.append(True)

        cluster_members=[] # list of lists : one list per cluster
        cluster_centers=[]


        while len(unclustered)>0:
            # get cluster with maximum weight
            max_neighbors=0
            currcenter=-1
            for eachu in unclustered:  # if multiple clusters have same maxweight this tie is broken arbitrarily! 
                if len(neighbors[eachu])>max_neighbors:
                    max_neighbors=len(neighbors[eachu])
                    currcenter=eachu   
       
            #form a new cluster with u and its neighbors
            cluster_centers.append(currcenter)
            cluster_members.append([n for n in neighbors[currcenter]]) 

            #update neighbors 
            for n in neighbors[currcenter]:
                #removes the neighbor from the pool
                unclustered.remove(n) #first occurence of n is removed. 
                boolUnclustered[n]=False # clustered

            for n in neighbors[currcenter]:
                for unn in neighbors[n]: #unclustered neighbor
                    if not boolUnclustered[unn]:
                        continue
                    neighbors[unn].remove(n)
        
        return cluster_centers, cluster_members

    def get_clusters(self, cutoffs_list):
        #Do Clustering on a Grid
        pvals=[]
        cvs=[]
        percents=[]

        f1=open("%s.ChiSquare_Grid_Stats.txt" % "Test", 'w+')

        bsm1, bsm2 = self.get_models(self.num_gsm)

        for c in cutoffs_list:
            cluster_centers,cluster_members = self.precision_cluster(c)

            ctable,retained_clusters = self.get_contingency_table(len(cluster_centers), cluster_members, bsm1+bsm2,
                                                           bsm1,bsm2)

            (pval,cramersv) = self.test_sampling_convergence(ctable, len(bsm1)+len(bsm2))
            percent_explained = self.percent_ensemble_explained(ctable, len(bsm1)+len(bsm2))

            pvals.append(pval)
            cvs.append(cramersv)
            percents.append(percent_explained)
            
            print >>f1, c, pval, cramersv, percent_explained

        return pvals, cvs, percents

    def get_cutoffs_list(self, gridSize):

        mindist = self.full_distance_matrix.min()
        maxdist = self.full_distance_matrix.max()

        print("Minimum and maximum pairwise model distances:",mindist, maxdist)
        cutoffs = numpy.arange(mindist,maxdist,gridSize)
        return cutoffs

    def get_contingency_table(self, num_clusters,cluster_members,all_models,run1_models,run2_models):
        full_ctable=numpy.zeros((num_clusters,2))

        for ic,cluster in enumerate(cluster_members):
            for member in cluster:
                    model_index=all_models[member]

                    if model_index in run1_models:
                        #print "run1", model_index
                        full_ctable[ic][0]+=1.0
                    elif model_index in run2_models:
                    #print "run2", model_index
                        full_ctable[ic][1]+=1.0

        ## now normalize by number of models in each run
        numModelsRun1 = float(numpy.sum(full_ctable,axis=0)[0])
        numModelsRun2 = float(numpy.sum(full_ctable,axis=0)[1])

        reduced_ctable=[]
        retained_clusters=[]
        
        for i in range(num_clusters):
            if full_ctable[i][0]<=10.0 or full_ctable[i][1]<=10.0:
                #if full_ctable[i][0]<=0.10*numModelsRun1 and full_ctable[i][1] <= 0.10*numModelsRun2:
                continue
            reduced_ctable.append([full_ctable[i][0],full_ctable[i][1]])
            retained_clusters.append(i)
        return numpy.array(reduced_ctable),retained_clusters


    def get_sampling_precision(self, cutoffs_list, pvals, cvs, percents):
        sampling_precision=max(cutoffs_list)
        pval_converged=0.0
        cramersv_converged=1.0
        percent_converged=0.0

        for i in range(len(cutoffs_list)):
            if percents[i]>80.0:
                if pvals[i]>0.05 or cvs[i]<0.10:
                    if sampling_precision>cutoffs_list[i]:
                        sampling_precision=cutoffs_list[i]
                        pval_converged=pvals[i]
                        cramersv_converged=cvs[i]
                        percent_converged=percents[i]
            else:
                sampling_precision=max(cutoffs_list)

        return sampling_precision,pval_converged,cramersv_converged,percent_converged

class ParseOutputFile(object):
    def __init__(self, output_file, state_name):
        self.output_file = output_file
        self.state_name = state_name
        self.datafiles = []
        self.sectors = []
        self.pf_grids = {}
        self.observed_residues = []
        self.parse_header()
        self.path = os.path.dirname(os.path.realpath(output_file))

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
            
            # #-symbol means datasets
            elif line[0:2]=="# ":
                self.datafiles.append( (line[2:].split("|")[0].strip(), float(line[2:].split("|")[2].strip())) )
            
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

    def cluster_models_kmeans(self, nmodels, nclust):
        # Uses kmeans to cluster models
        from sklearn.cluster import KMeans
        mods = [m[1] for m in self.get_best_scoring_models(N=nmodels)]
        models = numpy.array(mods)
        kmeans = KMeans(n_clusters=nclust).fit(models)

        #print(kmeans.labels_)

        for i in range(nclust):
            unique, counts = numpy.unique(kmeans.labels_, return_counts=True)
            print(nclust, " | ", i, dict(zip(unique, counts))[i] *1.0/nmodels)

        for c in range(len(kmeans.cluster_centers_)):
            for i in range(c,len(kmeans.cluster_centers_)):
                print(i, c, numpy.linalg.norm(kmeans.cluster_centers_[c]-kmeans.cluster_centers_[i]))
        print("###")


    def get_datasets(self):
        if len(self.datasets) == 0:
            self.generate_datasets()
        return self.datasets

    def generate_datasets(self):
        self.datasets=[]
        for f in self.datafiles:
            print(self.path, f[0])
            self.datasets.append(hxio.import_json(self.path+f[0]))

    def get_all_models(self, return_pf=False):
        f = open(self.output_file, "r")
        models = []
        # Cycle over all lines
        for line in f.readlines():       
            if line[0]==">":
                score = float(line.split("|")[1].strip())

                model_string = line[1:].split("|")[0].strip()
                model_list = []
                for m in model_string.split(" "):
                    model_list.append(int(m))
                if return_pf:
                    ml1 = model_list
                    model_list = self.models_to_protection_factors(model_list)
                    #print(score, model_list, ml1)
                models.append((score, model_list)) 

        self.models = models

        return models

    def get_scores(self):
        try:
            models = self.models
        except:
            models = self.get_all_models()

        scores = [mod[0] for mod in models]

        return numpy.sort(scores)


    def get_best_scoring_models(self, N="all", sigmas=False, return_pf=False, sort_sectors=False):
        ''' Get the N best scoring models from the output file
        Returns a list of tuples of best_scoring_models 
            [(score, [model])]
        and (if sigmas=True)
        a grid of the timepoint sigma values.
        '''
        # Model entries are marked with a > as the first character

        if N=="all":
            N=100000000000

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
                        #print(score, best_scoring_models[-1])
                    model_string = line[1:].split("|")[0].strip()
                    model_list = []
                    for m in model_string.split(" "):
                        model_list.append(int(m))
                    if sort_sectors:
                        model_list = self.sort_model_by_sector(model_list)
                    if return_pf:
                        ml1 = model_list
                        model_list = self.models_to_protection_factors(model_list)
                        #print(score, model_list, ml1)
                    best_scoring_models.append((score, model_list))
                    best_scoring_models = sorted(best_scoring_models, key=lambda x: x[0])

        self.best_scoring_models = best_scoring_models
        return best_scoring_models

    def models_to_protection_factors(self, models):
        # Input a list of list of integers.  
        # CHECK THAT THE MODEL SIZE IS CORRECT!
        if type(models[0]) != list:
            models = [models]
        #print(models)
        output = []
        for m in models:
            pf_model = []
            for res in range(len(m)):
                if res + 1 in self.observed_residues:
                    #print("RES", res+1, m[res], self.pf_grids[res+1][m[res]-1])
                    pf_model.append(float(self.pf_grids[res+1][m[res]-1]))

                else:
                    pf_model.append(numpy.nan)
            output.append(pf_model)

        return output

    def get_sectors(self):
        # returns the list of sectors
        return self.sectors

    def sort_model_by_sector(self, model):
        # Given an input of a single model (as a list of integers), and a list
        # of sectors (as a list of list of residue numbers), return the model
        # with the integers in each sector sorted in increasing order. 

        # check that the model and total length of the sectors is the same:
        s_len = 0
        s_len = sum([len(s) for s in self.sectors])
        #if s_len != len(model):
        #    print("ERROR: sector length", s_len, "is not equal to model length ", len(model))
        #    exit()
        out_model = numpy.zeros(len(model))
        # Loop over all sectors
        for s in self.sectors:
            if len(s) == 0:
                continue
            #print(s)
            sec_model = model[s[0]-1:s[-1]]
            sort_model = numpy.sort(sec_model) #sort indexes in increasing order

            # find indexes = 0
            idx_zero = [index for index, v in enumerate(sort_model) if v == 0]
            idx_nonzero = [index for index, v in enumerate(sort_model) if v != 0]
            for i in range(len(idx_nonzero)):
                out_model[idx_nonzero[i] - 1 + s[0]] = sort_model[i + len(idx_zero)]
            #print(s, sec_model, sort_model, out_model[s[0]-1:s[-1]])
        return out_model.astype(int)

    def sort_models(self):
        pass


    def generate_pf_probability_file(self, num_best_scoring_nmodels):
        bsm = self.get_best_scoring_models(num_best_scoring_models, return_pf=True)
        return 0

    def calculate_random_sample_convergence(self, replicates=100, pct_values=10):
        """ 
            For each sample of good scoring models (self.get_best_scoring_models)

            Get a list of top scores from each.

            Store the results as a list of two lists (one for each sample)
            with each sample list containing a tuple :: (pct, avg, stdev)
        """
        delpct = 1.0 / pct_values

        pct_grid = numpy.arange(delpct, 1.0, delpct)

        convergence_tuples = []

        scores = [g[0] for g in self.get_best_scoring_models()]

        pct_tuples = []

        for p in pct_grid:
            min_scores = []
            # For each percent value, take N samples and report avg/std of minimum values
            for r in range(replicates):
                subset = self._random_subset(scores, p)
                min_scores.append(min(subset))
            array = numpy.array(min_scores)

            pct_tuples.append((p, numpy.average(array), numpy.std(array)))

            print((p, numpy.average(array), numpy.std(array)), len(subset))

        convergence_tuples.append(pct_tuples)

        #print(convergence_tuples)
        del scores


    def _random_subset(self, models, pct):
        """
            Generate a random subset of pct percent of the given list of things
            @param models - python list of model dictionaries (really, could be anything)
            @param pct - the fraction of samples to return
        """
        import random
        if pct > 1.0:
            print("WARNING: percent value should be between 0 and 1.0. Diving by 100")
            pct = pct/100.0
        if pct > 1.0 or pct <= 0:
            print("WARNING: percent value is outside of bounds")

        num_samples = int(len(models) * pct)
        rand_smpl = [models[i] for i in random.sample(xrange(len(models)), num_samples)]
        #print(num_samples, len(models), len(rand_smpl), models)
        return rand_smpl

class OutputAnalysis(object):
    '''
    Class that:
        * Analyzes an output file(s)
        * Produces set of standard graphs
        * Outputs a file with Pf probability distributions
    '''
    def __init__(self, output):
        self.output = output


    def parse_output_files(self):
        pass


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
    #output is list of tuples with (value, frequency)
    unique, inverse = numpy.unique(a, return_inverse=True)
    count = numpy.zeros(len(unique), numpy.int)
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
