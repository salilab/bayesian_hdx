"""@namespace IMP.hdx.analysis
   Analysis functions for HDX simulations
"""
from __future__ import print_function
import hxio
import tools
#from scipy.stats import cumfreq
#from scipy.stats import chi2_contingency
import numpy
import math
import time
import scipy
from scipy import stats
from copy import deepcopy
import os.path
#from pylab import *
#from matplotlib import *
import sklearn
from sklearn.metrics.pairwise import pairwise_distances

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

        self.models1, self.sigmas1 = self.sample1.get_all_models()
        self.models2, self.sigmas2 = self.sample2.get_all_models()

        self.num_gsm = gsm

    def calc_num_of_gsm(self, models):
        # For each set of models (models1, models2), calculate the optimal number
        # of good scoring models.
        # NOT IMPLEMENTED!!
        for i in range(3,len(models)):
            j=0

    def get_scores(self, num_gsm="all"):
        #gsm1, gsm2 = self.get_gsm(num_gsm)
        return [g[0] for g in self.models1], [g[0] for g in self.models2]

    def get_models(self, num_gsm="all"):
        m1 = sorted(self.models1, key=lambda x: x[0])
        m2 = sorted(self.models2, key=lambda x: x[0])
        if num_gsm=="all":
            mod1 = [(g[0], g[1]) for g in m1]
            mod2 = [(g[0], g[1]) for g in m2]
        else:
            mod1 = [(g[0], g[1]) for g in m1[0:int(num_gsm)]]
            mod2 = [(g[0], g[1]) for g in m2[0:int(num_gsm)]]            
        return mod1, mod2         

    def total_score_pvalue_and_cohensd(self, num_gsm="all"):

        ts1, ts2 = self.get_scores(num_gsm)

        st, pvalue = stats.ttest_ind(ts1, ts2, equal_var=False)

        mean_diff = numpy.mean(ts1) - numpy.mean(ts2)

        vari1 = numpy.var(ts1)
        vari2 = numpy.var(ts2)

        #print(numpy.array(ts1).mean(), numpy.array(ts2).mean(), vari1, vari2, st)
        cohens_d = mean_diff / numpy.sqrt(len(ts1)*vari1 + len(ts2)*vari2) / (len(ts1)+len(ts2)-2)

        return pvalue, cohens_d

    def residue_pvalue_and_cohensd(self, num_gsm="all"):
        
        bsm1, bsm2 = self.get_models(num_gsm)
        output=[]

        # Loop over all residues
        for i in range(len(bsm1[0][1])):
            if bsm1[0][1][i] != 0:
                # now average over all
                b1 = []
                b2 = []
                for s in range(len(bsm1)):
                    b1.append(bsm1[s][1][i])
                for s in range(len(bsm2)):
                    b2.append(bsm2[s][1][i])

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

    def get_distance_matrix(self, num_models='all', njobs=-3):
        if num_models == 'all':
            num_models = self.num_gsm

        self.num_gsm = num_models

        self.bsm1, self.bsm2 = self.get_models(num_models)

        mod1 = [m[1] for m in self.bsm1]
        mod2 = [m[1] for m in self.bsm2]
        self.full_distance_matrix = (pairwise_distances(numpy.array(mod1+mod2), metric='cosine', n_jobs=njobs)* 10)**2
        return self.full_distance_matrix

    def precision_cluster(self, threshold):

        distmat = self.full_distance_matrix  #self.get_distance_matrix(num_models='all')
        num_models = len(distmat[0])

        # 1) get neighbors
        neighbors=[]
        for count in range(num_models):
            neighbors.append([count]) # model is a neighbor of itself

        for i in range(num_models-1):
            for j in range(i+1,num_models): 
                if distmat[i][j]<=threshold:
                    neighbors[i].append(j)
                    neighbors[j].append(i)
                #print(i,j,distmat[i][j],len(neighbors[i]),len(neighbors[j]))   


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
        # Do Clustering on a Grid
        pvals=[]
        cvs=[]
        percents=[]

        f1=open("%s.ChiSquare_Grid_Stats.txt" % "Test", 'w')
        f1.write("Threshold, Pvalue, CramersV, Pct_Ensemble_Explained\n")
        bsms1, bsms2 = self.get_models(self.num_gsm)

        bsm1 = [m[1] for m in bsms1]
        bsm2 = [m[1] for m in bsms2]

        for c in cutoffs_list:
            cluster_centers,cluster_members = self.precision_cluster(c)

            ctable,retained_clusters = self.get_contingency_table(len(cluster_centers), cluster_members, bsm1+bsm2,
                                                           bsm1,bsm2)

            (pval,cramersv) = self.test_sampling_convergence(ctable, len(bsm1)+len(bsm2))
            percent_explained = self.percent_ensemble_explained(ctable, len(bsm1)+len(bsm2))

            pvals.append(pval)
            cvs.append(cramersv)
            percents.append(percent_explained)
            
            f1.write(str(c)+", "+str(pval)+", "+str(cramersv)+", "+str(percent_explained)+"\n")

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
                        full_ctable[ic][0]+=1.0
                    elif model_index in run2_models:
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


    def test_sampling_convergence(self, contingency_table, total_num_models):

        if len(contingency_table)==0:
            return 0.0,1.0
        
        ct = numpy.transpose(contingency_table)
        [chisquare,pvalue,dof,expected]=scipy.stats.chi2_contingency(ct)
        if dof==0.0:
            cramersv=0.0
        else:
            cramersv=math.sqrt(chisquare/float(total_num_models))
            
        return(pvalue,cramersv)

    def percent_ensemble_explained(self, ctable,total_num_models):
        if len(ctable)==0:
            return 0.0
        percent_clustered=float(numpy.sum(ctable,axis=0)[0]+numpy.sum(ctable,axis=0)[1])*100.0/float(total_num_models)
        return percent_clustered

    def get_sampling_precision(self, cutoffs_list, pvals, cvs, percents):
        sampling_precision=max(cutoffs_list)
        pval_converged=0.0
        cramersv_converged=1.0
        percent_converged=0.0

        self.sampling_precision_stats = {
        "pvals" : [],
        "cvs" : [],
        "pcts" : []
        }

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

        self.sampling_precision = sampling_precision

        return sampling_precision,pval_converged,cramersv_converged,percent_converged

    def cluster_at_threshold_and_return_pofs(self, threshold):
        '''
        Cluster the gsms at a given threshold and return lists of POF objects 
        corresponding to all models from each cluster
        '''
        cluster_centers,cluster_members = self.precision_cluster(threshold)

        all_models = self.bsm1 + self.bsm2

        cluster_pofs = []
        print("Cluster Sizes:")
        print(" -- Cluster# number_of_models")
        for c in range(len(cluster_centers)):
            # Create a copy of the PO file
            new_po = self.sample1.clone_pof()
            #ParseOutputFile(self.sample1.output_file)
            new_po.clear_models()
            # Loop over all models in the cluster
            for m in cluster_members[c]:
                new_po.models.append(all_models[m])
            print(" --", c, len(new_po.models))
            cluster_pofs.append(new_po)

        return cluster_pofs


class DeltaHDX(object):
    '''
    A class that utilizes two POF objects to compute the mean (pof2 - pof1) and standard deviation 
    between the protection factors at each residue
    '''
    def __init__(self, pof1, pof2):
        self.pof1 = pof1
        self.pof2 = pof2

    def calculate_dhdx(self, weighted=False):
        # Calculates delta HDX
        #pf_models1 = numpy.array([m[1][0] for m in self.pof1.get_models(return_pf=True)])
        #pf_models2 = numpy.array([m[1][0] for m in self.pof2.get_models(return_pf=True)])
        #print(len(self.pof1.get_all_models(return_pf=True)), len(self.pof1.get_all_models(return_pf=True)[0]), self.pof1.get_all_models(return_pf=True)[0][0])
        pf_models1 = numpy.array([m[1][0] for m in self.pof1.get_all_models(return_pf=True)[0]])
        scores1 = numpy.array([m[0] for m in self.pof1.get_all_models(return_pf=True)])
        pf_models2 = numpy.array([m[1][0] for m in self.pof2.get_all_models(return_pf=True)[0]])
        scores2 = numpy.array([m[0] for m in self.pof2.get_all_models(return_pf=True)])

        if weighted:
            weights1 = numpy.exp(-1*(scores1-min(scores1)))
            mean1 = numpy.average(pf_models1, axis=0, weights=weights1)
            mean2 = numpy.average(pf_models1, axis=0, weights=weights2)
        else:
            mean1 = numpy.mean(pf_models1, axis=0)
            mean2 = numpy.mean(pf_models2, axis=0)

        std1 = numpy.std(pf_models1, axis=0)
        std2 = numpy.std(pf_models2, axis=0)

        diff = mean2 - mean1
        Z = (mean2 - mean1)/numpy.sqrt(numpy.square(std1)+numpy.square(std2)+0.0000001)

        return diff, Z, mean1, mean2, std1, std2

    def write_dhdx_file(self, prefix=None, resrange=None):
        '''
        Calculate the dhdx file and write a file of format:
        Res# Res state_name dhdx dhdx_z avg_lig avg_apo sd_lig sd_apo flag
        '''
        state_name1 = self.pof1.state_name
        state_name2 = self.pof2.state_name
        molecule_name = self.pof1.molecule_name
        sequence = self.pof1.get_sequence()

        filename = molecule_name+"_"+state_name1 + "_"+state_name2+".dhdx"

        if prefix is not None:
            filename = prefix + filename

        f = open(filename, "w")
        f.write("Res# Res XX dhdx dhdx_z avg_state2 avg_state1 sd_state2 sd_state1 flag\n")

        if resrange is None:
            minres = min(min(self.pof1.observed_residues), min(self.pof2.observed_residues))
            maxres = max(max(self.pof1.observed_residues), max(self.pof2.observed_residues))

            resrange = range(minres, maxres)
        else:
            resrange = resrange

        diff, Z, mean1, mean2, sd1, sd2 = self.calculate_dhdx()

        for r in resrange:
            f.write(str(r)+" "+str(sequence[r-1])+" XX %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f \n" % (diff[r-1], Z[r-1], mean2[r-1], mean1[r-1], sd2[r-1], sd1[r-1]))

        f.close()

    def write_pymol_file(self, pdb_file, pymol_file, pdb_offset=0):
        '''
        Write a pymol file that incorporates DHDX information onto a crystal structure
        '''
        
        pass

    def write_volcano_plot(self, plot_file):
        '''
        write out a volcano plot for all residues
        '''
        pass


class ParseOutputFile(object):
    '''
    An object that stores all of the information from an hxio.Output file.
    Also can be used as a general tool for analyzing sets of data.

    To create a copy of this object, do >new_pof = deepcopy(pof) to retain all
    of the header information and then clear the models
    '''
    def __init__(self, output_file, all_observed=True):
        self.header = {}
        self.output_file = output_file
        self.datafiles = []
        self.datasets = []
        self.sectors = []
        self.models=[]
        self.pf_grids = {}
        self.parse_header()
        if all_observed:
            self.observed_residues = list(set([item for sublist in self.sectors for item in sublist]))
        else:
            self.observed_residues = range(1,len(self.get_sequence()+1))
        self.path = os.path.dirname(os.path.realpath(output_file))
        self.generate_datasets()

    def clone_pof(self):
        new_pof = ParseOutputFile(self.output_file)
        new_pof.datafiles = self.datafiles
        new_pof.sequence = self.sequence
        new_pof.grid_size = self.grid_size
        new_pof.molecule_name = self.molecule_name
        new_pof.pf_grids = self.pf_grids
        new_pof.state_name = self.state_name
        new_pof.sectors = self.sectors
        new_pof.models = self.models
        return new_pof

    def get_sequence(self):
        '''
        Return the macromolecule sequence from the output file
        '''
        if len(self.get_datasets()) > 0:
            return self.get_datasets()[0].sequence
        else:
            return self.seq

    def compute_residue_resolved_moves(self, magnitude=False, model_range=(0,-1)):
        '''
        For this POF, compute the number of moves for each residue over all models

        This function loops over all models.

        A move for a single residue is logged when the HDX value is different
        '''
        moves = numpy.zeros((len(self.models[0][1]), len(self.pf_grids[list(self.pf_grids.keys())[0]])))

        modm1 = self.models[0][1]
        # Loop over all models
        for m in range(1,len(self.models)):
            if m%10000==0:
                print("Model", m, "of",len(self.models))
            mod = self.models[m][1]
            # Loop over all residues
            for r in range(len(mod)):
                mod_diff = abs(mod[r]-modm1[r])
                moves[r][mod_diff]+=1
            modm1 = mod

        if magnitude:
            return moves
        else:
            return numpy.sum(moves[:][1:], axis=1)


    def get_all_peptides(self):
        '''
        Get all the peptides from all datasets
        '''
        peptides = []
        for d in self.datasets:
            peptides += d.get_peptides()
        return peptides

    def parse_header(self):
        '''
        Function that parses the header of output files.
        Stores the logk grid, along with other experimental information
        '''
        f = open(self.output_file, "r")
        for line in f.readlines():
            
            # > means model data (so header is over.)
            if line[0]==">":
                break

            elif line.split(":")[0].strip() == "Sequence":
                self.sequence = line.split(":")[1].strip()
                self.seq = self.sequence
            
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
                            #self.observed_residues.append(int(r))
                    self.sectors.append(sector)

            elif line[0:2]=="$ ":
                if line[2:].split("|")[0].strip() != "Residue_number":
                    res = int(line[2:].split("|")[0].strip())
                    grid = []
                    for pf in line[2:].split("|")[1].strip().split(" "):

                        grid.append(float(pf))
                    self.pf_grids[res] = grid

            elif line.split(":")[0].strip() == "grid_size":
                self.grid_size = int(line.split(":")[1].strip())

            elif line.split(":")[0].strip() == "State":
                self.state_name = line.split(":")[1].strip()

            elif line.split(":")[0].strip() == "Molecule_Name":
                self.molecule_name = line.split(":")[1].strip()
        f.close()

    def clear_models(self):
        self.models = []

    def cluster_models_kmeans(self, nmodels, nclust):
        # Uses kmeans to cluster models
        from sklearn.cluster import KMeans
        mods = [m[1] for m in self.get_best_scoring_models(N=nmodels)]
        models = numpy.array(mods)
        kmeans = KMeans(n_clusters=nclust).fit(models)

        for i in range(nclust):
            unique, counts = numpy.unique(kmeans.labels_, return_counts=True)
            print(nclust, " | ", i, dict(zip(unique, counts))[i] *1.0/nmodels)
        '''
        for c in range(len(kmeans.cluster_centers_)):
            for i in range(c,len(kmeans.cluster_centers_)):
                print(i, c, numpy.linalg.norm(kmeans.cluster_centers_[c]-kmeans.cluster_centers_[i]))
        '''

    def get_datasets(self):
        if len(self.datasets) == 0:
            self.generate_datasets()
        return self.datasets

    def generate_datasets(self):
        self.datasets=[]
        for f in self.datafiles:
            self.datasets.append(hxio.import_json(self.path+f[0]))

    def get_models(self, return_pf=False):
        # returns all models stored in the POF (not the data file)
        if return_pf:
            return self.models_to_protection_factors(self.models)
        else:
            return self.models

    def get_all_models(self, return_pf=False):
        f = open(self.output_file, "r")
        models = []
        sigmas = []
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
                models.append((score, model_list)) 

                sigma_string = line[1:].split("||")[1].strip()
                sigma_list = []
                for s in sigma_string.split(" "):
                    sigma_list.append(s)
                sigmas.append(sigma_list)

        self.models = models

        return models, sigmas

    def get_pf_avg_sd(self, weighted=True):

        models = [m[1] for m in self.get_models()]     
        scores = numpy.array([m[0] for m in self.get_models()])

        if weighted:
            weights = numpy.exp(numpy.array(-1*(scores-min(scores))))
            mean_model = numpy.average(models, axis=0, weights=weights)
            stdev=[]
            
            for r in range(len(mean_model)):
                resi_vals = numpy.array([m[r] for m in models])
                stdev.append(numpy.average((resi_vals-mean_model[r])**2, weights=weights))
        else:
            mean_model = numpy.mean(models, axis=0)
            stdev = numpy.std(models, axis=0)

        return mean_model, stdev

    def calculate_shannon_bits(self, models=None, weighted=False):
        '''
        Calculate model probability distributions.
        Returns a len(residues) array of values
        '''
        if models==None:
            models = self.get_models()
        
        modvals = [m[1] for m in models]
        scores = [m[0] for m in models]

        bins = numpy.arange(0,self.grid_size)+0.5

        bits = []

        if weighted:
            weights = numpy.exp(numpy.array(-1*(scores-min(scores))))

        for r in range(len(modvals[0])):
            r_vals = [m[r] for m in modvals]
            hist = numpy.histogram(numpy.array(r_vals), bins=bins, density=True)[0]
            if sum(hist)==0:
                bits.append(0)

            nbins = len(hist)
            base_info = numpy.log(nbins)
            hist_info = 0

            for p in hist:
                if p != 0:
                    hist_info += p * numpy.log(1/p)

            bits.append(base_info - hist_info)
        return bits

    def calculate_accuracy(self, real_values, use_stdevs=False, weighted=True):
        '''
        real_values must be a list or array of len(residues). Any residue missing a value should be given a value of NaN.
        Returns cosine similarity for mean of the Pfs vs the real data
        '''

        meanmod, stdev = self.get_pf_avg_sd(weighted=weighted)

        if len(real_values) != len(meanmod):
            raise Exception("real_pf_value list must be equal to the number of residues,", len(meanmod))

        new_list_mod = []
        new_list_real = []
        for i in range(len(meanmod)):
            if real_values[i] != -1:
                new_list_mod.append(meanmod[i])
                new_list_real.append(real_values[i])

        return pairwise_distances([new_list_mod],[new_list_real], metric='cosine') * 10

    def get_scores(self, sorted=False):
        try:
            models = self.models
        except:
            models = self.get_all_models()

        scores = [mod[0] for mod in models]

        if sorted:
            return numpy.sort(scores)
        else:
            return scores


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
                    model_string = line[1:].split("|")[0].strip()
                    model_list = []
                    for m in model_string.split(" "):
                        model_list.append(int(m))
                    if sort_sectors:
                        model_list = self.sort_model_by_sector(model_list)
                    if return_pf:
                        ml1 = model_list
                        model_list = self.models_to_protection_factors(model_list)
                    best_scoring_models.append((score, model_list))
                    best_scoring_models = sorted(best_scoring_models, key=lambda x: x[0])

        self.best_scoring_models = best_scoring_models

        return best_scoring_models

    def models_to_protection_factors(self, models):
        # Input a list of list of integers.  
        # CHECK THAT THE MODEL SIZE IS CORRECT!

        if type(models[0]) is tuple:
            new_models = []
            for m in models:
                new_models.append(m[1])
            models = new_models
        elif type(models[0]) is not list:
            models = [models]

        protection_factor_models = []

        for m in models:
            # Model might be a tuple with a list
            pf_model = []
            for res in range(len(m)):
                rp1 = res+1
                if rp1 in self.observed_residues:
                    pf_model.append(float(self.pf_grids[res+1][m[res]-1]))
                else:
                    pf_model.append(numpy.nan)
            protection_factor_models.append(pf_model)

        return protection_factor_models

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

        out_model = numpy.zeros(len(model))
        # Loop over all sectors
        for s in self.sectors:
            if len(s) == 0:
                continue

            sec_model = model[s[0]-1:s[-1]]
            sort_model = numpy.sort(sec_model) #sort indexes in increasing order

            idx_zero = [index for index, v in enumerate(sec_model) if v == 0]
            idx_nonzero = [index for index, v in enumerate(sec_model) if v != 0]
            offset=0
            for i in range(len(s)):
                if i in idx_zero:
                    out_model[i] = 0
                    offset += 1
                else:
                    out_model[i+offset] = sort_model[i + len(idx_zero)-offset]
                out_model[idx_nonzero[i] - 1 + s[0]] = sort_model[i + len(idx_zero)]
        return out_model.astype(int)

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

        convergence_tuples.append(pct_tuples)

        del scores

    def calculate_peptide_scores(self, num_models = 1000000):
        '''
        Given a set of models, return the models a list of peptides and the score against each peptide
        '''
        # Get the peptides from the datasets
        peptides = self.get_all_peptides()

        pep_scores = []

        # Get the protection factors from the models and Pf grids, sigmas
        models = self.get_best_scoring_models(num_models, return_pfs=True)

        # Calculate fits 

        # Return a list of data peptides and their scores
        pass

    def calculate_peptide_chis(self, models="all"):

        chis = numpy.zeros(len(self.get_all_peptides()))
        models = self.get_models(return_pf=True)
        for m in models:
            chis += tools.get_peptides_chi_to_models(self.get_all_peptides(), m)
        return chis / len(models)

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
        return rand_smpl

class OutputAnalysis(object):
    '''
    Class that takes a list of output files for a single state and:
        * Concatenates all models into two ParseOutput objects

        * Runs Convergence
        * Runs Precision
        * Outputs a file with Pf probability distributions
    '''
    def __init__(self, output_files, analysis_output="analysis", shuffle=True, randint=89104):
        self.output_files = output_files
        self.get_output_file_consistency(output_files)
        self.output_directory = analysis_output # directory to place this analysis
        self.pof1, self.pof2 = self.split_into_two_POFs(shuffle, randint)
        self.sequence = self.pof1.get_sequence()

    def get_output_file_consistency(self, output_files):
        # For now, just ensure the sequence is the same and 
        pof_ref = ParseOutputFile(self.output_files[0])
        seq_ref = pof_ref.get_sequence()

        for of in output_files[1:]:
            seq = ParseOutputFile(of).get_sequence()
            if seq != seq_ref:
                raise Exception("All output files do not have the same sequence:", output_files)


    def calculate_rhat(self):
        # rhat = B/W
        # B = stdev of pooled sample of all MC iterations (between-chain variability)
        # W = within-chain variability

        # Returns a numpy array of rhat values for all residues

        # 1) get list of pofs:

        pofs = []
        for of in self.output_files:
            pofs.append(ParseOutputFile(of))

        n_resis = len(pofs[0].seq)
        all_models = []
        all_p_models = []
        self.all_models = all_models
        self.all_p_models = all_p_models

        for p in pofs:
            all_p_models.append([m[1] for m in p.get_all_models()[0]])
            all_models = all_models + all_p_models[-1]

        n = len(all_p_models[0])
        M = len(pofs)

        win_chain_var = numpy.zeros(n_resis)
        btw_chain_var = numpy.zeros(n_resis)
        # Get the average
        for i in range(n_resis):
            resi_means = []
            for p_models in all_p_models:
                resi_indices = [m[i] for m in p_models]
                win_chain_var[i] += numpy.var(resi_indices)
                resi_means.append(numpy.mean(resi_indices))

            btw_chain_var[i] = numpy.var(resi_means) * len(p_models)

        win_chain_var /= len(pofs)

        V_hat = (len(p_models)-1)/len(p_models) * win_chain_var + (len(pofs)+1)/(len(p_models)*len(pofs)) * btw_chain_var

        # return a list of psrf values for all residues
        psrf = numpy.zeros(n_resis)
        for i in range(n_resis):
            if win_chain_var[i] == 0:
                psrf[i] = 0.0
            else:
                psrf[i] = math.sqrt(V_hat[i]/win_chain_var[i])
            #resi_indices = numpy.array([m[i] for m in p_models])
            #var = resi_indices.var()
            #x = resi_indices - resi_indices.mean()
            #r = numpy.correlate(x, x, model='full')[-n:]

        # Just return the PRSF for now.
        # R_hat is a slight variation of this...I think?
        return psrf

    def calculate_neff(self):
        pass

    def split_into_two_POFs(self, shuffle, randint):
        '''
        Loop through the list of self.output_files and return two POF objects,
        each containing half of the models of interest
        '''

        import random
        # Split list of output files into two sets

        if len(self.output_files)==2:
            pof1 = ParseOutputFile(self.output_files[0])
            pof2 = ParseOutputFile(self.output_files[1])
            pof1.get_all_models()
            pof2.get_all_models()
            return pof1, pof2

        n_output_files = list(range(len(self.output_files)))

        random.seed(randint)    
        if shuffle:
            random.shuffle(n_output_files)

        if len(self.output_files)==3:
            pof1 = ParseOutputFile(self.output_files[n_output_files[0]])
            pof1.get_all_models()
            new_pof = ParseOutputFile(self.output_files[n_output_files[1]])
            concatenate_pofs(pof1, new_pof)
            pof2 = ParseOutputFile(self.output_files[n_output_files[2]])
            pof2.get_all_models()
            return pof1, pof2

        pof1 = ParseOutputFile(self.output_files[n_output_files[0]])
        pof1.get_all_models()

        for i in n_output_files[1:int(len(n_output_files)/2)]:
            new_pof = ParseOutputFile(self.output_files[n_output_files[i]])
            concatenate_pofs(pof1, new_pof)

        pof2 = ParseOutputFile(self.output_files[n_output_files[int(len(n_output_files)/2)]])
        pof2.get_all_models()
        for i in n_output_files[int(len(n_output_files)/2)+1:]:
            new_pof = ParseOutputFile(self.output_files[n_output_files[i]])
            concatenate_pofs(pof2, new_pof)

        return pof1, pof2

    def get_all_scores(self, sorted=True):
        # Plot the scores for all the models
        all_models = self.pof1.models + self.pof2.models
        scores = [mod[0] for mod in all_models]
        if sorted:
            return numpy.sort(scores)
        else:
            return scores

    def get_best_scoring_models(self, num):
        # Get the best scoring models from both independent sets
        new_pof = ParseOutpuFile(self.pof1.output_file)
        pof_all = concatenate_pofs(new_pof, self.pof2)
        return pof_all.get_best_scoring_models(num)

    def get_convergence(self, num_models="all"):
        return Convergence(self.pof1, self.pof2, num_models)

def concatenate_pofs(pof1, pof2):
    '''
    Function that concatenates the models and output_files from a list of pofs into pof1.

    Returns pof1
    '''
    if len(pof1.models)==0:
        all_models = pof1.get_all_models()
    else:
        all_models = pof1.models

    if len(pof2.models)==0:
        all_models += pof2.get_all_models()
    else:
        all_models += pof2.models

    pof1.models = all_models

    return pof1

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

    best_models = numpy.array(best_models)

    sorted_sectors = sector_sort(sectors)

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


def get_average_sector_values(exp_models, state):
    sector_model = get_sector_averaged_models(exp_models, state)
    avg = numpy.average(sector_model,0)
    std = numpy.std(sector_model,0)
    return (avg, std)

def get_cdf(exp_models):
    """ Takes a list of 1D exp_models and returns a sorted numpy array
    equivalent to the empirical density function for each residue. """
    exp_model_edf=numpy.empty((len(exp_models),len(exp_models[0])))
    A = numpy.array(exp_models)
    y = numpy.linspace(1./len(exp_models),1,len(exp_models))
    for i in range(len(exp_models[0])): 
        counts, edges = numpy.histogram(A[:,i], len(A), range=(-6,0), density=False) 
        exp_model_edf[:,i] = numpy.cumsum(counts*1.0/len(A))
    return exp_model_edf

def get_chisq(exp_models1, exp_models2, nbins):
    """ Takes two lists of exp_models and returns the chi2 value along the second axis """
    A = numpy.array(exp_models1)
    B = numpy.array(exp_models2)

    for i in range(269,len(exp_models1[0])):
        meanA = numpy.mean(A[:,i])
        ssd = numpy.std(A[:,i])**2 + numpy.std(B[:,i])**2 
        sstdev = numpy.sqrt( ssd / 5000 )
        meanB = numpy.mean(B[:,i])
        t = 1.96
        ci = t * sstdev
        dm = meanA - meanB

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

def calculate_autocorrelation(arr):
    n = len(arr)
    var = arr.var()
    arr = arr-arr.mean()
    r = np.correlate(arr, arr, mode = 'full')[-n:]
    assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(np.arange(n, 0, -1)))
    return result

def calculate_neff(models):
    # For a set of models, calculate the effective sample size per residue
    # Neff = n * det(S_cov)^(1/p) / det(MC_cov)^(1/p)
    #

    pass

