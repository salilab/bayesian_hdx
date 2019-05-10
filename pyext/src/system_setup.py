"""
   Classes to handle the HDX data hierarchy
"""
from __future__ import print_function
import hdx_models
import analysis
import input_data
import numpy
import scipy
import scipy.special
from numpy import linalg
import sys
from copy import deepcopy
import os.path

class HDXModel(object):
    '''
     Top-level HDX hierarchy representing a single macromolecular system

    HDXModel has (multiple) child objects HDXState, each representing different
    macromolecular states
    '''
    def __init__(self, name, inseq, offset, number_of_back_exchanged_amides = 2):
        '''
        Experimental variables contained in HDXModel:
        @param name - Unique identifier for this macromolecule
        @param inseq - FASTA sequence of macromolecule
        @param offset - offset to match FASTA sequence index to residue number in data
        @param number_of_back_exchanged_amides - The estimate for number of N-terminal amides that
                        will back exchange during analysis
        '''
        self.target_name=name
        self.states=[]
        self.score=0
        #offset is only used to provide offset between input sequence and fragment data
        #all calculations in this module will be based off of the sequence index of inseq, however
        #fragment data may be offset by some number
        self.offset=offset-1
        self.num_states=0
        self.seq=inseq
        self.num_res=len(self.seq)
        self.total_amides=self.calc_total_amides(self.seq)
        self.back_exchanged_amides = number_of_back_exchanged_amides

    def calc_total_amides(self, inseq):
        numP=inseq.count('P',2) + inseq.count('p',2)
        return len(self.seq)-numP

    def add_state(self, sname, mole_frac_liganded=1):
        self.num_states=self.num_states+1
        new_state=HDXState(self,sname,self.seq, self.offset, mole_frac_liganded)
        self.states.append(new_state)
        return new_state

    def calc_num_amides(self, inseq):
        numP=inseq.count('P',2) + inseq.count('p',2)
        return len(self.seq)-numP

    def get_apo_state(self):
        return self.states[0]



class HDXState(object):
    """
    Class representing a single protein state.
    Contains a list of Fragment objects unique to the state.
    """
    def __init__(self, model, state_name, seq, offset, mole_frac_liganded=0):
        '''
        @param model - The HDXModel instance
        @param state_name - A string for the state model
        @param seq - String containing the FASTA sequence
        @param offset - The offset between the PDB numbering and beginning of @seq string
        @param mole_frac_liganded - The mol fraction of the bound state in the liganded data
        '''
        self.state_name=state_name
        self.offset=offset
        self.seq=seq
        self.num_res=len(seq)
        self.frags=[]                    #List of fragments associated with this HDXState
        self.mole_fraction=mole_frac_liganded
        if mole_frac_liganded==0 and state_name!="apo" and state_name != "Apo":
            print("WARNING: Mole Fractions for state ",state_name," is zero.")

    def add_fragment(self, frag):
        '''
        Takes a fragment object and adds it to the HDXState fragment list
        Returns the fragment object
        '''
        if self.frag_seq_consistency(frag):
            self.frags.append(frag)
            return self.frags[-1]
        else:
            return None

    def create_fragment(self, seq, start, end):
        '''
        Manually creates a fragment object given fragment sequence, start and end
        residue and adds it to the end of HDXState fragment list

        Returns the fragment object
        '''
        frag_id=len(self.frags)
        new_frag=Fragment(seq, start, end, frag_id)
        if self.frag_seq_consistency(new_frag)==True:
            self.frags.append(new_frag)
            return new_frag
        else:
            print("This sequence", seq, "is inconsistent with model sequence in range", start, end)
            return None

    def remove_fragment(self, seq, start):
        '''
        Removes a fragment from the state's fragment list given sequence and start
        residue
        '''
        index=0
        for f in self.frags:
            if f.seq==seq and f.start_res==start:
                del self.frags[index]
                print(f.seq, "deleted")
            index+=1
        return f

    def frag_seq_consistency(self, frag):
        '''
        Returns True if fragment sequence and start residue aligns
        with macromolecule sequence

        Returns False with a warning if there is an inconsistency
        '''
        n_frag=0
        for n_seq in range(frag.start_res+self.offset,frag.end_res+self.offset):

            #print(frag.seq, frag.start_res-self.offset, n_seq)

            if frag.seq[n_frag]!=self.seq[n_seq]:
                print("Fragment ",frag.seq," does not match Sequence")
                print("Fragment position",n_frag+1," is ",frag.seq[n_frag])
                print("Sequence position",n_seq," is ",self.seq[n_seq])
                print("Offset is",self.offset)
                return False
            n_frag=n_frag+1
        return True

    def use_subset_of_fragments(self, start_frag, end_frag):
        '''
        Small hack function that pares the HDXState fragment list to
        a range of fragment indexes.  Useful for testing.
        '''
        new_frag_list=self.frags[start_frag:end_frag]
        self.frags=new_frag_list
        self.get_coverage(self.frags)
        self.get_sectors(self.frags)
        return self.frags

    def get_model(self):
        '''
        Returns model parent to this HDXState
        '''
        return self.model

    def get_coverage(self, frags=[]):
        '''
        Given a list of fragments, returns a vector of length
        len(seq) containing the per-residue coverage,
        defined as the number of times residue n
        is observed in the set of HDX fragments
        '''
        if len(frags)==0:
            frags = self.frags
        if len(frags)==0:
            print("No fragments imported into this state:", self.state_name)
        self.coverage=[0]*self.num_res
        for n in range(self.num_res):
            for f in frags:
                #cannot observe first two amides, do not count them in coverage
                if n >= f.start_res-self.offset and n < f.end_res-self.offset-1:
                    self.coverage[n]=self.coverage[n]+1
                    #print n, f.seq, f.start_res, f.end_res
            #print n,n-self.start_res,self.coverage[n-self.start_res]
        return self.coverage

    def get_all_sectors(self):
        return self.sectors

    def get_sectors(self, frags=[]):
        '''
        Given a list of fragments, returns a list of sector
        objects corresponding to the uniquely sampled
        segments in the sequence due to overlapping fragments
        '''
        if len(frags)==0:
            frags = self.frags
        if len(frags)==0:
            print("No fragments imported into this state:", self.state_name)
        self.sectors=[]
        sector_number=0
        start_res=0
        fraglist_m1=[]
        used_frags=set() # for debugging purposes
        for n in range(len(self.seq)+1):
            fraglist=[]
            # For this residue...what fragments cover it?
            for f in frags:
                if n >= f.start_res+self.offset + 2 and n <= f.end_res+self.offset:
                    fraglist.append(f)
                    used_frags.add(f)
            if len(fraglist)==0 and len(fraglist_m1) == 0:
                coverage=0
                #print("NO COVERAGE:", n, self.seq[n])#, fraglist
            elif len(fraglist)==0 and len(fraglist_m1) != 0:
                end_res=n-1
                #print("END_SECTOR:", n, self.seq[n], start_res, end_res, [(f.seq, f.start_res, f.end_res) for f in fraglist_m1])#, fraglist, fraglist_m1
                self.sectors.append(Sector(start_res, end_res, fraglist_m1, self.seq))
            elif fraglist==fraglist_m1:
                coverage=len(fraglist)
                #print("SAME_SECTOR:", n, self.seq[n])#, fraglist, fraglist_m1
            elif fraglist!=fraglist_m1 and fraglist_m1==[]:
                start_res=n
                sector_number+=1
                #print("START_SECTOR:", n, self.seq[n])#, fraglist, fraglist_m1
            elif fraglist!=fraglist_m1 and fraglist_m1!=[]:
                end_res=n-1
                #print("END/BEG_SECTOR", n, self.seq[n], start_res, end_res, [(f.seq, f.start_res, f.end_res) for f in fraglist_m1])#, fraglist, fraglist_m1
                self.sectors.append(Sector(start_res, end_res, fraglist_m1, self.seq))
                sector_number+=1
                start_res=n
            elif n==len(self.seq):
                self.sectors.append(Sector(start_res, n-1, fraglist_m1, self.seq))
            #print("XXX", n, self.seq[n], sector_number, len(used_frags), [(f.seq, f.start_res, f.end_res) for f in fraglist])#print f.seq, f.start_res, f.end_res, n, self.seq[n]
            fraglist_m1=fraglist

        #print("NSEC", len(self.sectors), len(self.frags), len(frags), len(used_frags))
        #print("FRAGS", len(frags), [f.seq for f in frags])
        return self.sectors

    def remove_experimental_outliers(self, frags, sigma, numsig=3):
        '''
        Given a list of fragments, and an experimental sigma,
        removes outlier %D values from the experimental data at the given
        numsig value.
        '''
        outliers=1
        while outliers > 0:
            outliers=0
            for f in frags:
                for t in f.timepoints:
                    outliers=outliers+t.remove_outliers(sigma, numsig)
                    if len(t.deut)==0:
                        f.timepoints.remove(t)
                        print(f.seq, "Timepoint: ", t.time, " removed")

    def get_frag(self, seq, start_res):
        ''' returns a fragment given sequence and start residue
        '''
        for f in self.frags:
            if f.seq==seq and f.start_res==start_res:
                frag = f
                break
            else:
                frag=""
        if frag == "":
            print("Frag", f.seq, "with start_res", start_res, "not found in state", self.state_name)
        return frag

    def get_all_frags(self):
        return self.frags

    def add_exp_model(self, exp_model):
        # Adds exp_model to state
        self.exp_model=exp_model


class Fragment(object):
    """ Class that stores a list of Timepoint data for each experimental HDX peptide fragment.
        Each Fragment object is bound to a single HDXState.
    """
    def __init__(self, inseq, start_res, end_res, frag_id = 0, sigma=1.0):
        #self.state=state
        self.id=frag_id
        self.seq=inseq
        self.timepoints=[]
        if end_res-start_res+1 != len(inseq):
            i=1
            #print "Length of fragment %s, %i, and residue numbers, %i, %i do not match" % (inseq, len(inseq), start_res, end_res)
        self.start_res=int(start_res)
        self.end_res=int(end_res)
        self.num_observable_amides=self.calc_num_observable_amides(inseq)
        self.sigma=sigma

    def calc_num_observable_amides(self, inseq):
        #number of observable amides is equal to fragment length - 2, minus remaining prolines
        num_amides=inseq.count('P',2) + inseq.count('p',2)
        return len(self.seq)-num_amides-2

    def get_num_observable_amides(self):
        return self.num_observable_amides

    def add_timepoint(self, time):
        tp=Timepoint(time)
        self.timepoints.append(tp)
        return tp

    def get_chi_value(self, sig):
        chi=0
        totd=0
        nrep=0
        for t in self.timepoints:
            t.get_model_avg()
            t.get_model_sd()
            for r in t.replicates:
                if t.model_avg is not None:
                    if sig=="model":
                        sig=t.sigma
                    chi+=(r.deut-t.model_avg)**2/sig
                    nrep+=1
                #print f.seq, chi, t.time, t.model_avg, t.deut
        if len(self.timepoints) > 0:
            self.chi = chi/nrep
        else:
            self.chi = 0
        #print self.seq, self.chi
        return self.chi

    def calculate_frag_score_freq_grid(self, freq_grid, exp_grid, sig, save=False, force=True):
        score=0
        for tp in self.timepoints:
            tp.calc_model_deut(freq_grid, exp_grid, self.num_observable_amides)
            score += tp.calculate_tp_score(grid, exp_grid, sig, self.num_observable_amides, force_calc=True)
        return score

    def calculate_frag_score(self, grid, exp_grid, sig, save=False, force=False):
        score = 0
        for tp in self.timepoints:
            tp.calc_model_deut(grid, exp_grid, self.num_observable_amides)
            score += tp.calculate_tp_score(grid, exp_grid, sig, self.num_observable_amides, force_calc=True)
        return score

    def set_tp_sigmas(self,sigma):
        for tp in self.timepoints:
            tp.sigma=sigma


class Sector(object):
    """ Each Sector object represents a portion of the peptide sequence with differential
        overlap by the MS peptide fragments.
        Sectors contain a list of Fragment objects which are contained in the Sector.
    """
    def __init__(self, start_res, end_res, fraglist, seq):
        self.fragments=fraglist
        self.id=[]
        self.start_res=start_res
        self.end_res=end_res
        self.coverage=len(fraglist)
        self.seq=seq[start_res:end_res+1]
        self.num_amides=self.calc_total_amides(self.seq)

    def calc_total_amides(self, seq):
        numP=seq.count('P') + seq.count('p')
        return len(self.seq)-numP


class Timepoint(object):
    '''
    Timepoint objects are bound to a specific Fragment, where Timepoint represents
    a regular observation time of the HDX experiment.
    The class contains both the experimental data (as a list of Replicate objects)
    as well as a list of calculated values from the forward model (self.models)
    '''
    def __init__(self, time, sigma0=5.0):
        '''
        @param time - Time in seconds
        @param sigma0 - Initial estimate of timepoint error sigma in pctD units.
        '''
        self.time=time
        self.models=[]
        self.num_replicates=0
        self.replicates=[]
        self.sigma=sigma0

    def add_replicate(self, deut=None, sat=1.0, recovery=1.0, temp=293.15):
        self.num_replicates+=1
        self.replicates.append(Replicate(deut, sat, recovery, temp))

    def calculate_tp_score(self, freq_grid, exp_grid, sig, num_amides, force_calc=False):
        """ given the exp_grid of exponential rates and number of observations
        in each grid, calculates the Bayesian score """
        score = 0
        for r in self.replicates:
            if force_calc:
                self.calc_model_deut(freq_grid, exp_grid, num_amides)
                rep_score = r.calculate_replicate_score(self.model_deut, sig)
            else:
                try:
                    rep_score = r.calculate_replicate_score(self.model_deut, sig)
                except:
                    self.calc_model_deut(freq_grid, exp_grid, num_amides)
                    rep_score = r.calculate_replicate_score(self.model_deut, sig)

            # This if statement prevents log overflows when score ~0
            if r.calculate_replicate_score(self.model_deut, sig) == 0.0:
                score += 10000
            else:
                score += -1.0*numpy.log(rep_score)
        return score / len(self.replicates)


    def calc_model_deut(self, freq_grid, exp_grid, num_amides):
        # Given an exp_frequency grid and the exp_grid, calculate the deuteration of the model
        # at this timepoint.
        deut=0
        for n in range(len(exp_grid)):
            #print(n, freq_grid)
            exchange_rate = 10**exp_grid[n]
            num_amides_at_this_rate = freq_grid[n]
            deut += num_amides_at_this_rate*(1-numpy.exp(-exchange_rate*self.time))
        self.model_deut = deut / num_amides * 100
        #print(self.time, grid, exp_grid, self.model_deut)
        return self.model_deut

    def get_avg_sd(self):
        if len(self.replicates) < 1:
            #print "No replicates in timepoint, ", self.time
            self.avg=None
            self.sd=None
        elif len(self.replicates) <= 2:
            #print "Not enough replicates, ", self.time
            self.avg=sum(r.deut for r in self.replicates)/len(self.replicates)
            self.sd=None
        else:
            sum_deut=0
            sumsq_deut=0
            for r in self.replicates:
                sum_deut=sum_deut+float(r.deut)
                sumsq_deut=sumsq_deut+r.deut**2
            self.avg=sum_deut/self.num_replicates
            self.sd=numpy.sqrt((self.num_replicates*sumsq_deut-sum_deut**2)/self.num_replicates**2)
        return self.avg, self.sd

    def get_model_avg(self):
        if len(self.models) < 1:
            print("No models imported into timepoint ", self.time)
            self.model_avg=None
            return None
        else:
            #print "Not enough replicates, ", self.time
            self.model_avg=numpy.average(self.models)
        return self.model_avg

    def get_model_sd(self):
        if len(self.models) <= 2:
            print("Not enough models to calculate SD ", self.time)
            self.model_sd=None
        else:
            #print "Not enough replicates, ", self.time
            self.model_sd=numpy.std(self.models)
        return self.model_sd

    def add_model_deuteration_value(self, deut):
        self.models.append(deut)

    def clear_model_deuteration_values(self):
        self.models=[]

    def remove_outliers(self, sigma, numsig=3):
        mean,sd=self.get_avg_sd()
        outliers=0
        for d in self.deut:
            if abs(d-mean)/sigma>numsig:
                self.deut.remove(d)
                print("outlier removed: ", d, mean, sigma)
                outliers=outliers+1
                self.num_replicates=self.num_replicates-1
        return outliers

    def set_sigma(self, sigma):
        '''
        @param sigma - float to set timepoint sigma; or "exp_sd" to use SD from
             experimental data.
        '''
        if sigma=="exp_sd":
            if self.num_replicates > 2:
                avg, self.sigma=self.get_avg_sd()
            else:
                raise Exception("Not enough models to calculate SD")
        else:
            self.sigma=sigma


class Replicate(object):
    """ Replicate objects store the individual data observations and have the ability to
        contain experiment-specific data, in the case of merging multiple data sets under
        different conditions
    """
    def __init__(self,deut, sat=1.0, recovery=1.0, temp=293):
        self.deut=deut
        self.sat=sat
        self.nearest_gridpoint=-1
        self.recovery=recovery
        self.temp=temp

    def calculate_replicate_score(self, model_deut, sig):
        return self.gaussian_model(self.deut, model_deut, sig)

    def lognormal_model(self,exp,model,sig):
        return numpy.exp(-( (numpy.log(model)-numpy.log(exp) )**2)/(2*sig**2))/(sig*numpy.sqrt(numpy.pi)*model)

    def gaussian_model(self,exp,model,sig):
        return numpy.exp(-((model-exp)**2)/(2*sig**2))/(sig*numpy.sqrt(2*numpy.pi))

    def calc_nearest_gridpoint(self, grid):
        self.nearest_gridpoint = (numpy.abs(grid-self.deut)).argmin()

    def get_nearest_gridpoint(self, grid):
        if self.nearest_gridpoint==-1:
            self.calc_nearest_gridpoint(grid)

        return self.nearest_gridpoint
