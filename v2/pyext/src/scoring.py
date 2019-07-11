"""
   Scoring functions for HDX models
"""
from __future__ import print_function
import numpy
import tools
import scipy
import scipy.stats
import math


class ScoringFunction(object):
    '''
    A scoring function consists of a forward model + noise model along with
    priors for each parameter within the model and the data.

    It depends on a model representation (from bayesian_hdx::model -- should be remaned Representation)
    a Forward model, a Noise model and sets of priors for each type of parameter

    '''
    def __init__(self, fwd_model):
        '''
        Input is a model, which must have a function get_model(), which returns the current model
        self.representation 
        '''
        self.forward_model = fwd_model
        self.priors = []

    def add_prior(self, prior):
        self.priors.append(prior)


    def evaluate(self, model_values=None, peptides=None):
        '''
        Evaluates an entire model based on a single set of model values.
        '''
        total_score = 0

        if model_values is None:
            model_values = self.representation.get_current_model()

        total_score = self.forward_model.evaluate(model_values, peptides)
        #print(" FM:",total_score, model_values)
        for p in self.priors:
            #print("  ",p.evaluate(model_values))
            total_score += p.evaluate(model_values)

        #print("   ", model_values[2:-1])
        return total_score



class ProtectionFactorNaturalAbundancePrior(object):
    '''
    A prior on the likelihood of observing a given protection factor given no other information

    Can choose among built-in functions or supply your own function.

    Functions are evaluated by numpy.interp()
    '''
    def __init__(self, prior_type, input_bins=None, input_pfs=None, gaussian_parameters=[(1,0.8,3), (5,1.6,80)]):
        '''
        prior_type can be "bmrb", "knowledge", "user" or "uninformative"
        '''
        self.bins = []
        self.probs = []
        self.prior_type = prior_type
        if prior_type == "bmrb":
            for i in sorted(bmrb_pf_histograms.keys()):
                self.bins.append(i)
                self.probs.append(bmrb_pf_histograms[i])
        elif prior_type == "knowledge":
            for i in sorted(knowledge_pf_histograms.keys()):
                self.bins.append(i)
                self.probs.append(knowledge_pf_histograms[i])   
        elif prior_type == "uninformative":
            for i in range(-2,14):
                self.bins.append(i)
                self.probs.append(1.0) 
        elif prior_type == "Gaussian":
            self.gaussians = []
            for g in gaussian_parameters:
                self.gaussians.append((scipy.stats.norm(g[0], g[1]),g[2]))

        elif prior_type == "user":
            if input_bins is not None and input_pfs is not None:
                if len(input_bins) == len(input_pfs):
                    self.bins = input_bins
                    self.probs = input_pfs
                else:
                    raise Exception("scoring.SingleProtectionFacotrPrior.set_prior: Length of bins and probabilities is not the same")
            else:
                raise Exception("scoring.SingleProtectionFacotrPrior.set_prior: Must supply input_bins and input_pfs for a user-supplied prior")
        # make sure we normalize!
        self.probs = numpy.linalg.norm(self.probs)

    def evaluate_pf_likelihood(self, pf):
        return numpy.interp(pf, self.bins, self.probs)

    def evaluate(self, model):
        score = 0
        if self.prior_type == "Gaussian":
            for m in model:
                prob = 1
                for g in self.gaussians:
                    prob *= g[0].pdf(m)*g[1]
                    #print(g[0].pdf(m), g[1], m, prob)
                score += -1*numpy.log(prob +0.0000000001)
        else:
            for m in model:
                score += -1*numpy.log(self.evaluate_pf_likelihood(m))
        return score


class SigmaPriorCauchy(object):
    '''
    A prior that is applied to the sigma for each timepoint or data
    '''
    def __init__(self, state, sigma_estimate=1.0, prior_scale=20.0):
        self.state = state
        self.scale = prior_scale

    def evaluate(self):
        # Prior on the sigma value. Long tailed to allow for outlier values.
        sigma_prior_score = 0
        for d in self.state.data:
            for p in d.get_peptides():
                for tp in p.get_timepoints():
                    sigma_prior_score += (1 / tp.sigma**2) * math.exp(-tp.sigma0**2 / tp.sigma**2)

        return sigma_prior_score * self.prior_scale


class LogSigmaPriorNormal(object):
    '''
    A prior that is applied to the log of the sigma.
    Keeps the value of sigma positive and 
    '''
    def __init__(self, state, log_sigma_estimate=0.0, sd=1.0, prior_scale=1.0):

        self.state = state
        self.prior_scale = prior_scale
        self.point_estimate = log_sigma_estimate
        self.sd = sd
        self.distribution = scipy.stats.norm(log_sigma_estimate, sd)

    def evaluate(self, model):
        # Prior on the sigma value. Long tailed to allow for outlier values.
        sigma_prior_score = 0
        for d in self.state.data:
            for p in d.get_peptides():
                for tp in p.get_timepoints():
                    sigma_prior_score += -1*numpy.log(self.distribution.pdf(numpy.log(tp.sigma)))

        return sigma_prior_score * self.prior_scale


class FlatDistribution(object):
    '''
    Simple object to return a flat prior.

    Compatible with scipy.stats distributions that contain a pdf() function

    Returns the same value regardless of input.
    '''
    def __init__(self, scale=0.1):
        self.scale = scale

    def pdf(self, value):
        return self.scale


class ResiduePfPrior(object):
    '''
    Given a set of estimates for the protection factor of each residue, apply a 
    Gaussian prior on that residue.

    Ideally determined from an MD simulation or, perhaps, NMR data of certain residues
    '''
    def __init__(self, pf_estimates, scale=1.0):
        '''
        pf_estimates is in the form (mean, sd).
        For residues without an estimate, denoted by an SD < 0, the prior is flat for any value
        '''
        self.prior_scale = scale
        self.priors = []
        for p in range(len(pf_estimates)):
            prior = pf_estimates[p]
            if prior[1] < 0:
                self.priors.append(FlatDistribution())
            else:
                self.priors.append(scipy.stats.norm(prior[0], prior[1]))

    def evaluate(self, model):
        if len(model) != len(self.priors):
            raise Exception("ERROR scoring.ResiduePfPrior: The length of the Pf prior is not th same as the model")

        score = 0
        for p in range(len(model)):
            #print(p, model[p], -1*numpy.log(self.priors[p].pdf(model[p])))
            score += -1*numpy.log(self.priors[p].pdf(model[p]))

        return score * self.prior_scale



class GaussianNoiseModel(object):
    '''
    This is a Forward Model.  Actually a Forward Model plus noise model.
    It converts a model (set of protection factors) into an expected value for each piece of data 
    and then derives a likelihood for the data:model pair using the noise model.

    This model gathers its standard deviation parameters from the individual timepoint objects.
    '''
    def __init__(self, state, truncated=False, bounds=(None, None)):
        self.truncated=truncated

        self.state = state

        if truncated:
            # 10 and -10 are numerically equivalent to no truncation factor when
            # data is in the range of 0-->1
            if bounds[1] is None:
                self.upper_bound = 10
            else:
                self.upper_bound = bounds[1]

            if bounds[0] is None:
                self.lower_bound = -10
            else:
                self.lower_bound = bounds[0]

    def replicate_score(self, model, exp, sigma):
        # Forward model
        #priors = self.model_prior(protection_factor) * self.exp_prior(exp) * self.sigma_prior()
        raw_likelihood = math.exp(-((model-exp)**2)/(2*sigma**2))/(sigma*math.sqrt(2*numpy.pi))
        if self.truncated:
            raw_likelihood *= 1/ ( 0.5 * ( scipy.special.erf( (self.upper_bound-exp)/sigma * math.sqrt(3.1415) ) - scipy.special.erf( (self.lower_bound-exp)/sigma * math.sqrt(3.1415) ) )  )
        return raw_likelihood

    def peptide_confidence_score(self, peptide):
        # User-definable function for converting peptide confidence into a likelihood.
        # The function must be evaluatable between 0 and 1. 
        pass 


    def calculate_peptides_score(self, peptides, protection_factors):
        '''
        Will deprecate calculate_dataset_score. Given a list of peptides,
        calculate the score. Useful for calculating changes that only affect
        a susbset of peptides.
        '''
        self.state.calculate_residue_incorporation(protection_factors)
        
        total_score = 0
        if peptides is None:
            # peptides are the ones that we recalculate
            peptides = self.state.get_all_peptides()
            non_peptides = []
        else:
            # non_peptides are the ones where we simply read the score from last time (didn't change)
            non_peptides = list(set(self.state.get_all_peptides())-set(peptides))
        #print("PEP:", len(peptides), len(non_peptides))
        for pep in peptides:
            peptide_score = 0

            d = pep.get_dataset()

            observable_residues = pep.get_observable_residue_numbers()

            # Cycle over all timepoints
            for tp in pep.get_timepoints():
                # initialize tp score to the sigma prior
                tp_score = 0
                #model_tp_raw_deut = self.sum_incorporations(self.calculate_residue_incorporation(self.output_model.model_protection_factors)[d], observable_residues, tp.time)
                model_tp_raw_deut = 0

                #try:
                #    self.state.residue_incorporations[d][0][tp.time]
                #except:


                for r in observable_residues:
                    model_tp_raw_deut += self.state.residue_incorporations[d][r][tp.time]

                #------------------------------
                # Here is where we would add back exchange estimate
                #------------------------------

                # Convert raw deuterons into a percent
                model_tp_deut = float(model_tp_raw_deut)/pep.num_observable_amides * 100

                # Calculate a score for each replicate
                for rep in tp.get_replicates():
                    replicate_likelihood = self.replicate_score(model=model_tp_deut, exp=rep.deut, sigma=tp.sigma) 
                    if replicate_likelihood <= 0:
                        rep.set_score(10000000000)
                    else:
                        rep.set_score(-1*math.log(replicate_likelihood))

                    tp_score += rep.get_score()

                # Set the timepoint score
                tp.set_score(tp_score)

                peptide_score += tp_score

            total_score += peptide_score
            #print(pep.sequence, peptide_score, protection_factors)

        for pep in non_peptides:
            #print(pep.sequence, pep.get_score(), protection_factors)
            total_score += pep.get_score()

        self.total_score = total_score
        return total_score

    def evaluate(self, model, peptides):
        return self.calculate_peptides_score(peptides, model)



