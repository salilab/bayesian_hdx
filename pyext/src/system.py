"""
   Classes to handle the HDX data hierarchy
"""
from __future__ import print_function
#import hdx_models
#import analysis
from scoring import GaussianNoiseModel
import hxio
import numpy
import scipy
import tools
#from model import ResidueGridModel
#import scipy.special
from numpy import linalg
import sys
from copy import deepcopy
import os.path



class System(object):
    """ Class defining a set of macromolecular objects
    within an HX Experiment
    """
    def __init__(self, output_dir=None, noclobber=True):
        self.macromolecules = []
        if output_dir is not None:
            self.output = hxio.Output(self, output_directory=output_dir, noclobber=noclobber)
        else:
            self.output = None

    def add_macromolecule(self, sequence, name=None):
        """ add a Macromolecule to the experiment. Sequence can either be
        a string or a filename.
        @param sequence - FASTA string OR filepath
        @param name - name for FASTA string OR list of fasta IDs in fastafile 
            (first word after '>' up until separator ' ', ':', ';')
        """

        # Check that sequence is a string
        if not isinstance(sequence, basestring):
            raise Exception("Please input a string or filepath as the sequence input")

        # if there is a period, then we will assume that it is a filepath
        if "." in sequence:
            try: 
                open(sequence, "r")
                seqs = io.read_fasta(sequence)
            except:
                raise Exception("Cannot open filename " + sequence)

            for s in seqs:
                if s[0] in name or s[0] == name:
                    self.macromolecules.append(s[0], s[1])

        else:
            if name is None:
                raise Exception("Please give a name to molecule with sequence" + sequence)
            self.macromolecules.append(Macromolecule(self, name, sequence))

        return self.macromolecules[-1]

    def get_molecules(self):
        return self.macromolecules 

    def get_output(self):
        return self.output    



class Macromolecule(object):
    def __init__(self, system, name, sequence, initialize_apo=True):
        self.name = name
        self.system = system
        self.sequence = sequence
        self.states = []

        if initialize_apo:
            self.add_state("Apo", None)

    def get_sequence(self):
        return self.sequence

    def add_state(self, name, perturbations):
        new_state = State(self, name, perturbations)
        self.states.append(new_state)
        return new_state

    def get_state(self, state_number=0, name=None):
        if name is not None:
            if name in [s.get_name() for s in self.states]:
                for s in self.states:
                    if s.get_name() == name:
                        return s
            else:
                raise Exception(name + " is not in the list of states:" + str([s.get_name() for s in self.states]))

        else:
            if state_number > len(self.states):
                raise Exception("There are only " + str(len(self.states)) + " states for molecule " + self.name)

            return self.states[0]

    def get_states(self):
        return self.states

    def get_apo_state(self):
        return self.states[0]




class HDXModel(object):
    '''
     Top-level HDX hierarchy representing a single macromolecular system
        
    HDXModel has (multiple) child objects HDXState, each representing different 
    macromolecular states
    '''
    def __init__(self, name, inseq, offset, number_of_back_exchanged_amides=2):
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
        self.sectors = []

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


class State(object):
    """ add one or more perturbations to the base system.
    A perturbation could be a small molecule
    a point mutation or a new complex
    """
    def __init__(self, mol, name, perturbations=None, 
                    output_model=None,
                    scoring_function=GaussianNoiseModel()):
        self.name = name
        self.macromolecule = mol
        self.perturbations = []
        self.sequence = self.macromolecule.get_sequence()

        if perturbations is not None:
            self.add_perturbations(perturbations)

        self.data = []
        self.sectors = []
        self.output_model = output_model
        self.scoring_function = scoring_function
        self.residue_sector_dictionary = {}

    def get_output(self):
        return self.macromolecule.system.get_output()

    def get_name(self):
        return self.name

    def set_output_model(self, model):
        '''adds a model object
        '''
        self.output_model = model
        self.has_model=True

    def set_scoring_function(self, scoring_function=GaussianNoiseModel):
        self.scoring_function = scoring_function
        self.has_data = True

    def get_output_model(self):
        return self.output_model

    def has_data(self):
        if len(self.data) == 0:
            return False
        else:
            return True

    def get_datasets(self):
        return self.data

    def add_perturbation(perturbation, clear=False):
        if type(perturbations) is tuple:
            perturbations = [perturbations]
        self.perturbations = perturbations
        for pert in self.perturbations:
            if pert[0]=="mutation":
                AAi = pert[1][0]
                AAm = pert[1][-1]
                resnum = pert[1][1:-1]
                if self.seq[int(resnum)-1] == AAi:
                    self.seq[int(resnum)-1] = AAm
                else:
                    return Exception("Mutation" + pert[1] + "cannot be applied to residue" + self.seq[resnum-1] + resnum)

    def get_sequence(self):
        return self.sequence

    def add_dataset(self, dataset):
        for p in dataset.get_peptides():
            if not self.peptide_sequence_consistency(p):
                raise Exception("Exiting at State.add_dataset. Peptide " + p.sequence + " does not match the sequence")
        self.data.append(dataset)
        dataset.set_state(self)

    def peptide_sequence_consistency(self, peptide):
        '''
        Returns True if peptide sequence and start residue aligns 
        with macromolecule sequence

        Returns False with a warning if there is an inconsistency
        '''
        peptide_position = 0

        for residue_number in range(peptide.start_residue, peptide.start_residue + len(peptide.sequence)):
            peptide_position = residue_number - peptide.start_residue
            if peptide.sequence[peptide_position] != self.sequence[residue_number-1]:
                print("Peptide ", peptide.sequence, " does not match Sequence")
                print("Peptide position ", peptide_position+1, " is ", peptide.sequence[peptide_position])
                print("Sequence position ", residue_number, " is ", self.sequence[residue_number])
                return False
        return True

    def get_macromolecule(self):
        return self.macromolecule

    def get_coverage(self, peptides=None):
        '''
        Given a list of peptides, returns a vector of length 
        len(seq) containing the per-residue coverage,
        defined as the number of times residue n
        is observed in the set of HDX fragments
        '''
        if peptides is None:
            peptides = self.peptides

        if len(peptides)==0:
            print("No peptides imported into this state:", self.state_name)

        #initialize to zero coverage for 
        self.coverage = numpy.zeros(len(self.sequence))

        for n in range(len(self.sequence)):
            for p in peptides:
                #cannot observe first two amides, do not count them in coverage
                if n >= f.start_res and n < f.end_res-1:
                    self.coverage[n]+=1
        return self.coverage

    def get_all_peptides(self):
        peptides = []
        for d in self.data:
            peptides += d.get_peptides()
        return peptides

    def get_sector_consolidated_model(self, model):
        '''Given a residue model, return the sector-consolidated
        model, represented as a tuple of two lists.
        ([residues], [values])
        '''
        if self.sectors == []:
            self.calculate_sectors()
        sector_consolidated_model = []
        for s in self.sectors:
            resis = s.get_residues()
            sector_values = []
            for r in resis:
                if int(model[r-1])==0:
                    resis.remove(r)
                else:
                    sector_values.append(int(model[r-1]))
            sector_consolidated_model.append((sorted(resis, key=lambda x: x[0]), sector_values))

        return sector_consolidated_model

    def calculate_sectors(self, peptides=None):
        '''
        Given a list of peptides or datasets, returns a list of sector
        objects corresponding to the uniquely sampled
        segments in the sequence due to overlapping fragments
        '''

        if peptides is None:
            peptides = self.get_all_peptides()
        #print(peptides)
        if len(peptides)==0:
            print("No peptides imported into this state:", self.name)

        # Overwrite the sector list each time
        self.sectors=[]

        residue_peptide_sets = dict({})
        self.residue_peptide_sets = residue_peptide_sets

        # First, determine all observable peptides that cover each residue
        # Store unique peptide ids in a dictionary of sets: residue_peptide_ids
        # with the dictionary key = residue number
        for i in range(1, 1+len(self.sequence)):
            peptide_ids = set()
            for p in peptides:
                if i in p.get_observable_residue_numbers():
                    peptide_ids.add(p)
            residue_peptide_sets[i] = peptide_ids

        # Second, collect all residues with the same
        # peptide_id sets in the sector_dictionary
        sector_number = 0
        self.sector_dictionary = dict()

        while len(residue_peptide_sets) > 0:
            #print(x, residue_peptide_sets)
            # Get the first element in the dictionary
            first_key = list(residue_peptide_sets)[0]
            sector_peptide_ids = residue_peptide_sets[first_key]
            residue_peptide_sets.pop(first_key)
            if len(sector_peptide_ids) > 0:
                sector_residues = set()
                sector_residues.add(first_key)

                # Compare peptide ids for each residue.
                # Same peptide ids = same sector
                for j in list(residue_peptide_sets):
                    peptide_ids = residue_peptide_sets[j]
                    if peptide_ids == sector_peptide_ids:
                        # Add residue to the sector list and remove it from the dictionary
                        sector_residues.add(j)
                        residue_peptide_sets.pop(j)

                self.sector_dictionary[sector_number] = sector_residues
                new_sector = Sector(self, sector_residues, sector_peptide_ids, sector_number)
                for s in sector_residues:
                    self.residue_sector_dictionary[s] = new_sector
                self.sectors.append(new_sector)
                sector_number+=1

        # Determine the list of observed_residues
        resis = set()
        for s in self.sectors:
            for r in s.residues:
                resis.add(r)
        self.observed_residues = list(resis)

        return self.sectors

    def get_sectors(self):
        return self.sectors

    def create_simulated_data(self):
        data = data.Dataset()

    def calculate_score(self, model):
        # @param model is a proposed model for the system
        '''
        if not self.has_model:
            raise Exception("System.calculate_score: Cannot calculate score without a model!!")
        
        if not self.has_data():
            raise Exception("System.calculate_score: Cannot calculate score without data!!")

        if not self.has_scoring_function:
            raise Exception("System.calculate_score: Cannot calculate score without a scoring function!!")
        '''
        # First, use the output_model to convert the proposed model into residue resolved protection factors
        protection_factors = self.output_model.convert_model_to_protection_factors(model)
        # and calculate the prior on these proposed observations

        total_score = 0
        # Second, pass these protection factors to each dataset object
        for d in self.data:
            dataset_score = self.calculate_peptides_score(d.get_peptides(), protection_factors)
            #print("----", d, dataset_score)
            total_score += dataset_score

        return total_score

    def initialize(self, init_model="random"):
        # Ensure that all components of the state are initialized
        # State
        # Datasets
        # output_model
        # Sampler
        self.calculate_sectors()
        for d in self.data:
            d.calculate_observable_rate_bounds()
            # Calculate the prior on the deuterium incorporation for each replicate
            #d.calculate_replicate_priors(self.scoring_function)

        if init_model=="random":
            self.output_model.generate_model(initialize=True)
        else:
            self.output_model.generate_model(random=False, value=init_model, initialize=True)

        self.calculate_score(self.output_model.model)

    def calculate_peptides_score(self, peptides, protection_factors):
        '''
        Will deprecate calculate_dataset_score. Given a list of peptides,
        calculate the score.  Useful for calculating changes that only affect
        a susbset of peptides.
        '''
        tools.get_residue_peptide_deuteration_at_each_timepoint(peptides, protection_factors)

        total_score = self.scoring_function.protection_factor_prior(protection_factors)

        scoring_function = self.scoring_function

        for pep in peptides:
            sigma0 = pep.dataset.sigma_estimate
            peptide_score = 0
            for tp in pep.get_timepoints():
                # initialize tp score to the sigma prior
                tp_score = -1*numpy.log(scoring_function.experimental_sigma_prior(tp.sigma, sigma0))

                # Get deuteration percent of this timepoint with the given Pfs
                model_tp_deut = tp.get_model_deuteration()
                # Convert raw deuterons into a percent
                model_tp_deut=float(model_tp_deut)/pep.num_observable_amides * 100

                # Calculate a score for each replicate
                for rep in tp.get_replicates():
                    #####
                    replicate_likelihood = scoring_function.replicate_score(model=model_tp_deut, exp=rep.deut, sigma=tp.sigma) 
                    rep.set_score(-1*numpy.log(replicate_likelihood))

                    tp_score += rep.get_score()

                tp.set_score(tp_score)
                #print("TP_SCORE", pep.sequence, tp.time, "||", model_tp_deut, tp.get_replicates()[0].deut, "||", tp.get_score(), tp.sigma)#, -1*numpy.log(scoring_function.experimental_sigma_prior(tp.sigma, sigma0)))
                peptide_score += tp_score
                #print(tp.time, tp_score, len(pep.get_timepoints()))

            total_score += peptide_score
            #print(pep.sequence, peptide_score, tp_score, tp.time)

        self.total_score = total_score
        return total_score

    def get_observed_residues(self):
        '''
        Returns a list of residues that are observed in any dataset
        '''
        resis = set()
        for s in self.sectors:
            for r in s.get_residues():
                resis.add(r)
        self.observed_residues = list(resis)
        return self.observed_residues

    def consolidate_model_to_sectors(self, model):
        # Given a model, return a list of lists with 
        pass


class Sector(object):
    """ Each Sector object represents a portion of the peptide sequence with differential
        overlap by the MS peptide fragments.
        @param state - The macromolecule state to which this sector belongs
        @param residues - list or set of residue numbers in the sector
        @param peptide_ids - list or set of peptide ids used in this sector
        @param sector_number - the unique numerical id for this sector

    """
    def __init__(self, state, residues, peptides, sector_number):
        self.state = state
        self.residues = residues
        self.peptides = peptides
        self.id = sector_number
        self.num_amides = len(self.residues) # self.get_number_of_amides()
        self.length = len(self.residues)

    def get_number_of_residues(self):
        return len(self.residues)

    def get_number_of_amides(self):
        amides = 0
        for i in self.residues:
            if self.state.sequence[i-1]!="P":
                amides += 1
        return amides

    def get_coverage(self):
        return len(peptide_ids)  

    def get_length(self):
        return self.length

    def calculate_sector_score(self, protection_factors):
        self.score = self.state.calculate_peptides_score(self.peptides, protection_factors)
        return self.score

    def set_score(self, score):
        self.score = score

    def get_score(self):
        return self.score

    def get_residues(self):
        return self.residues

    def get_peptides(self):
        return self.peptides


class Residue(object):
    # Each residue can hold the score for itself
    def __init__(self, state, number, residue_type):
        self.state = state
        self.number = number
        self.residue_type = residue_type

    def set_intrinsic_rate(self, rate):
        self.intrinsic = rate


def setup_single_state(sequence, name):
    '''Simple function that sets up a System with one Macromolecule
    and one State.  Returns the State object.
    '''
    sys = System()
    mol = sys.add_macromolecule(sequence, name)
    return mol.get_apo_state()




