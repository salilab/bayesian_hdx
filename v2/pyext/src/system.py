"""
   Classes to handle the HDX data hierarchy
"""
from __future__ import print_function
import io
from scoring import ScoringFunction, GaussianNoiseModel
import hxio
import numpy
import math
import tools
#from model import ResidueGridModel
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


    def add_macromolecule(self, sequence, name=None, initialize_apo=True):
        """ add a Macromolecule to the experiment. Sequence can either be
        a string or a filename.
        @param sequence - FASTA string OR filepath
        @param name - name for FASTA string OR list of fasta IDs in fastafile 
            (first word after '>' up until separator ' ', ':', ';')
        """

        # Check that sequence is a string
        if not isinstance(sequence, str):
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
            self.macromolecules.append(Macromolecule(self, name, sequence, initialize_apo))

        return self.macromolecules[-1]

    def get_macromolecules(self):
        return self.macromolecules 

    def get_output(self):
        return self.output  

    def initialize_output(self):
        for m in self.macromolecules:
            for s in m.get_states():
                self.output.initialize_output_model_file(s, s.output_model.pf_grids)

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

    def add_state(self, name, perturbations=None):
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

class State(object):
    """ add one or more perturbations to the base system.
    A perturbation could be a small molecule
    a point mutation or a new complex
    """
    def __init__(self, mol, name, perturbations=None, 
                    output_model=None,
                    scoring_function=None):
        self.name = name
        self.macromolecule = mol
        self.perturbations = []
        self.has_model = False
        self.has_scoring_function = False
        self.has_model = False
        self.observable_residue_numbers = []
        self.sequence = self.macromolecule.get_sequence()

        if perturbations is not None:
            self.add_perturbations(perturbations)

        self.data = []
        self.sectors = []
        if output_model is not None:
            self.set_output_model(output_model)
        else:
            self.output_model = None

        if scoring_function is None:
            self.scoring_function = ScoringFunction(GaussianNoiseModel(self))
        self.residue_sector_dictionary = {}
        self.score=10000000

    def get_output(self):
        return self.macromolecule.system.get_output()

    def get_scoring_function(self):
        return self.scoring_function

    def get_name(self):
        return self.name

    def set_output_model(self, model):
        '''adds a model object
        '''
        self.output_model = model
        self.has_model = True
        return self.output_model

    def get_exchanging_residues(self):
        '''
        Returns a list of residue numbers that exchange 
        (simply the list of non-proline residues)
        '''
        exchanging_residues = []
        for i in range(len(self.sequence)):
            if self.sequence[i] != "P":
                exchanging_residues.append(i+1)

        return exchanging_residues

    '''
    def set_scoring_function(self, scoring_function=GaussianNoiseModel):
        self.scoring_function = scoring_function
        self.has_scoring_function = True
    '''

    def get_output_model(self):
        return self.output_model

    def has_data(self):
        if len(self.data) == 0:
            return False
        else:
            return True

    def get_datasets(self):
        return self.data

    def add_perturbation(self, perturbation, clear=False):
        if type(perturbation) is tuple:
            perturbations = [perturbation]
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
            peptides = self.get_all_peptides()

        #initialize to zero coverage for all residues
        self.coverage = numpy.zeros(len(self.sequence))

        if len(peptides)==0:
            print("No peptides imported into state", self.name)
        else:
            for n in range(len(self.sequence)):
                for p in peptides:
                    if n+1 in p.get_observable_residue_numbers():
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
                #print(p, p.get_observable_residue_numbers())
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

    def calculate_residue_incorporation(self, protection_factors, change_tp_deut=True):
        # The self.residue_incorporations dictionary holds the per-residue
        # deuteration level for a given protection factor for each dataset
        self.residue_incorporations = {}

        for d in self.data:

            timepoints = set([tp.time for tp in d.get_all_timepoints()])
            self.residue_incorporations[d] = tools.calculate_incorporation(numpy.ones(len(protection_factors))*d.get_intrinsic_rates(), protection_factors, timepoints)
            #print(protection_factors, self.residue_incorporations[d])
            d.sum_residue_incorporations(self.residue_incorporations[d])

        return self.residue_incorporations


    def change_single_residue_incorporation(self, residue_number, new_val, change_tp_deut=True):
        '''
        Upon changing a protection factor, update the 2D incorporation for that residue
        among all peptides in all datasets and the output_model. 

        @param residue_number - the residue number to be changed
        @param new_val - the new protection factor model value to change it to
        @param change_tp_deut - True if you want to change the data incorporation values
        '''
        old_val = self.output_model.model[residue_number-1]

        #print("PF", residue_number, new_val, old_val)
        # Change residue already bakes in the -1
        self.output_model.change_residue(residue_number, new_val)

        for d in self.data:
            delta = {}
            new_rate = d.intrinsic[residue_number-1] - self.output_model.model_protection_factors[residue_number-1]

            times = d.get_all_times() 

            for time in times:
                new_deut = tools.calculate_simple_deuterium_incorporation(new_rate, time)
                old_deut = self.residue_incorporations[d][residue_number][time]
                self.residue_incorporations[d][residue_number][time] = new_deut
                delta[time] = new_deut - old_deut
                #print("  Change r#", residue_number, "from PF=", old_pf, "to", self.output_model.model_protection_factors[residue_number-1], "Deuts:", old_deut, new_deut)
                #print(residue_number, new_val, self.output_model.pf_grids[residue_number-1][new_pf-1], time, new_deut, old_deut)

            if change_tp_deut:   
                for pep in d.get_peptides_with_residue(residue_number):
                    for tp in pep.get_timepoints():
                        old_inc = tp.model_deuteration
                        #tp.model_deuteration += delta[tp.time]
                        tp.model_deuteration += delta[tp.time]
                        #print("  >>", pep.sequence, new_val, tp.time, "|||", old_inc, delta[tp.time], tp.model_deuteration)

    def calculate_score(self, model, calc_incorporations=False):
        # @param model is a proposed model for the system

        protection_factors = self.output_model.convert_model_to_protection_factors(model)
        total_score = 0

        if calc_incorporations:
            self.calculate_residue_incorporation(protection_factors)

        # Second, pass these protection factors to each dataset object
        #for d in self.data:
        #    dataset_score = self.calculate_peptides_score(d.get_peptides(), protection_factors)
        #    total_score += dataset_score

        return self.calculate_peptides_score(self.get_all_peptides(), protection_factors)

    def get_observable_residue_numbers(self):
        if self.observable_residue_numbers == []:
            coverage = self.get_coverage()
            orn = []

            for i in range(len(coverage)):
                if coverage[i] > 0:
                    orn.append(i+1)

            self.observable_residue_numbers = orn

        return self.observable_residue_numbers

    def get_observed_residues(self):
        '''
        Returns a list of residues that are observed in any dataset
        '''
        resis = set()
        for s in self.sectors:
            #print(s, s.get_residues())
            for r in s.get_residues():
                resis.add(r)
        self.observed_residues = list(resis)
        return self.observed_residues

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

        self.calculate_residue_incorporation(self.output_model.model_protection_factors)
        self.calculate_peptides_score(self.get_all_peptides(), self.output_model.model_protection_factors)

        for d in self.data:
            d.collect_times()

    def sum_incorporations(self, incorp, residues, time):
        deut = 0
        for r in residues:
            deut += incorp[r][time]
        return deut

    def calculate_peptides_score(self, peptides, protection_factors):
        '''
        Will deprecate calculate_dataset_score. Given a list of peptides,
        calculate the score.  Useful for calculating changes that only affect
        a susbset of peptides.
        '''
        return self.scoring_function.evaluate(protection_factors, peptides)

    def consolidate_model_to_sectors(self, model):
        # Given a model, return a list of lists with 
        pass

    def collect_score(self):
        score=0
        for d in self.data:
            score += d.get_score()
        self.score = score
        return score

    def get_score(self):
        try:
            return self.score
        except:
            return -1

    def set_score(self, score):
        self.score = score

    def calculate_residue_information_content(self):
        '''
        The set of peptides and timepoints can be used to calculate the amount of
        information at each residue for each protection factor in the model.

        For each peptide, the information at Pf = P for a residue is defined as
        1/(number of observable amides) * Resolving_power(P)

        Where resolving power is the sum over all timepoints, tp, of: 
            10^[(k_ex/P)*exp(-10^(k_ex/P)*tp)
        where k_ex is the intrinsic rate of the amide at the dataset conditions (T, pH)
        
        residue_information is a 2D grid of length=# of residues and height of the number of grid
        points in the hdx_model
        '''

        residue_information = numpy.zeros((len(self.sequence), self.output_model.grid_size))
        protection_factors = self.output_model.pf_grids[0]
        for d in self.data:
            residue_information += d.get_residue_information_content(protection_factors)

        return residue_information



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
        self.score = 100000000

    def get_number_of_residues(self):
        return len(self.residues)

    def get_number_of_amides(self):
        amides = 0
        for i in self.residues:
            if self.state.sequence[i-1]!="P":
                amides += 1
        return amides

    def get_coverage(self):
        return len(self.peptide_ids)

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

    def set_model_pf(self, pf, update=True):
        self.model_pf = pf
        # Change the residue incorporation values in the model
        if update:
            self.state.change_single_residue_incorporation(self.number, self.model_pf, True)

    def get_log_kex(self):
        return self.model_pf + self.intrinsic


def setup_single_state(sequence, name):
    '''Simple function that sets up a System with one Macromolecule
    and one State.  Returns the State object.
    '''
    sys = System()
    mol = sys.add_macromolecule(sequence, name)
    return mol.get_apo_state()




