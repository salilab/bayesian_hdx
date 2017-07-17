"""
   Classes to handle the HDX data hierarchy
"""
from __future__ import print_function
#import hdx_models
#import analysis
import io
import numpy
import scipy
#import scipy.special
from numpy import linalg
import sys
from copy import deepcopy
import os.path


class System(object):
    """ Class defining a set of macromolecular objects
    within an HX Experiment
    """
    def __init__(self):
        self.macromolecules = []

    def add_molecule(self, sequence, name=None):
        """ add a Molecule to the experiment. Sequence can either be
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
            self.macromolecules.append(Macromolecule(name, sequence))



class Macromolecule(object):
    def __init__(self, name, sequence, initialize_apo=True):
        self.name = name
        self.sequence = sequence
        self.states = []
        if initilize_apo:
            self.add_state("Apo", None)

    def get_sequence(self):
        return self.sequence

    def calc_intrinsic_rates(self):
        self.intrinsic = tools.get_sequence_intrinsic_rates(self.sequence)

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
                return Exception(name + " is not in the list of states:" + str(s.get_name() for s in self.states]))

        else:
            if state_number > len(self.states):
                return Exception("There are only " + str(len(self.states)) + "for molecule" + self.name)

            return self.states[0]

    def get_states(self):
        return self.states



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
    def __init__(self, mol, name, perturbations=None):
        self.name = name
        self.macromolecule = mol
        self.perturbations = []
        self.sequence = self.macromolecule.get_sequence()

        if perturbations is not None:
            self.add_perturbations(perturbations)

        self.data = []
        self.model = None

        self.intrinsic = self.calculate_intrinsic_rates()

    def has_model(self):
        if self.model is None:
            return False
        else:
            return True

    def has_data(self):
        if len(self.data) == 0:
            return False
        else:
            return True

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
                raise Exception("Exiting at State.add_dataset. Peptide " + p.sequence "does not match")
        self.data.append(dataset)
        dataset.set_state(self)

    def peptide_sequence_consistency(self, peptide):
        '''
        Returns True if fragment sequence and start residue aligns 
        with macromolecule sequence

        Returns False with a warning if there is an inconsistency
        '''
        peptide_position = 0

        for residue_number in range(peptide.start_residue, frag.end_residue):
            peptide_position = residue_number - peptide.start_residue
            if peptide.sequence[peptide_position] != self.sequence[residue_number]:
                print("Peptide ", peptide.sequence, " does not match Sequence")
                print("Peptide position", peptide_position+1, " is ", peptide.sequence[peptide_position])
                print("Sequence position", residue_number, " is ", self.sequence[residue_number])
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
            print("No fragments imported into this state:", self.state_name)

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
        for d in datasets:
            peptides += d.get_peptides()

    def calculate_sectors(self, peptides=None):
        '''
        Given a list of peptides, returns a list of sector
        objects corresponding to the uniquely sampled
        segments in the sequence due to overlapping fragments
        '''

        if peptides is None:
            peptides = self.get_all_peptides()
        if len(peptides)==0:
            print("No peptides imported into this state:", self.name)

        # Overwrite the sector list each time
        self.sectors=[]

        residue_peptide_sets = dict({})
        self.residue_peptide_sets = residue_peptide_sets

        # First, determine all peptides that cover each residue
        # Store unique peptide ids in a dictionary of sets: residue_peptide_ids
        # with the dictionary key = residue number
        for i in range(1, 1+len(self.seq)):
            peptide_ids = set()
            for p in peptides:
                if i in p.get_residue_numbers():
                    peptide_ids.add(p.get_id())
            residue_peptide_sets[i] = peptide_ids

        # Second, collect all residues with the same
        # peptide_id sets in the sector_dictionary
        sector_number = 0
        self.sector_dictionary = dict({})

        for x in residue_peptide_sets: # x is residue number in PDB terms
            sector_residues = set(x)
            sector_peptide_ids = residue_peptide_sets[x]
            # Compare peptide ids for each residue.
            # Same peptide ids = same sector
            for j in residue_peptide_sets:
                peptide_ids = residue_peptide_sets[j]

                if peptide_ids == sector_peptide_ids:
                    sector_residues.add(j)

            sector_dictionary[sector_number] = sector_residues

            self.sectors.append(Sector(self, sector_residues, sector_peptide_ids, sector_number))

            # Remove all residues contained in this sector
            for n in sector_residues:
                residue_peptide_sets.pop(n)

        return self.sectors


    def get_sectors(self):
        return self.sectors

    def create_simulated_data(self):
        io.Dataset()

    def calculate_score(self):
        #Ensure that we have data
        if not self.has_model:
            raise Exception("System.calculate_score: Cannot calculate score without a model!!")
        if not self.has_data:
            raise Exception("System.calculate_score: Cannot calculate score without data!!")



 class Sector(object):
    """ Each Sector object represents a portion of the peptide sequence with differential
        overlap by the MS peptide fragments.
        @param state - The macromolecule state to which this sector belongs
        @param residues - list or set of residue numbers in the sector
        @param peptide_ids - list or set of peptide ids used in this sector
        @param sector_number - the unique numerical id for this sector

    """
    def __init__(self, state, residues, peptide_ids, sector_number):
        self.state = state
        self.residues = residues
        self.peptide_ids = peptide_ids
        self.id = sector_number
        self.num_amides = self.calc_total_amides(self.seq)

    def get_number_of_residues(self):
        return len(self.residues)

    def get_number_of_amides(self):
        amides = 0
        for i in self.residues:
            if self.state.seq[i-1]!="P":
                amides += 1
        return amides

    def get_coverage(self):
        return len(peptide_ids)  
        

class Residue(object):
    def __init__(self, state, number, residue_type):
        self.state = state
        self.number = number
        self.residue_type = residue_type

    def set_intrinsic_rate(self, rate):
        self.intrinsic = rate




