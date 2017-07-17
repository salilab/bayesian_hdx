"""@namespace IMP.hdx.input
   Utility Classes for reading various HDX file formats
"""

from __future__ import print_function
import System
import re
from itertools import groupby


def read_fasta(fasta_file):
    """
    Given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_file)
    # Headers for all entries begin with ">"
    entries = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in entries:
        # drop the ">" and grab only the first word.
        # Allow for dashes, underscore and colons (PDB)
        fields = re.findall(r'[A-Za-z0-9-_:]+', header.next())

        # Some fasta files have the header form "">gi|fasta_id|gb"
        if fields[0] == 'gi' and fields[2]=='gb':
            header = fields[1]
        else:
            header = fields[0]

        # join all subsequent sequence lines into one string.
        seq = "".join(s.strip() for s in entries.next())
        yield header, seq


class Conditions(object):
    """ Class containing a set of conditions for an HDX experiment """
    def __init__(self, temperature, pH, deuterium_concentration):
        self.temp = 0
        self.theta = 0
        self.perturbation = None



class Peptide(object):
    """ Class that stores a list of Timepoint data for each experimental HDX peptide fragment.
        Each Peptide object is bound to a single State.
        Contains experiment-specific data
        @param sequence - sequence of the peptide
        @param start_residue - starting residue of the peptide
        @param peptide_id - a unique peptide ID
        @param sigma - an error 

    """
    def __init__(self, sequence, start_residue, peptide_id, sigma=1.0, charge_state=None, retention_time=None, pH=7.0, temp=283):
        self.id = peptide_id
        self.sequence = sequence
        self.timepoints = []
        self.start_residue = int(start_residue)
        self.num_observable_amides = self.calc_num_observable_amides(inseq)
        self.sigma = sigma
        self.charge_state = charge_state
        self.retention_time = retention_time
        self.conditions = conditions
        self.mass = tools.calc_mass(inseq)

    def set_state(self, state):
        self.state = state

    def get_state(self):
        return self.state

    def set_sigma(self, sigma):
        self.sigma = sigma

    def get_residue_numbers(self):
        #return a list of the residue numbers in the peptide
        return range(self.start_residue, self.start_residue + len(self.sequence))

    def get_id(self):
        return self.id

    def calculate_number_of_observable_amides(self, inseq):
        #number of observable amides is equal to fragment length - 2, minus remaining prolines
        num_amides = inseq.count('P',2) + inseq.count('p',2)
        return len(self.seq)-num_amides-2

    def get_num_observable_amides(self):
        return self.num_observable_amides

    def add_timepoint(self, time):
        tp = Timepoint(time)
        self.timepoints.append(tp)
        return tp

    def get_timepoints(self):
        return self.timepoints

    def get_chi_value(self, sig):
        '''
        Return chi score of model compared to experimental data
        '''
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

    def calculate_peptide_score_freq_grid(self, freq_grid, exp_grid, sig, save=False, force=True):
        score=0
        for tp in self.timepoints:
            tp.calc_model_deut(freq_grid, exp_grid, self.num_observable_amides)
            score += tp.calculate_tp_score(grid, exp_grid, sig, self.num_observable_amides, force_calc=True)
        return score        

    def calculate_peptide_score(self, grid, exp_grid, sig, save=False, force=False):
        score = 0     
        for tp in self.timepoints:
            tp.calc_model_deut(grid, exp_grid, self.num_observable_amides)
            score += tp.calculate_tp_score(grid, exp_grid, sig, self.num_observable_amides, force_calc=True)
        return score

    def set_timepoint_sigmas(self, sigma=None):
        if sigma = None:
            sigma = self.sigma
        for tp in self.timepoints:
            tp.sigma = sigma 

    def get_peptides(self):
        return self.peptides


class Dataset(object):
    """ Characterizes a single HDX dataset for a macromolecular system.
    Requires a Macromolecule object
    Consists of a list of Fragments with associated timepoints, along with experimental parameters
    Can contain a perturbation
    """   

    def __init__(self, state, conditions):
        self.state = state
        self.conditions = conditions
        self.peptides = []

    def get_peptides(self):
        return self.peptides

    def get_perturbation(self):
        return self.conditions.perturbation

    def set_perturbation(self, perturbation):
        self.conditions.set_perturbation()

    def ensure_unique_peptide_id(self, peptide, modify_if_not_unique=True):
        '''
        Given a new peptide, make sure the id of this peptide is not already
        in the peptide list
        '''

        id_list = [p.id for p in self.peptides]

        if peptide.id in id_list:
            if modify_if_not_unique:
                peptide.id = str(peptide.id) + "_X"

            else:
                raise Exception("Peptide id" + str(peptide.id) + "is already in list")

    def add_peptide(self, peptide):
        '''
        Takes a peptide object and adds it to the State peptide list
        Returns peptide object
        '''
        if self.peptide_sequence_consistency(peptide):
            self.peptides.append(peptide)
            peptide.aset_state(self)
            return self.peptides[-1]
        else:
            raise Exception("Exiting at Dataset.add_peptide")


    def create_peptide(self, sequence, start_residue, peptide_id=None, sigma=1.0, charge_state=None, retention_time=None):
        '''
        Manually creates a peptide object
        adds it to the end of State peptide list

        Returns the fragment object
        '''
        if peptide_id is None:
            peptide_id = str(len(self.peptides))
        new_peptide = Peptide(seq, start_residue, peptide_id, sigma, charge_state, retention_time)
        if self.peptide_sequence_consistency(new_peptide): 
            self.peptides.append(new_peptide)
            new_peptide.set_state(state)
            return new_peptide
        else:
            raise Exception("Exiting at Dataset.create_peptide")


def import_MSStudio(system, infile, state=None):
    """ Function for converting a MS Studio file
    into a Dataset object
    """
    self.infile = infile
    conditions = Conditions()
    data = Dataset()

    def get_conditions(self):
        return 0



def import_csv(system, infile, state=None):
    """ Function for converting a csv file
    into a Dataset object
    """
    data = Dataset()

def import_Waters(system, infile, state=None):
    """ Function for converting a Waters file
    into a Dataset object
    """
    data = Dataset()



class HDXWorkbench(object):
    """ Collects HDX data from a standard HDX workbench file.
        File header must consist of two comma separated values.
        Data column headers must be equal to those in the example data folder and comma-delimited
        Outputs a file with a list of fragments (protein_name.ligand.frags.dat) and stores timepoint info in hdx.representation.fragment classes
    """
    def __init__(self, system, infile, sigma0=1.0):
        self.system = system
        self.temp=0
        self.theta=0.1
        self.file=infile
        self.conditions = Conditions()

        f=open(infile,"r")
        line=f.readline()

        column_headers=None
        # Get Header values
        # All are 2 value lines
        while len(line.split(','))==2:
            self.get_header_value(line)
            line=f.readline()
            #print line

        for line in f:

            if len(line.split(','))==0:
                continue

            if line.split(',')[0]=="peptide" and column_headers is None:
                column_headers=self.get_column_headers(line)
                continue

            if len(line.split(',')) >= 1 and column_headers != None and line.split(',')[column_headers.index("percentd_replicate")]=="NA":
                state=[]
                #This is a consolidated fragment entry. Just grab the state names for now (only two states in current format) and add to model if not already present.
                state.append(line.split(',')[column_headers.index("sample1_name")].replace(" ","_"))
                state.append(line.split(',')[column_headers.index("sample2_name")].replace(" ","_"))

                for s_data in state:
                    add_state_to_model=True
                    for s_mod in model.states:
                        if s_mod.state_name==s_data:
                            add_state_to_model=False
                    if add_state_to_model==True:
                        model.add_state(s_data) # Need some mechanism to add mole_fraction_liganded if weak binding ligand. Current default is 100% binding
                continue

            if len(line.split(',')) >= 1 and column_headers != None and line.split(',')[column_headers.index("percentd_replicate")]!="NA":
                #This is replicate data. First see if it was discarded:
                discarded=line.split(',')[column_headers.index("discarded_replicate")]
                if discarded==True or discarded=="T" or discarded == "true" or discarded == "True" or discarded == "TRUE":
                    continue
                #Get fragment seq / start / end
                frag_seq=line.split(',')[column_headers.index("peptide")]
                start_res=int(line.split(',')[column_headers.index("start")])
                end_res=int(line.split(',')[column_headers.index("end")])
                state_name=line.split(',')[column_headers.index("sample")].replace(" ","_")

                #First, see if this fragment is already in the model states
                for s in model.states:
                    if s.state_name==state_name:
                        if self.is_this_fragment_in_state(s, frag_seq, start_res)==False:
                            state_frag=s.create_fragment(frag_seq, start_res, end_res)
                            if state_frag is not None:
                                print("Fragment ", frag_seq, "created for state", s.state_name)
                            else:
                                print("Skipping this fragment")
                        else:  #if it is not, 
                            state_frag=next((f for f in s.frags if f.seq==frag_seq and f.start_res==start_res), None)
                            #print "State not created", s.state_name, "created for fragment", frag_seq
                   # elif state_frag==None:
                    #    print("Error finding fragment",frag_seq," in State,", s.state_name)
                        #exit()
                if state_frag is not None:
                    self.add_timepoint_to_frag(state_frag, column_headers, line, default_sig = sigma0, empirical_sig=True)
                    #print state_frag.timepoints[-1].time, state_frag.timepoints[-1].replicates[-1].deut
        for s in model.states:
            s.get_sectors(s.frags)
            s.get_coverage(s.frags)

        # Make dataset objects:


    def get_header_value(self, line):
        param=line.split(',')[0]
        value=line.split(',')[1]
        if param=="Used fully deuterated controls":
            if value=="true" or value=="True" or value=="T":
                self.conditions.fully_deuterated=True
            elif value=="false" or value=="False" or value=="F":
                self.conditions.fully_deuterated=False
        if param=="Temperature":
            self.conditions.temp=float(value)+273.15
        if param=="Offset":
            self.offset=int(value)
        if param=="Deuterium solution concentration":
            self.conditions.theta=float(value)
        if param=="Recovery estimate":
            self.conditions.recovery=float(value)
        if param=="Experiment Protein Sequence":
            self.seq=str(value)
        if param=="Experiment name":
            self.name=str(value).replace(" ","_")

    def is_this_fragment_in_state(self, state, frag_seq, start_res):
        answer=False
        for f in state.frags:
            if frag_seq==f.seq and start_res==f.start_res:
                answer=True
        return answer

    def get_column_headers(self, line):
        return line.split(",")

    def add_timepoint_to_frag(self, frag, column_headers, line, default_sig=1.0, empirical_sig=False):
        time=float(line.split(',')[column_headers.index("timepoint")].replace("s",""))
        if time==0.0:
            return 0
        deut=float(line.split(',')[column_headers.index("percentd_replicate")])
        #if deut > 100:
            #print(frag.seq, time, deut)
        if deut < 105:  # XXX Hard coded cap on %D value
            new_timepoint=True
            # If timepoint is already in fragment, grab that timepoint
            for t in frag.timepoints:
                if t.time==time:
                    new_timepoint=False
                    tp=t
            #If not, add a new timepoint at this time.
            if new_timepoint==True:
                tp=frag.add_timepoint(time)
 
            # add the deuteration value.
            # Any other replicate information from the file should be added at this step.
            tp.add_replicate(deut)#, temp=self.temp, sat=self.theta, recovery=self.recovery)
            #print len(tp.replicates)
            if empirical_sig==True and len(tp.replicates) > 2:
                avg, sd = tp.get_avg_sd()
                #print(len(tp.replicates), sd)
                tp.sig = sd          
            else:
                tp.sig = default_sig

class HDXColumns(object):
    '''
    Class for inputting HDX data from a simple file for a single state 
    Takes as input a file with columns:
    # peptide_seq, start_res, end_res, time, D_inc
    '''
    def __init__(self, model, infile, state_name, default_sigma=1.0, offset=0, temp=298.15, saturation=1.0, percentD=False):
        self.file=infile
        self.temp=temp
        self.sat=saturation
        self.offset=offset

        s = model.add_state(state_name)


        f = open(infile,"r")

        f.readline()

        for line in f.readlines():
            #self.add_timepoint_to_frag(line, percentD, self.sat, self.offset)

            fields=line.split(",")
            frag_seq=str(fields[0].strip())
            start_res=int(fields[1].strip())
            end_res=int(fields[2].strip())
            time=float(fields[3].strip())
            deut=float(fields[4].strip())

            if self.is_this_fragment_in_state(s, frag_seq, start_res)==False:
                state_frag=s.create_fragment(frag_seq, start_res, end_res)
                if state_frag is not None:
                    print("Fragment", frag_seq, "created for state", s.state_name)
                else:
                    print("Skipping this fragment")
            else:  #if it is not, 
                state_frag=next((f for f in s.frags if f.seq==frag_seq and f.start_res==start_res), None)

            # If this is not percentD, we want to divide D incorporation by number of amides in fragment
            # and multiply by 100
            if not percentD:
                deut = deut / float(state_frag.get_num_observable_amides()) * 100

            if state_frag is not None:
                self.add_timepoint_to_frag(state_frag, time, deut, default_sigma)

        s.get_sectors(s.frags)
        s.get_coverage()

    def add_timepoint_to_frag(self, frag, time, deut, sigma):

        if deut < 105:  # XXX Hard coded cap on %D value
            new_timepoint=True
            #print(frag, type(frag))
            # If timepoint is already in fragment, grab that timepoint
            for t in frag.timepoints:
                if t.time==time:
                    new_timepoint=False
                    tp=t
            #If not, add a new timepoint at this time.
            if new_timepoint==True:
                tp=frag.add_timepoint(time)
 
            # add the deuteration value.
            # Any other replicate information from the file should be added at this step.
            tp.add_replicate(deut)#, temp=self.temp, sat=self.theta, recovery=self.recovery)
            #print len(tp.replicates)

            tp.sig = sigma


    def is_this_fragment_in_state(self, state, frag_seq, start_res):
        answer=False
        for f in state.frags:
            if frag_seq==f.seq and start_res==f.start_res:
                answer=True
        return answer



class OutputFile(object):
    ''' Class that reads and stores the elements of an modeling output file
    '''
    def __init__(self, infile):
        self.output_file = infile

        for l in self.output_file.readlines():
            if l[0] == "#":
                self.read_header(line)


    def read_header(self, line):
        fields = line[1:].strip().split()

        if fields[0] == "logk_grid":
            self.logk_grid = numpy.linspace(fields[1], fields[2], fields[3])
            return False
        elif fields[0] == "peptides":
            # Each peptide is a tuple

        elif fields[0]=="models":
            mstring=""
            for l in self.output_file.readlines():
                models+=l
            self.models = numpy.genfromtxt(models)






