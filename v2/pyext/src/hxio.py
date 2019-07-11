"""@namespace IMP.hdx.io
   Utility Classes for reading various HDX file formats
"""

from __future__ import print_function
from data import Conditions, Dataset
import system
import re
import os
import tools
from itertools import groupby



def add_peptide_to_dataset(dataset, peptide_sequence, start_residue, charge_state=None):
    '''
    if peptide is in a dataset, returns that peptide, otherwise, creates a new peptide
    and returns that
    '''
    if not dataset.is_this_peptide_in_dataset(peptide_sequence, start_residue, charge_state):
        print(peptide_sequence, start_residue, charge_state, dataset.is_this_peptide_in_dataset(peptide_sequence, start_residue, charge_state))
        new_peptide = dataset.create_peptide(sequence=peptide_sequence, start_residue=start_residue, charge_state=charge_state)
        #if new_peptide is not None:
        #    print("Peptide ", peptide_sequence, "chg=",charge_state, "created for state", state_name)
        #else:
        #    print("Skipping this peptide")
                    # If it is there...grab that peptide
    else:   
        new_peptide = next((pep for pep in dataset.get_peptides() if pep.sequence==peptide_sequence and pep.start_residue==start_residue and pep.charge_state==charge_state), None)

    return new_peptide

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


def import_MSStudio(infile, state=None):
    """ Function for converting a MS Studio file
    into a Dataset object
    """
    conditions = Conditions()

    data = Dataset(conditions = conditions, input_file=infile)


def import_HXcolumns(infile, sequence, name="Data", percentD=False, conditions=None, error_estimate=5.0, n_fastamides=2, offset=0):
    """ Function for converting a set of HX columns intoa Dataset object
    """

    if type(conditions) is not Conditions:
        print("Standard Conditions used.  Please modify these in the script if you are not at 283K and pH=7")
        conditions = Conditions()

    f = open(infile,"r")
    line = f.readline()

    dataset = Dataset(name=name, sequence=sequence, conditions=conditions, error_estimate=error_estimate, offset=offset, input_file=infile)

    column_headers = line.rstrip().split(",")  # Pre-set to None.  See #XXXX

    for line in f.readlines():


        fields = line.split(",")
        sequence = fields[column_headers.index("peptide_seq")]
        start_res = int(fields[column_headers.index("start_res")]) + offset
        time = float(fields[column_headers.index("time")])
        deut = float(fields[column_headers.index("D_inc")])
        #print("-----", start_res, column_headers.index("start_res"), offset, line)
        if not percentD:
            deut = deut / tools.calculate_number_of_observable_amides(sequence, n_fastamides) * 100
        score = float(fields[column_headers.index("score")])

        #print("IO", sequence, start_res)

        new_peptide = dataset.create_peptide(sequence, start_res)

        if new_peptide is not None:
            # If the time is 0.0, that's weird.  Ignore that.
            if time==0.0:
                continue

            if deut < 105:  # XXX Hard coded cap on %D value.  Not ideal.
                new_timepoint=True

                # If timepoint is already in fragment, grab that timepoint
                if time in [tp.time for tp in new_peptide.get_timepoints()]:
                    tp = new_peptide.get_timepoint_by_time(time)
                #If not, add a new timepoint at this time.  What to put for sigma??
                else:
                    tp = new_peptide.add_timepoint(time)
     
                # add the deuteration value as a replicate.
                # Any other replicate information from the file should be added at this step.
                tp.add_replicate(deut, score=score)

            new_peptide.add_timepoint(time)

    dataset.calculate_observable_rate_bounds()
    dataset.calculate_observable_protection_factors()

    return dataset

def import_json(infile, name=""):
    # Imports a json written by Dataset.write_to_file()
    import json
    print(infile)
    with open(infile) as json_data:
        d = json.load(json_data)

    cond_dict = d["conditions"]
    conditions = Conditions()
    for a in cond_dict:
        conditions.a = cond_dict[a]

    dataset = Dataset(d["name"], conditions, d["sequence"],
        input_file = d["raw_data_file"],
        error_estimate = d["error_estimate"],
        offset = d["offset"],
        number_fast_exchanging_amides = d["num_fast_exchanging_amides"],
        percent_deuterium = d["percent_deuterium"])

    pep_dict = d['peptides']

    for pep in pep_dict:
        p = pep_dict[str(pep)]
        peptide = dataset.create_peptide(p["sequence"], p["start_residue"], 
                            peptide_id=pep, 
                            charge_state=p["charge_state"], 
                            retention_time=p["retention_time"],
                            sigma=p["sigma"])
        tp_dict = p["timepoints"]
        for tp in tp_dict:
            timepoint = peptide.add_timepoint(tp_dict[str(tp)]["time"], tp_dict[str(tp)]["sigma"])
            reps = tp_dict[str(tp)]["replicates"]
            for rep in reps:
                timepoint.add_replicate(reps[str(rep)]["deut"], experiment_id=reps[str(rep)]["experiment_id"], 
                                score=reps[str(rep)]["reliability"], rt=reps[str(rep)]["rt"])

    return dataset


def import_HDXWorkbench(infile, macromolecule=None, name="Data", sequence=None, error_estimate=5.0, n_fastamides=2, offset=0,
                    max_t=36000):
    '''
    HDXWorbench files are different because they contain information from two experiments.
    They recieve a macromolecule rather than a single state and will create
    each state within that macromolecule

    Will also create a list of two Dataset objects if no macromolecule is passed
    '''
    # set up the dataset conditions
    conditions = Conditions()

    if macromolecule is not None and sequence is None:
        sequence = macromolecule.get_sequence()

    column_headers = None  # Pre-set to None.  See #XXXX

    f = open(infile,"r")
    line = f.readline()
    # First, let's get the header information from the workbench file
    # All lines are param, value pairs separated by a comma
    while len(line.split(','))==2:
        param=line.split(',')[0]
        value=line.split(',')[1]
        if param=="Used fully deuterated controls":
            if value=="true" or value=="True" or value=="T":
                conditions.fully_deuterated = True
            elif value=="false" or value=="False" or value=="F":
                conditions.fully_deuterated=False
        if param=="Temperature":
            conditions.temperature = float(value)+273.15
        if param=="Offset":
            conditions.offset = int(value)
        if param=="Deuterium solution concentration":
            conditions.saturation = float(value)
        if param=="Recovery estimate":
            conditions.recovery = float(value)
        if param=="Experiment Protein Sequence":
            file_sequence = str(value)
            if file_sequence != sequence and sequence is not None:
                print("WARNING: Sequence in HDXWorkbench file does not match inputted sequence")
            sequence = file_sequence.strip()
        if param=="Experiment name":
            name = str(value).replace(" ","_")
        line=f.readline()
    states = set()
    for line in f:

        # For nothing, just keep going
        if len(line.split(','))==0:
            continue

        #XXXX - the first line beginning with "peptide" is the list of column headers.
        if line.split(',')[0]=="peptide" and column_headers is None:
            column_headers = line.split(",")
            continue

        # This is a consolidated fragment entry. Just grab the state names for now 
        # (only two states in current format).
        if len(line.split(',')) >= 1 and column_headers != None and line.split(',')[column_headers.index("percentd_replicate")]=="NA" and len(states)==0:
            
            
            states.add(line.split(',')[column_headers.index("sample1_name")].replace(" ","_"))
            states.add(line.split(',')[column_headers.index("sample2_name")].replace(" ","_"))

            # Now create the different dataset
            datasets=[]
            for s in states:
                d = Dataset(name=s, sequence=sequence, conditions=conditions, error_estimate=error_estimate, input_file=infile, percent_deuterium=True)
                datasets.append(d)

        #############################################################
        #  This is replicate data. 
        #############################################################
        if len(line.split(',')) >= 1 and column_headers != None and line.split(',')[column_headers.index("percentd_replicate")]!="NA":
            # First see if it was discarded:
            discarded = line.split(',')[column_headers.index("discarded_replicate")]
            if discarded==True or discarded=="T" or discarded == "true" or discarded == "True" or discarded == "TRUE":
                continue

            # If not, get the peptide seq / start 
            # ********** Also, potentially other parameters *************
            fields = line.split(",")
            charge_state = int(fields[column_headers.index("charge")])
            peptide_sequence = fields[column_headers.index("peptide")]
            start_residue = int(fields[column_headers.index("start")])
            state_name = fields[column_headers.index("sample")].replace(" ","_")
            replicate_id = int(fields[column_headers.index("replicate")])

            try:
                replicate_score = float(fields[column_headers.index("score_replicate")])
            except ValueError:
                try:
                    replicate_score = float(fields[column_headers.index("score_replicate")+1])
                    #print("DHJSDNSLCNSIOACJSNCSAJ", line)
                except ValueError:
                    replicate_score = 0
                    #print("DHJSDNSLCNSIOACJSNCSAJ", 0, line)

            # retention time is a tuple with the start and end
            replicate_rt = (float(fields[column_headers.index("rt_start_replicate")]), float(line.split(',')[column_headers.index("rt_end_replicate")]))

            for data in datasets:
                if state_name == data.name:
                    # If the peptide is not there...create the peptide
                    is_pep_in, new_peptide = data.is_this_peptide_in_dataset(peptide_sequence, start_residue, charge_state)
                    
                    if new_peptide is None:
                        new_peptide = data.create_peptide(sequence=peptide_sequence, start_residue=start_residue, charge_state=charge_state)
            #print("XX", state_name, replicate_id, peptide_sequence, charge_state, new_peptide)            
            #print("Replicate", new_peptide)
            # Now, once the we know the peptide is there (or has just been created), add the data
            if new_peptide is not None:
                time = float(line.split(',')[column_headers.index("timepoint")].replace("s",""))

                # If the time is 0.0, that's weird.  Ignore that.
                if time==0.0 or time > max_t:
                    continue

                deut=float(line.split(',')[column_headers.index("percentd_replicate")])

                if deut < 155:  # XXX Hard coded cap on %D value.  Not ideal.
                    new_timepoint=True

                    # If timepoint is already in fragment, grab that timepoint
                    if time in [tp.time for tp in new_peptide.get_timepoints()]:
                        tp = new_peptide.get_timepoint_by_time(time)
                    #If not, add a new timepoint at this time.  What to put for sigma??
                    else:
                        tp = new_peptide.add_timepoint(time)
         
                    # add the deuteration value as a replicate.
                    # Any other replicate information from the file should be added at this step.
                    tp.add_replicate(deut, experiment_id=replicate_id, score=replicate_score, rt=replicate_rt)

                new_peptide.add_timepoint(time)

    if macromolecule is not None:
        for d in datasets:
            d.calculate_observable_rate_bounds()
            d.calculate_observable_protection_factors()
            s = macromolecule.add_state(d.name)
            s.add_dataset(d)

    return datasets


def import_csv(infile, state=None):
    """ Function for converting a csv file
    into a Dataset object
    """
    data = Dataset()

def import_Waters(infile, state=None):
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

class Output(object):
    '''
    Class for reading and writing output files from the method
    '''
    def __init__(self, system, output_directory, prefix="", noclobber=True):
        self.system = system
        self.output_directory = output_directory
        self.noclobber = noclobber
        if os.path.exists(output_directory):
            if noclobber:
                raise Exception("Output directory " + output_directory + " already exists. Exiting.")
            #else:
            #    raise Warning("Output directory " + output_directory + " already exists. Results will be overwritten")
        else:
            os.makedirs(output_directory)

    def change_output_directory(self, dir):
        self.output_directory = dir
        if os.path.exists(self.output_directory):
            if self.noclobber:
                raise Exception("Output directory " + output_directory + " already exists. Exiting.")
            #else:
            #    raise Warning("Output directory " + output_directory + " already exists. Results will be overwritten")
        else:
            os.mkdir(self.output_directory)

    def initialize_output_model_file(self, state, grid):
        # The output model file should contain:
        # list of values that relate grid values to 
        # sigma_sampler peptide
        # 0 0 2 3 4 4 9 4 9 | score || 

        output_model_file = self.output_directory+"/models_scores_sigmas-" + state.name +".dat"

        f = open(output_model_file, "w")
        f.write("&&*&%&*&& HX Modeling Output File &&*&%&*&&\n")
        f.write("Molecule_Name : " + state.macromolecule.name + "\n")
        f.write("State : " + state.name + "\n")
        f.write("\n")
        # Datasets
        f.write("### Datasets used | error_estimate\n")
        for d in state.data:
            if d.raw_data_file is not None:
                print(state.macromolecule.name, state.name, d.name, d.raw_data_file, d.sigma0)
                f.write("# " +"/datasets/"+ state.macromolecule.name +"_" + state.name + "_" + d.name + ".hxd" + " | "+ d.raw_data_file + " | " + str(d.sigma0) + "\n")
        f.write("\n")
        # Sectors
        f.write("@@@ Sector Residues and Coverage\n")
        sec_string = "@ "
        for s in state.calculate_sectors():
            resis = list(s.get_residues())
            resis.sort()
            for res in resis:
                sec_string+= str(res)+" "
            sec_string += "| "

        cov_string = "@C "
        for s in state.get_coverage():
            cov_string+=str(int(s))+ " "

        f.write(sec_string +"\n")
        f.write(cov_string +"\n")
        f.write("\n")
        f.write("grid_size : " + str(state.output_model.grid_size) + "\n")
        # Protection Factor Grids
        f.write("$$$ Residue PF Grids\n")     
        f.write("$ Residue_number | protection_factors\n") 
        pf_grids = state.output_model.pf_grids
        for i in state.get_observed_residues():
            string = "$ "+str(i)+" | "
            for j in state.output_model.pf_grids[i-1]:
                string += str(j)+" "
            f.write(string + "\n")
        f.write("\n")

        f.close()

    def write_datasets(self):
        print("DS", self.output_directory)
        outdir = self.output_directory + "/datasets/"
        try:
            os.mkdir(outdir)
        except:
            pass
        for mmols in self.system.get_macromolecules():
            for state in mmols.get_states():
                for dataset in state.get_datasets():
                    outfile = outdir + mmols.name +"_" + state.name + "_" + dataset.name + ".hxd"
                    dataset.write_to_file(outfile)


    def write_model(self, state, model, score, acceptance, sigmas=True):

        output_model_file = self.output_directory+"/models_scores_sigmas-" + state.name +".dat"
        f = open(output_model_file, "a")
        outstring = "> "
        for i in model:
            outstring += str(int(i))+" "
        outstring += "| " + str(score) + " | " + str(acceptance) + " || "

        if sigmas:
            for d in state.data:
                for tp in d.get_all_timepoints():
                    outstring += str(tp.get_sigma()) + " "

        f.write(outstring + "\n")

        f.close()


    def get_output_file(self, state):
        return self.output_directory+"/models_scores_sigmas-" + state.name +".dat"

    def write_model_to_file(self, f, state, model, score, acceptance, sigmas=True):

        outstring = "> "
        for i in model:
            outstring += str(int(i))+" "
        outstring += "| " + str(score) + " | " + str(acceptance) + " || "

        if sigmas:
            for d in state.data:
                for tp in d.get_all_timepoints():
                    outstring += str(tp.get_sigma()) + " "

        f.write(outstring + "\n")


