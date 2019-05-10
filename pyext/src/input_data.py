"""@namespace IMP.hdx.input
   Utility Classes for reading various HDX file formats
"""

from __future__ import print_function
import system_setup

class HDXWorkbench(object):
    """ Collects HDX data from a standard HDX workbench file.
        File header must consist of two comma separated values.
        Data column headers must be equal to those in the example data folder and comma-delimited
        Outputs a file with a list of fragments (protein_name.ligand.frags.dat) and stores timepoint info in hdx.representation.fragment classes
    """
    def __init__(self, model, infile, sigma0=1.0):
        self.temp=0
        self.theta=0.1
        self.file=infile

        state_frag=None

        f=open(infile,"r")
        line=f.readline()

        column_headers=None
        # Get Header values
        # All are 2 value lines
        while len(line.split(','))==2:
            self.get_header_value(model, line)
            line=f.readline()
            #print line

        for line in f:
            #print line
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

    def get_header_value(self, model, line):
        param=line.split(',')[0]
        value=line.split(',')[1]
        if param=="Used fully deuterated controls":
            if value=="true" or value=="True" or value=="T":
                self.fully_deuterated=True
            elif value=="false" or value=="False" or value=="F":
                self.fully_deuterated=False
        if param=="Temperature":
            self.temp=float(value)+273.15
        if param=="Offset":
            self.offset=int(value)
        if param=="Deuterium solution concentration":
            self.theta=float(value)
        if param=="Recovery estimate":
            self.recovery=float(value)
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
