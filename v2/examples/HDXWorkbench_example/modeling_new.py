import scoring
import sampling
import system
import model
import hxio
import tools
import argparse
import math
import cProfile

##########################################
###    File/Directory Setup 
#
#    HDX data file location
workbench_file="Data_Export_HDX_Workbench.csv"

f=open(workbench_file, "r")
for line in f.readlines():
    if line.split(",")[0].strip()=="Experiment Protein Sequence":
        inseq=line.split(",")[1].strip()
        continue
if inseq=="":
    raise Exception("HDX Workbench file does not contain FASTA sequence. Please manually add the sequence to the command line using the flag -s or --inseq")


###  output directory for this simulation.
outputdir = "./output_v2/"
###   
##########################################
###    Experimental  Input Parameters 
  
offset = 0                  # offset between fragment start/end values and FASTA sequence.
sigma0 = 5                  # Estimate for experimental error in %D Units 

###########################################
###    Simulation Parameters
 
num_exp_bins=20           # Number of log(kex) values for sampling. 20 is generally sufficient. 
init = "enumerate"        # How to initialize - either "random" or "enumerate". Enumerate is slower but sampling will converge faster
annealing_steps=20       # steps per temperature in annealing - 100-200 sufficient

nsteps=1000                # equilibrium steps. 5000 to 10000

###############################
###   System Setup:
###############################


# Initialize model
sys = system.System(output_dir=outputdir, noclobber=False)
mol = sys.add_macromolecule(inseq, "VDR", initialize_apo=False)

# Import data
datasets = hxio.import_HDXWorkbench(workbench_file,       # Workbench input file
                          macromolecule=mol,
                          sequence=None,      # FASTA sequence string
                          sigma0=sigma0) # The initial estimate for experimental SD in % deuterium untis.
sigs = []
for d in datasets:
  sigs+= d.get_all_tp_avg_sigmas()

# Add data to molecule states and initialize models
for s in range(len(mol.get_states())):
    state = mol.get_states()[s]
    state.add_dataset(datasets[s])
    output_model = state.set_output_model(model.ResidueGridModel(state, grid_size=num_exp_bins))
    sys.output.initialize_output_model_file(state, output_model.pf_grids)

sys.get_output().write_datasets()

sampler = sampling.MCSampler(sys, sigma_sample_level="timepoint")

# First, run a short minimization step
sampler.run(50, 0.0001, write=True)

'''
for dataset in datasets:
  for pep in dataset.get_peptides():
    for tp in pep.get_timepoints():
        #try:
        i = tp.get_replicates()[0]
        rep_score = -1*math.log(state.scoring_function.replicate_score(tp.get_model_deuteration()/pep.num_observable_amides*100, tp.get_replicates()[0].deut, tp.get_sigma()))
        print pep.sequence, tp.time, tp.get_model_deuteration()/pep.num_observable_amides*100, tp.get_replicates()[0].deut, tp.get_score(), rep_score, "|", tp.get_sigma(), -1*math.log(state.scoring_function.experimental_sigma_prior(tp.get_sigma(), sigma0))
        #except:
        #    pass
'''
# Slowly cool system
sampler.run(annealing_steps, 3)

# Optional modification of the number of grid points that the residue can shift (default is 5)
# (False) allows for any bin in the sampling space to be chosen at random
#sampler.residue_sampler.set_adjacency(True, 4)

sampler.run(annealing_steps, 2)
#sampler.residue_sampler.set_adjacency(True, 3)
sampler.run(annealing_steps, 1)


# This temperature tends to sit around 15% MC acceptance rate, which seems to be good.
sampler.run(nsteps, 1, write=True)


