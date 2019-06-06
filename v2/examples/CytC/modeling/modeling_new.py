import scoring
import sampling
import system
import model
import hxio
import tools
import argparse
import math
import cProfile
import data
import analysis

##########################################
###    File/Directory Setup 
#
#    HDX data file location
workbench_file="../data/cytc_workbench.csv"

# Quick hack for retrieving sequence
f=open(workbench_file, "r")
for line in f.readlines():
    if line.split(",")[0].strip()=="Experiment Protein Sequence":
        inseq=line.split(",")[1].strip()
        continue
if inseq=="":
    raise Exception("HDX Workbench file does not contain FASTA sequence. Please manually add the sequence to the command line using the flag -s or --inseq")

###  output directory for this simulation.
outputdir = "./test"
###   
##########################################
###    Experimental  Input Parameters 
  
offset = 0                  # offset between fragment start/end values and FASTA sequence.
sigma0 = 3                  # Estimate for experimental error in %D Units 

###########################################
###    Simulation Parameters
 
num_exp_bins=50           # Number of log(kex) values for sampling. 20 is generally sufficient. 
init = "enumerate"        # How to initialize - either "random" or "enumerate". Enumerate is slower but sampling will converge faster
annealing_steps=100       # steps per temperature in annealing - 100-200 sufficient

nsteps=1000              # equilibrium steps. 5000 to 10000

###############################
###   System Setup:
###############################


# Initialize model
sys0 = system.System(output_dir=outputdir, noclobber=False)
mol0 = sys0.add_macromolecule(inseq, "CytC", initialize_apo=False)

# Import data
datasets = hxio.import_HDXWorkbench(workbench_file,       # Workbench input file
                          macromolecule=mol0,
                          sequence=None,      # FASTA sequence string
                          error_estimate=sigma0) # The initial estimate for experimental SD in % deuterium untis.


sys = system.System(output_dir=outputdir, noclobber=False)
mol = sys.add_macromolecule(inseq, "CytC", initialize_apo=True)

state = mol.get_apo_state()

conditions = [data.Conditions(294, 6.5, 0.8), data.Conditions(294, 7.4, 0.8)]
datasets[0].conditions = conditions[1]
datasets[1].conditions = conditions[0]

ds = datasets[0]

#state.add_dataset(ds)

ds.calculate_intrinsic_rates()
ds.calculate_observable_rate_bounds()
ds.max_rate = ds.get_max_rate()

# Add data to molecule states and initialize models
for s in range(len(mol.get_states())):
    state = mol.get_states()[s]
    mod = model.ResidueGridModel(state, grid_size=num_exp_bins)
    output_model = state.set_output_model(mod)
    state.initialize()

    sys.output.initialize_output_model_file(state, output_model.pf_grids)

sys.get_output().write_datasets()

sampler = sampling.MCSampler(sys, sigma_sample_level=None, pct_moves=100)

# Slowly cool system
print "Anneal 1"
sampler.run(annealing_steps, 100, find_temperature=False)
print "Anneal 2"
# Optional modification of the number of grid points that the residue can shift (default is 5)
# (False) allows for any bin in the sampling space to be chosen at random

for i in range(5):
  sampler.run(annealing_steps, 50, find_temperature=False)
  #sampler.residue_sampler.set_adjacency(True, 3)
  sampler.run(annealing_steps, 20, find_temperature=False)
  sampler.run(annealing_steps, 10, find_temperature=False)
  sampler.residue_sampler.set_adjacency(True, 4)
  sampler.run(annealing_steps, 5, find_temperature=False)
  # This temperature tends to sit around 15% MC acceptance rate, which seems to be good.
  sampler.run(nsteps, 2, write=True, find_temperature=False)
  sampler.run(20, 0.1, write=True, find_temperature=False)



