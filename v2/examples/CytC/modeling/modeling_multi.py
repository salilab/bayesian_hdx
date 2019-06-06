import scoring
import sampling
import system
import model
import hxio
import tools
import math
import data
import analysis
import sys
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
outputdir = sys.argv[1]
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
sys = system.System(output_dir=outputdir, noclobber=False)
mol = sys.add_macromolecule(inseq, "CytC", initialize_apo=False)

# Import data
datasets = hxio.import_HDXWorkbench(workbench_file,       # Workbench input file
                          macromolecule=mol,
                          sequence=None,      # FASTA sequence string
                          error_estimate=sigma0) # The initial estimate for experimental SD in % deuterium untis.

#sys = system.System(output_dir=outputdir, noclobber=False)
#mol = sys.add_macromolecule(inseq, "CytC", initialize_apo=True)

state = mol.get_apo_state()

conditions = [data.Conditions(294, 6.5, 0.8), data.Conditions(294, 7.4, 0.8)]
datasets[0].conditions = conditions[1]
datasets[1].conditions = conditions[0]

ds = datasets[0]
state.add_dataset(ds)

#ds.calculate_intrinsic_rates()

#ds.calculate_observable_rate_bounds()
#ds.max_rate = ds.get_max_rate()
'''
for i in range(len(datasets)):
  #sigs+= d.get_all_tp_avg_sigmas()
  for p in datasets[i].get_peptides():
    print p.start_residue, p.sequence


  peps = datasets[i].get_peptides()
  #tp = peps[0].get_timepoints()[0]
  #peps[0].timepoints = [tp]
  #datasets[i].peptides = [peps[5]]
  datasets[i].conditions = conditions[i]
  datasets[i].calculate_intrinsic_rates()
  datasets[i].calculate_observable_rate_bounds()
  datasets[i].max_rate = datasets[i].get_max_rate()
  datasets[i].calculate_observable_protection_factors()
  print "Dataset", i, "::", len(datasets[i].get_all_timepoints()), "total timepoints", len(datasets[i].get_peptides())

  print datasets[i].intrinsic
'''

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
sampler.run(annealing_steps, 100, find_temperature=False)

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


'''
pof = analysis.ParseOutputFile(outputdir + "/models_scores_sigmas-CytC_pH_7.4.dat", mol.get_states()[0])
pof2 = analysis.ParseOutputFile(outputdir + "/models_scores_sigmas-CytC_pH_6.5.dat", mol.get_states()[1])

plots.plot_residue_protection_factors([pof, pof2], num_best_models=num_best_models, 
    resrange=resrange, true_vals=res_pfs, sort_sectors=True, outputdir=outputdir)
plots.plot_residue_protection_factors([pof, pof2], num_best_models=num_best_models, 
    resrange=resrange, true_vals=res_pfs, sort_sectors=False, outputdir=outputdir)

plots.plot_incorporation_curve_fits(pof, num_best_models, write_plots=True, single_plot=False, output_directory=outputdir)
plots.plot_incorporation_curve_fits(pof2, num_best_models, write_plots=True, single_plot=False, output_directory=outputdir)

'''