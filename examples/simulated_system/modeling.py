## Example of a simulated system
## USAGE:
# From the directory where this file is located, run:
# python modeling.py

import sys
sys.path.append( "../../../bayesian_hdx/pyext/src" )
import input_data
import sampling
import system_setup
import hdx_models
import plots
import cProfile


####################
###  What kind of run do you want to do?

run_type = "sampling"   # Use this to do the full sampling/analysis
#run_type = "benchmark" # Use this to run a quick simulation to see how long a full sampling run will take

##########################################
###    File/Directory Setup 

outputdir="./output2/" # output directory for the simulation results. 
 
#  HDX data file location
datadir="/Users/saltzberg/salilab/bayesian_hdx/examples/simulated_system/data/"
input_file_apo=datadir + "apo_hdx.dat"
input_file_lig=datadir + "lig_hdx.dat"

##########################################
###    Experimental  Input Parameters 
inseq="AAMNSTQRLLVAGGA"   # FASTA sequence
offset=0                  # offset between fragment start/end values and FASTA sequence
sigma0=5                  # Estimate for experimental error in %D Units 
saturation = 1            # Deuterium saturation in experiment
percentD=True             # Is the data in percent D (True) or Deuterium units?

###########################################
###    Simulation Parameters
 
num_exp_bins=20           # Number of log(kex) values for sampling. 20 is generally sufficient. 
init = "enumerate"        # How to initialize - either "random" or "enumerate". Enumerate is slower but sampling will converge faster
annealing_steps=200       # steps per temperature in annealing - 100-200 sufficient
nsteps=20000              # equilibrium steps. 5000 to 10000

num_best_models=1000      # Number of best models to consider for analysis


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###   Here the real work begins....
###   You should not have to change anything beneath this line.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############################################
###   System Setup:

# Initialize model
model=system_setup.HDXModel("Simulated",     # Name for this system
                              inseq,         # Fasta sequence string
                              offset=offset) # Offset from FASTA sequence to HDX data

# Add data to model. 
input_data.HDXColumns(model, input_file_apo, "Apo", 
                        default_sigma=sigma0, 
                        offset=offset,        
                        temp=298.15,          # temperature of experiment - currently not used, so no need to change.
                        saturation=saturation, 
                        percentD=percentD)

input_data.HDXColumns(model, input_file_lig, "Lig", 
                        default_sigma=sigma0, 
                        offset=offset, 
                        temp=298.15, 
                        saturation=saturation, 
                        percentD=percentD)



# Initialize a calculation model for each state (Only MultiExponentialModel available now)

for state in model.states:
        hdxm = hdx_models.MultiExponentialModel(model = model,
                                                 state = state,
                                                 sigma = sigma0,
                                                 num_exp_bins = num_exp_bins,
                                                 init = init)

###############################
###   Sampling:
###

# If benchmark is set to true, run a short simulation to estimate runtime
if run_type=="benchmark":
    sampling.benchmark(model, sample_sigma=True)
    exit()

#   Simulated Annealing macro runs high temperature dynamics and relaxes
#   to low temperature, followed by an equilibration run of "nsteps"

if run_type=="sampling":
    sampling.simulated_annealing(model, sigma=sigma0, equil_steps=nsteps, sample_sigma=True, annealing_steps=annealing_steps, 
                            outfile_prefix="Macro", outdir=outputdir)

###############################
###   Analysis:

# Use this line if analyzing data post-run
# hdxm.import_model_deuteration_from_file(model.states[1].frags, model.states[1].modelfile)


plots.plot_2state_fragment_avg_model_fits(model.states[0], model.states[1], 
                                          sig=5.0, 
                                          num_best_models=num_best_models, 
                                          write_file=False, 
                                          outdir=outputdir, 
                                          show_plot=False)

plots.plot_apo_lig_dhdx(model, show_plot=True, save_plot=False, 
                                          outfile="dhdx.png", 
                                          outdir=outputdir, 
                                          noclobber=False)

plots.plot_fragment_chi_values(model.states[0], sig="model", outdir=outputdir, show_plot=True)
plots.plot_fragment_chi_values(model.states[1], sig="model", outdir=outputdir, show_plot=True)
