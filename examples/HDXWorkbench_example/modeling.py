import sys
sys.path.append( "../../../bayesian_hdx/pyext/src" )
import input_data
import sampling
import system_setup
import hdx_models
import plots
import argparse



####################
###  What kind of run do you want to do?

#run_type = "benchmark" # Use this to run a quick simulation to see how long a full sampling run will tak

benchmark=False   # Use this to do the full sampling/analysis

##########################################
###    File/Directory Setup 
#
#    HDX data file location
workbench_file="Data_Export_HDX_Workbench.csv"

f=open(workbench_file, "r")
for line in f.readlines():
    if line.split(",")[0]=="Experiment Protein Sequence":
        inseq=line.split(",")[1].strip()
        continue
if inseq=="":
    raise Exception("HDX Workbench file does not contain FASTA sequence. Please manually add the sequence to the command line using the flag -s or --inseq")


###  output directory for this simulation.
outputdir = "./test_output2/"
###   
##########################################
###    Experimental  Input Parameters 
  
offset=0                  # offset between fragment start/end values and FASTA sequence.
sigma0=5                  # Estimate for experimental error in %D Units 
saturation=1            # Deuterium saturation in experiment
percentD=True             # Is the data in percent D (True) or Deuterium units?

###########################################
###    Simulation Parameters
 
num_exp_bins=20           # Number of log(kex) values for sampling. 20 is generally sufficient. 
init = "enumerate"        # How to initialize - either "random" or "enumerate". Enumerate is slower but sampling will converge faster
annealing_steps=200       # steps per temperature in annealing - 100-200 sufficient

nsteps=50000               # equilibrium steps. 5000 to 10000

num_best_models=5000       # Number of best models to consider for analysis


###############################
###   System Setup:
###############################


# Initialize model  (name, FASTA sequence, offset)
model=system_setup.HDXModel("Sample",
                            inseq,
                            offset=offset)

# Add data to model (model, filename)
input_data.HDXWorkbench(model, workbench_file)


#Initialize a sampling model for each state (Multiexponential in this case)

for state in model.states:
        hdxm = hdx_models.MultiExponentialModel(model = model,
                                         state = state,
                                         sigma=sigma0,
                                          init = init)

###############################
###   Sampling:
###

# If benchmark is set to true, run a short simulation to estimate runtime
if benchmark:
    sampling.benchmark(model, sample_sigma=True)
    exit()

#   Simulated Annealing macro runs high temperature dynamics and relaxes
#   to low temperature, followed by an equilibration run of "nsteps"
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

plots.plot_apo_lig_dhdx(model, show_plot=True, save_plot=True, 
                                          outfile="dhdx.png", 
                                          outdir=outputdir, 
                                          noclobber=False)

plots.plot_fragment_chi_values(model.states[0], sig="model", outdir=outputdir, show_plot=True)
plots.plot_fragment_chi_values(model.states[1], sig="model", outdir=outputdir, show_plot=True)

