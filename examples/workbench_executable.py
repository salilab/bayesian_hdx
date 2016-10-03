import sys
sys.path.append( "../../../bayesian_hdx/pyext/src" )
import input_data
import sampling
import system_setup
import hdx_models
import argparse


# Here, set the default values for variables that can be changed via command line options
inseq=""

# Arguments: HDXWorkbench file, output_dir, run_type, inseq
parser = argparse.ArgumentParser(description='Analysis of HDXWorkbench DHDX')
parser.add_argument('input_file', metavar='i', type=str,
                    help='Path to HDXWorkbench .csv file')
parser.add_argument('output_dir', metavar='o', type=str,
                    help='Output directory')
parser.add_argument('--benchmark', help='Run a short benchmark instead of full sampling run',
                      action="store_true")
parser.add_argument('--inseq', metavar='s', type=str,
                    help='Input sequence string', required=False)
parser.add_argument('--NSTEPS', metavar='n', type=str,
                    help='Number of steps', required=False)


args = parser.parse_args()
print args

exit()

####################
###  What kind of run do you want to do?

if benchmark:
    run_type = "benchmark" # Use this to run a quick simulation to see how long a full sampling run will tak
else:
    run_type = "simulation"   # Use this to do the full sampling/analysis

##########################################
###    File/Directory Setup 
###
###  output directory for this simulation.
outputdir = args.output_dir
###   
###  HDX data file location
workbench_file=args.input_file

f=open(workbench_file, "r")
for line in f.readlines():
    if line.split(",")[0]=="Experiment Protein Sequence":
        inseq=line.split(",")[1]
        continue

if inseq="":
    raise Exception("HDX Workbench file does not contain FASTA sequence. Please manually add the sequence to the command line using the flag -s or --inseq")


##########################################
###    Experimental  Input Parameters 
  
offset=0                  # offset between fragment start/end values and FASTA sequence.
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

##########################################
###     Input Parameters 
# offset between fragment start/end values and FASTA sequence
offset=0
# Estimate for experimental error in %D Units
sigma0=1    
# Number of exp     
num_exp_bins=20   

###############################
###   System Setup:
###############################


# Initialize model  (name, FASTA sequence, offset)
model=system_setup.HDXModel("ERa",
                            inseq,
                            offset=offset)

# Add data to model (model, filename)
input_data.HDXWorkbench(model, workbench_file)


#Initialize a sampling model for each state (Multiexponential in this case)

for state in model.states:
        hdxm = hdx_models.MultiExponentialModel(model = model,
                                         state = state,
                                         sigma=sigma0,
                                          init = "enumerate")
        '''
        for f in state.frags:
          noa = f.get_num_observable_amides()
          exp_seq = hdxm.exp_seq[f.start_res-1:f.end_res-1]
          #print f.start_res, f.end_res, len(hdxm.exp_seq), exp_seq
          for t in f.timepoints:
            avg, sd = t.get_avg_sd()
            s_1 = t.calculate_tp_score(exp_seq, hdxm.exp_grid, 1.0, noa)
            s_5 = t.calculate_tp_score(exp_seq, hdxm.exp_grid, 5.0, noa)
            s_10 = t.calculate_tp_score(exp_seq, hdxm.exp_grid, 10.0, noa)
            s_20 = t.calculate_tp_score(exp_seq, hdxm.exp_grid, 20.0, noa)
            print noa, t.time, avg, sd, s_1, s_5, s_10, s_20
        '''

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

