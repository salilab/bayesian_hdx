import sys
import argparse

inseq=""
benchmark=False
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
parser.add_argument('--NSTEPS', metavar='n', type=int,
                    help='Number of steps', required=False)
parser.add_argument('--path', metavar='p', type=str,
                    help='Path to python code (if not already in PYTHONPATH)', required=False)
parser.add_argument('--init', metavar='z', type=str,
                    help='Initialize model by "random" or "enumerate". Random is faster, enumeration is slower, but will speed up sampling convergence',
                    required=False)
parser.add_argument('--bins', metavar='b', type=int,
                    help='Integer number of rate constant bins to use for sampling. Default is 20.',
                    required=False)

args = parser.parse_args()

sys.path.append( args.path )
try:
  import input_data
except:
  raise Exception("bayesian_hdx code not in PYTHONPATH. Use --path 'path/to/code/pyext/src' ")
import sampling
import system_setup
import hdx_models
import plots

if args.inseq is not None:
    inseq=args.inseq


####################
###  What kind of run do you want to do?

if args.benchmark:
    run_type = "benchmark" # Use this to run a quick simulation to see how long a full sampling run will tak
else:
    run_type = "sampling"   # Use this to do the full sampling/analysis


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
    if line.split(",")[0]=="Offset":
        offset=int(line.split(",")[1].strip())
    if line.split(",")[0]=="Experiment Protein Sequence":
        inseq=line.split(",")[1]
    if line.split(",")[0]=="Deuterium solution concentration":
        saturation=float(line.split(",")[1].strip())
    if line.split(",")[0]=="Experiment name":
        name=line.split(",")[1].strip().replace(" ","_")  

if inseq=="":
    raise Exception("HDX Workbench file does not contain FASTA sequence. Please manually add the sequence to the command line using the flag -s or --inseq")

###########################################
###    Simulation Parameters
 

# User-controlled variables
if args.benchmark:
    benchmark = True

if args.bins is None:
    num_exp_bins = 20            # Number of log(kex) values for sampling. 20 is generally sufficient. 
else:
    num_exp_bins= args.bins 

if args.init is None:
  if benchmark:
      init = "random" 
  else:
      init = "enumerate"
else:
    init = args.init       # How to initialize - either "random" or "enumerate". Enumerate is slower but sampling will converge faster

if args.NSTEPS is None:
    nsteps=5000              # equilibrium steps. 5000 to 10000
else:
    nsteps=args.NSTEPS

num_best_models=1000      # Number of best models to consider for analysis


# Non user controlled vbl - for now. 
annealing_steps=200       # steps per temperature in annealing - 100-200 sufficient
sigma0=5                  # Estimate for experimental error in %D Units 
saturation = 1.0            # Deuterium saturation in experiment
percentD=True             # Is the data in percent D (True) or Deuterium units? - Always percentD for Workbench.
###############################
###   System Setup:
###############################


# Initialize model  (name, FASTA sequence, offset)
model = system_setup.HDXModel("name",
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

