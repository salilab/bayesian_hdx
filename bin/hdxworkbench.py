import argparse

from scipy.stats import sem

import analysis
import hdx_models
import input_data
import plots
import sampling
import system_setup

inseq = ""
benchmark = False
# Arguments: HDXWorkbench file, output_dir, run_type, inseq
parser = argparse.ArgumentParser(description='Analysis of HDXWorkbench DHDX')
parser.add_argument('-w', help='HDXWorkbench CSV file', required=True)
parser.add_argument('--benchmark', help='Run a short benchmark instead of full sampling run',
                    action="store_true")
parser.add_argument('--inseq', metavar='s', type=str,
                    help='Input sequence string', required=False)
parser.add_argument('--nsteps', help='Equilibrium steps. 5000 to 10000. '
                                     'Default: 1000',
                    default=1000,
                    type=int,
                    required=False)
parser.add_argument('--init', help='How to initialize - either "random" or "enumerate". '
                                   'Enumerate is slower but sampling will converge faster. '
                                   'Default: enumerate',
                    default="enumerate",
                    required=False)
parser.add_argument('--num_exp_bins', help='Number of log(kex) values for sampling. 20 is generally sufficient. '
                                           'Default: 20',
                    default=20,
                    type=int,
                    required=False)
parser.add_argument('--mol_name', help='Molecule name', type=str,
                    required=True)
parser.add_argument('-o', '--outputdir', help='Output directory.', required=True)
parser.add_argument('--annealing_steps', help='Steps per temperature in annealing - 100-200 sufficient. '
                                              'Default: 20',
                    default=20,
                    type=int,
                    required=False)
parser.add_argument('--sigma0', help='Estimate for experimental error in %%D Units. '
                                     'Default: 5',
                    default=5,
                    type=float,
                    required=False)
parser.add_argument('--saturation', help='Deuterium saturation in experiment. '
                                         'Default: 1.0',
                    default=1.0,
                    type=float,
                    required=False)
parser.add_argument('--offset', help='Offset between fragment start/end values and FASTA sequence. '
                                     'Default: 0',
                    default=0,
                    type=float,
                    required=False)

args = parser.parse_args()

if args.inseq is not None:
    inseq = args.inseq

####################
###  What kind of run do you want to do?

if args.benchmark:
    run_type = "benchmark"  # Use this to run a quick simulation to see how long a full sampling run will tak
else:
    run_type = "sampling"  # Use this to do the full sampling/analysis

##########################################
###    File/Directory Setup
###
###  output directory for this simulation.
outputdir = args.outputdir
###
###  HDX data file location
workbench_file = args.w

f = open(workbench_file, "r")
for line in f.readlines():
    if line.split(",")[0] == "Offset":
        offset = int(line.split(",")[1].strip())
    if line.split(",")[0] == "Experiment Protein Sequence":
        inseq = line.split(",")[1]
    if line.split(",")[0] == "Deuterium solution concentration":
        saturation = float(line.split(",")[1].strip())
    if line.split(",")[0] == "Experiment name":
        name = line.split(",")[1].strip().replace(" ", "_")

if inseq == "":
    raise Exception(
        "HDX Workbench file does not contain FASTA sequence. Please manually add the sequence to the command line using the flag -s or --inseq")

###########################################
###    Simulation Parameters


# User-controlled variables
if args.benchmark:
    benchmark = True

num_exp_bins = args.num_exp_bins  # Number of log(kex) values for sampling. 20 is generally sufficient.
init = args.init  # How to initialize - either "random" or "enumerate". Enumerate is slower but sampling will converge faster
nsteps = args.nsteps  # equilibrium steps. 5000 to 10000
num_best_models = 1000000  # Number of best models to consider for analysis
outputdir = args.outputdir

# Non user controlled vbl - for now.
offset = args.offset
annealing_steps = args.annealing_steps  # steps per temperature in annealing - 100-200 sufficient
sigma0 = args.sigma0  # Estimate for experimental error in %D Units
saturation = args.saturation  # Deuterium saturation in experiment
percentD = True  # Is the data in percent D (True) or Deuterium units? - Always percentD for Workbench.
###############################
###   System Setup:
###############################


# Initialize model  (name, FASTA sequence, offset)
model = system_setup.HDXModel(args.mol_name,
                              inseq,
                              offset=offset)

# Add data to model (model, filename)
input_data.HDXWorkbench(model, workbench_file)

# Initialize a sampling model for each state (Multiexponential in this case)
for state in model.states:
    hdxm = hdx_models.MultiExponentialModel(model=model,
                                            state=state,
                                            sigma=sigma0,
                                            init=init)
###############################
###   Sampling:
###

# If benchmark is set to true, run a short simulation to estimate runtime
if run_type == "benchmark":
    sampling.benchmark(model, sample_sigma=True)
    exit()

#   Simulated Annealing macro runs high temperature dynamics and relaxes
#   to low temperature, followed by an equilibration run of "nsteps"

if run_type == "sampling":
    sampling.simulated_annealing(model, sigma=sigma0, equil_steps=nsteps, sample_sigma=True,
                                 annealing_steps=annealing_steps,
                                 outfile_prefix="Macro", outdir=outputdir)

bsm, scores = analysis.get_best_scoring_models(model.states[0].modelfile, model.states[0].scorefile,
                                               num_best_models=num_best_models,
                                               prefix=model.states[0].state_name, write_file=False)

num_best_models_list = [100, 100]
minstderr2 = 1.0E+34
for i in range(100, len(scores)):
    subset = scores[0:i]
    stderr = sem(subset)
    if stderr < minstderr2:
        minstderr2 = stderr
        num_best_models_list[0] = i
print('No. model: {} stderr: {}'.format(num_best_models_list[0], minstderr2))

bsm, scores = analysis.get_best_scoring_models(model.states[1].modelfile, model.states[1].scorefile,
                                               num_best_models=num_best_models,
                                               prefix=model.states[1].state_name, write_file=False)

minstderr2 = 1.0E+34
for i in range(100, len(scores)):
    subset = scores[0:i]
    stderr = sem(subset)
    if stderr < minstderr2:
        minstderr2 = stderr
        num_best_models_list[1] = i
print('No. model: {} stderr: {}'.format(num_best_models_list[1], minstderr2))

###############################
###   Analysis:

plots.plot_2state_fragment_avg_model_fits_num_model_per_state(model.states[0], model.states[1],
                                                              sig=5.0,
                                                              num_best_models=num_best_models_list,
                                                              write_file=False,
                                                              outdir=outputdir,
                                                              show_plot=False)

plots.plot_apo_lig_dhdx(model, show_plot=False, save_plot=True,
                        outfile="dhdx.png",
                        outdir=outputdir,
                        noclobber=False)

plots.plot_fragment_chi_values(model.states[0], sig="model", outdir=outputdir, show_plot=False)
plots.plot_fragment_chi_values(model.states[1], sig="model", outdir=outputdir, show_plot=False)
