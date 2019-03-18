## Example of a simulated system
## USAGE:
# From the directory where this file is located, run:
# python modeling.py

import sys
sys.path.append( "../../pyext/src" )
import input_data
import sampling
import system
import hdx_models
import model
import plots
import cProfile
import tools
import numpy
import random

####################
###  What kind of run do you want to do?

run_type = "sampling"   # Use this to do the full sampling/analysis
#run_type = "benchmark" # Use this to run a quick simulation to see how long a full sampling run will take

##########################################
###    File/Directory Setup 

outputdir = "./sample_output2/" # output directory for the simulation results. 

##########################################
###    Experimental  Input Parameters 
inseq = "AAMNSTQRLL"   # FASTA sequence
pfs = [0,0,-1.2,-3.4,-3.5,-3.8,-5.8,-5.8,-1.5,-1.5] # log rate
peptides = [(1,5),(2,5),(2,8),(4,8),(4,10),(6,10)]
timepoints = [10,30,90,600,3600]
offset = 0                  # offset between fragment start/end values and FASTA sequence
sigma0 = 5                  # Estimate for experimental error in %D Units 
saturation_est = 0.7        # Deuterium saturation in experiment
sat_pct_error = 0.05
percentD = True             # Is the data in percent D (True) or Deuterium units?
time_error = 0.02
fwd_exchange_ph = 7.4
fwd_exchange_t = 293

back_exchange_est = 0.35        # Back exchange for the entire system
back_exchange_ph = 3.0          # 
back_exchange_t = 277

###########################################
####


# System Setup
sys = system.System("test_system", noclobber=False)
mol = sys.add_macromolecule(inseq, name="Test")
state = mol.get_apo_state()

# HDX model setup
hmod = model.ResidueGridModel(state, 20, protection_factors=True)

##############
# Get rates for each residue

act_sat = numpy.random.normal(saturation_est, sat_pct_error*saturation_est)
act_timepoints = [numpy.random.normal(tp,tp*time_error) for tp in timepoints]



fwd_rates = [-1]
back_rates = [-1]
for i in range(1, len(inseq)):
  ra2='A'
  la2='A'
  if i == 1:
    la2="NT"
  if i == len(inseq)-1:
    ra2 = "CT"

  fwd_rates.append(tools.calc_intrinsic_rate(inseq[i-1], inseq[i], fwd_exchange_ph, fwd_exchange_t, La2=la2, Ra2=ra2, log=False)) 
  back_rates.append(tools.calc_intrinsic_rate(inseq[i-1], inseq[i], back_exchange_ph, back_exchange_t, La2=la2, Ra2=ra2, log=False))
  print i, inseq[i], fwd_rates[-1], back_rates[-1]

for pep in peptides:
  full_deut = []
  for t in act_timepoints:
    deut = 0
    for i in range(pep[0]+1,pep[1]):
      k_ex = fwd_rates[i] * 10**pfs[i]
      k_back = back_rates[i]
      deut += (1-numpy.exp(-k_ex*t))*act_sat - (1-numpy.exp(-k_back*50))*(1-act_sat)
    full_deut.append(deut)
  print pep, full_deut

print act_sat, act_timepoints

exit()

###########################################
###    Simulation Parameters
 
num_exp_bins = 20         # Number of log(kex) values for sampling. 20 is generally sufficient. 
init = "random"           # How to initialize - either "random" or "enumerate". Enumerate is slower but sampling will converge faster
annealing_steps = 200     # steps per temperature in annealing - 100-200 sufficient
nsteps = 20000            # equilibrium steps. 5000 to 10000

num_best_models = 200     # Number of best models to consider for analysis










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
    sampling.simulated_annealing(model, sigma=sigma0, equil_steps=nsteps, sample_sigma=False, annealing_steps=annealing_steps, 
                            outfile_prefix="Macro", outdir=outputdir)

###############################
###   Analysis:

# Use this line if analyzing data post-run
# hdxm.import_model_deuteration_from_file(model.states[1].frags, model.states[1].modelfile)

'''
plots.plot_2state_fragment_avg_model_fits(model.states[0], model.states[1], 
                                          sig=5.0, 
                                          num_best_models=num_best_models, 
                                          write_file=True, 
                                          outdir=outputdir, 
                                          show_plot=False)

for s in model.states:
    plots.plot_fragment_avg_model_fits(s, 
                                          sig=5.0, 
                                          num_best_models=num_best_models, 
                                          write_file=True, 
                                          outdir=outputdir, 
                                          show_plot=False)



plots.plot_apo_lig_dhdx(model, show_plot=True, save_plot=False, 
                                          outfile="dhdx.png", 
                                          outdir=outputdir, 
                                          noclobber=False)
'''
plots.plot_fragment_chi_values(model.states[0], sig="model", outdir=outputdir, show_plot=True)
plots.plot_fragment_chi_values(model.states[1], sig="model", outdir=outputdir, show_plot=True)
