## Example of a simulated system
## USAGE:
# From the directory where this file is located, run:
# python modeling.py

import sys
sys.path.append( "../../pyext/src" )
import scoring
import sampling
import system
import model
import data
import plots
import cProfile


####################
###  What kind of run do you want to do?

run_type = "sampling"   # Use this to do the full sampling/analysis
#run_type = "benchmark" # Use this to run a quick simulation to see how long a full sampling run will take

##########################################
###    File/Directory Setup 

outputdir = "./d50_s5_t2_n10000_b20_r1/" # output directory for the simulation results. 

##########################################
###    Experimental  Input Parameters 
inseq = "AAA"   # FASTA sequence
offset = 0                  # offset between fragment start/end values and FASTA sequence
sigma0 = 5                  # Estimate for experimental error in %D Units 
saturation = 1.0            # Deuterium saturation in experiment
percentD = True             # Is the data in percent D (True) or Deuterium units?
tp_time = 300
rep_deut = 50
###########################################
###    Simulation Parameters
 
num_exp_bins = 20         # Number of log(kex) values for sampling. 20 is generally sufficient. 
init = "random"           # How to initialize - either "random" or "enumerate". Enumerate is slower but sampling will converge faster
annealing_steps = 200     # steps per temperature in annealing - 100-200 sufficient
nsteps = 10000            # equilibrium steps. 5000 to 10000

num_best_models = 200     # Number of best models to consider for analysis


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Input Data

sys = system.System()
mol = sys.add_macromolecule(name="Test", sequence=inseq)
state = mol.get_apo_state()

c = data.Conditions()
d = data.Dataset("Test", c, sequence=inseq)
pep = d.create_peptide("AAA", start_residue=1)
pep.add_timepoints([tp_time])
tp = pep.get_timepoint_by_time(tp_time)
tp.add_replicate(rep_deut)


state.add_dataset(d)

output_model = model.ResidueGridModel(state, grid_size=num_exp_bins)
state.set_output_model(output_model)
sampler = sampling.MCSampler(state, sigma_sample_level="timepoint")

state.initialize()
print(d.calculate_observable_rate_bounds())
print output_model.pf_grids[2]


sampler.run(10000, 2, write=True)

#for i in range(nsteps):
#    print(i, sampler.run_one_step(10000), output_model.model, tp.get_sigma())

plots.plot_residue_rate_distributions(output_model.output_model_file, rate_bins=rate_bins)

exit()



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
