'''
Usage: python cytc_multi.py run# output_prefix

This script initiates a single MCMC chain using simulated annealing.

User parameters that can be modified are explained in the header.

This script uses an HDX Workbench file.  Other file formats require their own
method to import. See hxio.py for instructions on how to import different types of HDX data.
'''

from __future__ import print_function
import sys

# Either use bayesian_hdx/setup.sh or modify this path to point to the python code
#sys.path.append("/Users/saltzberg/salilab/repositories/bayesian_hdx/pyext/src")

import scoring
import sampling
import system
import model
import hxio
import plots
import numpy
import tools


##########################################
###   Experimental Data Input Parameters

inseq = "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"
offset = 0                  # offset between fragment start/end values and FASTA sequence.
datafile = "../data/pH7.5.dat" # Experimental data file (HDXWorkbench)

plot_data_info = True      # Plot peptide overlap and data information figure

##########################################
###   Priors

# Create MD Prior
md_prior_file = "../data/cytc_md.dat" # MD estimated protection factor values
md_prior_scale = 0.2                  # Scale value for MD prior - recommended 0.1-0.2
pf_list = tools.open_md_prior_file(md_prior_file, inseq)
MD_prior = scoring.ResiduePfPrior(pf_list, scale=md_prior_scale)

# Create Natural Abundance Prior
nat_abun_prior_scale = 0.1            # Natural abundance prior scale - recommended 0.01-0.2
NA_prior = scoring.ProtectionFactorNaturalAbundancePrior(prior_type="Gaussian", scale=nat_abun_prior_scale)

# Add all priors to the prior_list
prior_list = [MD_prior, NA_prior]

plot_priors = True    # Plot prior distribution figure

###########################################
###   Simulation Parameters

num_exp_bins = 100      # Number of log(kex) values in grid - 50-100 recommended
annealing_steps = 50    # Number of temperature steps in exponential decay - 20-100 recommended
steps_per_anneal = 5    # Steps per temperature in annealing - 5-10 recommended
nsteps = 50             # equilibrium steps per anneal - 2000 to 10000 recommended
Ns = 2                  # Number of simulated anneals per chain

###########################################
###   Advanced Sampling Parameters
Tmax = 100.0              # Annealing maximum temp
Tmin = 2.0                # Annealining minumum temp and equilibration temp
anneal_pct_moves = 100    # % of residues to move during each annealing steps
equil_pct_moves = 50      # % of residues to move during each equilibration step
anneal_adjacency = int(0.4*num_exp_bins)  # Maximum move size (in grid units) for annealing
equil_adjacency = int(0.2*num_exp_bins)   # Maximum move size (in grid units) for equilibration


##-------------------------------------------
##   Below here the work



def run_sim(nrun, output_prefix):
    print("")
    ##########################################
    ###   File/Directory Setup  
    
    #  output directory for this simulation.
    outputdir = "./"+output_prefix+"_"+str(nrun)+"/"

    ###############################
    ###   System Setup:

    # Initialize model
    sys = system.System(output_dir=outputdir, noclobber=False)
    mol = sys.add_macromolecule(inseq, "CytC", initialize_apo=False)
    
    # Import datasets using hxio function
    datasets = hxio.import_HDXWorkbench(datafile, mol, offset=0)
    
    states = mol.get_states()

    print("Generated",len(states),"states of molecule",mol.name)
    
    if plot_data_info:
        protection_factors = numpy.arange(-2,14,0.1)
        for s in states:
            print("* Plotting peptide overlap and data information files for state", s.name)
            plots.plot_overlap_and_information(s, protection_factors, 
                                               outfile=outputdir+mol.name+"_"+s.get_name()+"_data-info.png", 
                                               figwidth_scale=2.0)
    ###############################
    ###  Scoring and Representation

    # Define the HDX representation model. 
    # ResidueGridModel models the protection factor at each residue along a finite grid of size num_exp_bins.
    output_models = []
    for s in states:
        
        output_models.append(s.set_output_model(model.ResidueGridModel(s, grid_size=num_exp_bins, sample_only_observed_residues=False)))

    ### Add Priors ###
    for s in states:
        sf = s.get_scoring_function()
        for p in prior_list:
            sf.add_prior(p)

    # Create priors plot
    if plot_priors:
        print("* Plotting prior probability distributions from", len(prior_list), "priors") 
        plots.plot_priors(sf, len(inseq), outfile=outputdir+"priors.png")

    # Initialize states/output and write out data in .hxd format.
    states[0].initialize()
    states[1].initialize()
    output = sys.get_output()
    output.write_datasets()

    output.initialize_output_model_file(states[0], output_models[0].pf_grids)
    output.initialize_output_model_file(states[1], output_models[1].pf_grids)

    ###############################
    ###  Sampling

    # Initilize MC Sampler
    sampler = sampling.MCSampler(sys)

    # Run a few simulated annealings
    for i in range(Ns): 
        sampler.pct_moves = anneal_pct_moves
        # Do an exponential temperature decay annealing from T=100 to T=2
        sampler.run_exponential_temperature_decay(tmax=Tmax, tmin=Tmin, 
                                            annealing_steps=annealing_steps, 
                                            steps_per_anneal=steps_per_anneal,   # 3-5 seems to work well. Number of steps to run at each 
                                            write=False,
                                            adjacency=anneal_adjacency) # How many grid values each residue can move. 10-25% of num_exp_bins works well

        # Set the adjacency and pct_moves for equilibrium sampling
        sampler.residue_sampler.set_adjacency(True, equil_adjacency)
        sampler.pct_moves = equil_pct_moves
        
        # Run burn-in steps. Do not write these models out.
        sampler.run(nsteps,Tmin, write=False)

        # Run equilibrium steps and write these to the output file
        sampler.run(nsteps,Tmin, write=True)


if __name__ == '__main__':
    run_sim(sys.argv[1], sys.argv[2])
