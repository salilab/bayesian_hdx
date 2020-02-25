import sys
sys.path.append("/Users/saltzberg/salilab/repositories/bayesian_hdx/pyext/src")
import scoring
import sampling
import system
import model
import hxio
import plots
import numpy
import tools

##########################################
###    Experimental Data Input Parameters
offset = 0                  # offset between fragment start/end values and FASTA sequence.
md_prior_file = "../data/cytc_md.dat"
datafile = "data/pH7.5.dat"
md_prior_scale = 0.2
nat_abun_prior_scale = 0.1


###########################################
###    Simulation Parameters
num_exp_bins=100         # Number of log(kex) values for sampling. 20 is generally sufficient.
annealing_steps=5       # steps per temperature in annealing - 100-200 sufficient
steps_per_anneal=5
nsteps=5000              # equilibrium steps. 5000 to 10000
Ns = 2


def run_sim(nrun, output_prefix):
    #nrun = sys.argv[1]
    #output_prefix = sys.argv[2]
    ##########################################
    ###    File/Directory Setup
    #
    inseq="GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"

    ###  output directory for this simulation.
    outputdir = "./"+output_prefix+"_"+str(nrun)+"/"

    ###############################
    ###   System Setup:
    ###############################

    # Initialize model
    sys = system.System(output_dir=outputdir, noclobber=False)

    mol = sys.add_macromolecule(inseq, "Molecule", initialize_apo=False)
    datasets = hxio.import_HDXWorkbench(datafile, mol, offset=0)

    states = mol.get_states()

    ###############################
    #  Scoring and Representation
    ###############################

    # Define the HDX representation model. 
    # ResidueGridModel models the protection factor at each residue along a finite grid of size num_exp_bins.
    output_models = []
    output_models.append(states[0].set_output_model(model.ResidueGridModel(states[0], grid_size=num_exp_bins, sample_only_observed_residues=False)))
    output_models.append(states[1].set_output_model(model.ResidueGridModel(states[1], grid_size=num_exp_bins, sample_only_observed_residues=False)))

    # Define other aspects of the scoring function
    sf0 = states[0].get_scoring_function()
    sf1 = states[1].get_scoring_function()

    pf_list = tools.open_md_prior_file(md_prior_file, inseq)
    sf0.add_prior(scoring.ResiduePfPrior(pf_list, scale=md_prior_scale))
    sf1.add_prior(scoring.ResiduePfPrior(pf_list, scale=md_prior_scale))
    sf0.add_prior(scoring.ProtectionFactorNaturalAbundancePrior(prior_type="Gaussian", scale=nat_abun_prior_scale))
    sf1.add_prior(scoring.ProtectionFactorNaturalAbundancePrior(prior_type="Gaussian", scale=nat_abun_prior_scale))

    #plots.plot_priors(sf0, len(inseq))

    states[0].initialize()
    states[1].initialize()


    # Initialize output and write out data in .hxd format.
    output = sys.get_output()
    output.write_datasets()

    output.initialize_output_model_file(states[0], output_models[0].pf_grids)
    output.initialize_output_model_file(states[1], output_models[1].pf_grids)

    # Initialize MonteCarlo sampler
    sampler = sampling.MCSampler(sys)

    # Run a few simulated annealing steps.
    
    for i in range(Ns): 
        # Do an exponential temperature decay annealing from T=100 to T=2
        sampler.run_exponential_temperature_decay(tmax=100.0, tmin=2.0, 
                                            annealing_steps=annealing_steps, 
                                            steps_per_anneal=steps_per_anneal,   # 3-5 seems to work well. Number of steps to run at each 
                                            write=False,
                                            adjacency=0.2*num_exp_bins) # How many grid values each residue can move. 10-25% of num_exp_bins works well

        # Set the adjacency to be smaller now for equilibrium sampling
        sampler.residue_sampler.set_adjacency(True, 5)

        # Run burn-in step. Do not write these models out.
        sampler.run(nsteps,2.0, write=False)

        # Run equilibrium simulations and write these to the output file
        sampler.run(nsteps,2.0, write=True)


if __name__ == '__main__':
    run_sim()