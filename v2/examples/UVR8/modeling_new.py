# Usage:  /package_home/setup.py python analyze_output.py
import sys
import cProfile
sys.path.append( "../../pyext/src" )
import scoring
import sampling
import system
import model
import hxio
import tools
import analysis
import plots

##########################################
###    File/Directory Setup 

outputdir="./XXtesting5000_15_20pct/" # output directory for the simulation results. 
 
#  HDX data file location
datadir="/Users/saltzberg/Dropbox/projects/salilab/bayesian_hdx_new/examples/UVR8/data/"
infile=datadir + "testdata_5.csv"

##########################################
###    Experimental Input Parameters 

sequence="GMAEDMAADEVTAPPRKVLIISAGASHSVALLSGDIVCSWGRGEDGQLGHGDAEDRPSPTQLSALDGHQIVSVTCGADHTVAYSQSGMEVYSWGWGDFGRLGHGNSSDLFTPLPIKALHGIRIKQIACGDSHCLAVTMEGEVQSWGRNQNGQLGLGDTEDSLVPQKIQAFEGIRIKMVAAGAEHTAAVTEDGDLYGWGWGRYGNLGLGDRTDRLVPERVTSTGGEKMSMVACGWRHTISVSYSGALYTYGWSKYGQLGHGDLEDHLIPHKLEALSNSFISQISGGWRHTMALTSDGKLYGWGWNKFGQVGVGNNLDQCSPVQVRFPDDQKVVQVSCGWRHTLAVTERNNVFAWGRGTNGQLGIGESVDRNFPKIIEALSVDGASGQHIESSNIDPSSGKSWVSPAERYAVVPDETGLTDGSSKGNGGDISVPQTDVKRVRI"   # FASTA sequence
offset = 1                  # offset between fragment start/end values and FASTA sequence
sigma0 = 3.0                # Estimate for experimental error in %D Units 
saturation = 1.0            # Deuterium saturation in experiment (from 0-1)
percentD = False            # Is the data in percent D (True) or Deuterium units (False)?

###########################################
###    Simulation Parameters
 
num_exp_bins = 71          # Number of log(kex) values for sampling. 15-20 is generally sufficient. 
init = "random"             # How to initialize - either "random" or "enumerate". Enumerate is slower but sampling will converge faster
nsteps = 5000              # equilibrium steps. 5000 to 10000

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###   Here the real work begins....
###   You should not have to change anything beneath this line.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############################################
###   System Setup:

# Initialize model
sys = system.System(output_dir = outputdir, noclobber=False)
mol = sys.add_macromolecule(sequence, "UVR8")
state = mol.get_apo_state()
state1 = mol.add_state(name="Apo2")
state2 = mol.add_state(name="Apo3")
# Alternatively, you can use this single line macro, which is equivalent
#state = system.setup_single_state(sequence, "UVR8", output_dir = "test_output")
# Import data
dataset = hxio.import_HXcolumns(infile,       # HX Columns input file
                          sequence,           # FASTA sequence string
                          name="Data",          # A string identifier for the dataset (optional)
                          percentD=False,     # Is the data in percent deuterium (True) or absolute deuterons (False)
                          conditions=None,    # A data.Conditions object specifying pH, Temp, etc...  None uses a standard set of conditions
                          error_estimate=5.0, # The initial estimate for experimental SD in % deuterium untis.
                          n_fastamides=2,     # Number of fast exchanging N-terminal amides
                          offset=offset)      # Numbering offset between data and FASTA file. (positive or negative integer)

# Add data to molecule state
state.add_dataset(dataset)
state1.add_dataset(dataset)
state2.add_dataset(dataset)

output_model = model.ResidueGridModel(state, grid_size=num_exp_bins)
state.set_output_model(output_model)
output_model1 = model.ResidueGridModel(state1, grid_size=num_exp_bins)
state1.set_output_model(output_model1)
output_model2 = model.ResidueGridModel(state2, grid_size=num_exp_bins)
state2.set_output_model(output_model2)

sampler = sampling.MCSampler(sys, pct_moves = 20)#, sigma_sample_level="timepoint")

sys.output.initialize_output_model_file(state, output_model.pf_grids)
sys.output.initialize_output_model_file(state1, output_model1.pf_grids)
sys.output.initialize_output_model_file(state2, output_model2.pf_grids)

sampler.run(nsteps, 2.0, write=True)

pof = analysis.ParseOutputFile(outputdir + "/models_scores_sigmas-Apo.dat", state)
pof2 = analysis.ParseOutputFile(outputdir + "/models_scores_sigmas-Apo2.dat", state1)
pof3 = analysis.ParseOutputFile(outputdir + "/models_scores_sigmas-Apo3.dat", state2)

#pof.calculate_random_sample_convergence()
#pof2.calculate_random_sample_convergence()

conv = analysis.Convergence(pof, pof2)

print(conv.total_score_pvalue_and_cohensd())
#print(conv.residue_pvalue_and_cohensd())

plots.plot_residue_protection_factors([pof, pof2, pof3], num_best_models=200, sort_sectors=True)

