# Usage:  /package_home/setup.py python analyze_output.py
import scoring
import system
import analysis
import plots


##########################################
###    File/Directory Setup 
outputdir = "./testing" # output directory for the simulation results. 
outputdir2 = "./testing2" # output directory for the simulation results.

sequence="GMAEDMAADEVTAPPRKVLIISAGASHSVALLSGDIVCSWGRGEDGQLGHGDAEDRPSPTQLSALDGHQIVSVTCGADHTVAYSQSGMEVYSWGWGDFGRLGHGNSSDLFTPLPIKALHGIRIKQIACGDSHCLAVTMEGEVQSWGRNQNGQLGLGDTEDSLVPQKIQAFEGIRIKMVAAGAEHTAAVTEDGDLYGWGWGRYGNLGLGDRTDRLVPERVTSTGGEKMSMVACGWRHTISVSYSGALYTYGWSKYGQLGHGDLEDHLIPHKLEALSNSFISQISGGWRHTMALTSDGKLYGWGWNKFGQVGVGNNLDQCSPVQVRFPDDQKVVQVSCGWRHTLAVTERNNVFAWGRGTNGQLGIGESVDRNFPKIIEALSVDGASGQHIESSNIDPSSGKSWVSPAERYAVVPDETGLTDGSSKGNGGDISVPQTDVKRVRI"   # FASTA sequence
resrange = (100,200) # Residue range is a tuple in pdb numbering (starts at 1).
num_best_models = 500


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###   Analysis.  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Initialize System
sys = system.System(output_dir = None)
mol = sys.add_macromolecule(sequence, "UVR8")
state = mol.get_apo_state()

pof = analysis.ParseOutputFile(outputdir + "/models_scores_sigmas-Apo.dat", state)
pof2 = analysis.ParseOutputFile(outputdir2 + "/models_scores_sigmas-Apo.dat", state)

plots.plot_residue_protection_factors([pof, pof2], num_best_models=num_best_models, resrange=(100,200))

