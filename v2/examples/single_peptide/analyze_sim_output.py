# Usage:  /package_home/setup.py python analyze_output.py
import scoring
import system
import analysis
import plots


##########################################
###    File/Directory Setup 
outputdir = "./sim_output1" # output directory for the simulation results. 
outputdir2 = "./sim_output2" # output directory for the simulation results.
outputdir3 = "./sim_output3" # output directory for the simulation results.
outputdir4 = "./sim_output4" # output directory for the simulation results.
outputdir5 = "./sim_output5" # output directory for the simulation results.
outputdir_enum = "./sim_output_enum" # output directory for the simulation results.


sequence="GMAEDMAADEVTAPPRKVLIISAGASHSVALLSGDIVCSWGRGEDGQLGHGDAEDRPSPTQLSALDGHQIVSVTCGADHTVAYSQSGMEVYSWGWGDFGRLGHGNSSDLFTPLPIKALHGIRIKQIACGDSHCLAVTMEGEVQSWGRNQNGQLGLGDTEDSLVPQKIQAFEGIRIKMVAAGAEHTAAVTEDGDLYGWGWGRYGNLGLGDRTDRLVPERVTSTGGEKMSMVACGWRHTISVSYSGALYTYGWSKYGQLGHGDLEDHLIPHKLEALSNSFISQISGGWRHTMALTSDGKLYGWGWNKFGQVGVGNNLDQCSPVQVRFPDDQKVVQVSCGWRHTLAVTERNNVFAWGRGTNGQLGIGESVDRNFPKIIEALSVDGASGQHIESSNIDPSSGKSWVSPAERYAVVPDETGLTDGSSKGNGGDISVPQTDVKRVRI"   # FASTA sequence
resrange = (100,200) # Residue range is a tuple in pdb numbering (starts at 1).
num_best_models = 20


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###   Analysis.  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Initialize System
sys = system.System(output_dir = None)
mol = sys.add_macromolecule(sequence, "ERa")
state = mol.get_apo_state()
#mol.add_state("088074")

pof = analysis.ParseOutputFile(outputdir + "/models_scores_sigmas-Apo.dat", state)
pof2 = analysis.ParseOutputFile(outputdir2 + "/models_scores_sigmas-Apo.dat", state)
pof3 = analysis.ParseOutputFile(outputdir3 + "/models_scores_sigmas-Apo.dat", state)
pof4 = analysis.ParseOutputFile(outputdir4 + "/models_scores_sigmas-Apo.dat", state)
pof5 = analysis.ParseOutputFile(outputdir5 + "/models_scores_sigmas-Apo.dat", state)
pof_e = analysis.ParseOutputFile(outputdir_enum + "/models_scores_sigmas-Apo.dat", state)
plots.plot_residue_protection_factors([pof_e], num_best_models=num_best_models)


