import sys
sys.path.append( "../../pyext/src" )
import plots


model_file = ["output_new/test_output.dat", "output_new/test_output2.dat"]

plots.plot_residue_rate_distributions(model_file, resrange=(1,400))

