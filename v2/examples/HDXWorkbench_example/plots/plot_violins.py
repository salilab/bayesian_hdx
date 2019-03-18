import sys
sys.path.append( "../../../pyext/src" )
import plots
import numpy
import hdx_models
import scipy
import scipy.cluster
from scipy.cluster.vq import vq, whiten, kmeans

ftp = 10
ltp = 3600
num_bins=19

rate_bins = numpy.linspace(numpy.log10(-numpy.log(0.99)/ltp), numpy.log10(-numpy.log(0.01)/ftp), num_bins)

model_file = "../test_output/wb_output_3339095_500.dat"
model_files = ["../test_output/wb_output_apo_500.dat", "../test_output/wb_output_3339095_500.dat"]


d = numpy.loadtxt(model_file)

white = whiten(d)

print(kmeans(white, 2))


plots.plot_residue_rate_distributions(model_files, rate_bins=rate_bins, resrange=(1,250))

