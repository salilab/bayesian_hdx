import os
import multiprocessing as mp
from subprocess import Popen
import sys
import cytc_multi
import glob
import numpy


running_script = "./cytc_multi.py"
excalibur_path = "/Users/saltzberg/salilab/repositories/bayesian_hdx/setup.sh"
output_prefix = "cytc_multi_500"

def run_one_simulation(i, output_prefix):
    Popen([excalibur_path, "python", running_script, str(i), output_prefix])

nruns = 4

manager = mp.Manager()

output = mp.Queue()

processes = [mp.Process(target=cytc_multi.run_sim, args=(x, output_prefix)) for x in range(nruns)]

for p in processes:
    p.start()

for p in processes:
    p.join()


sys.path.append("/Users/saltzberg/salilab/repositories/bayesian_hdx/pyext/src")
import analysis
import plots
psrfs = []

files = glob.glob("./"+output_prefix+"*/models_scores_sigmas-2+.dat")
oa = analysis.OutputAnalysis(files)
psrfs.append(numpy.array(oa.calculate_rhat()))


files = glob.glob("./"+output_prefix+"*/models_scores_sigmas-3+.dat")
oa = analysis.OutputAnalysis(files)
psrfs.append(numpy.array(oa.calculate_rhat()))
plots.plot_rhat(psrfs, show=True)