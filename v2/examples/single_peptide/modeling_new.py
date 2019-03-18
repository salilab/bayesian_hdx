## \example hdx/ERa/ERa_088074_modeling.py
import scoring
import sampling
import system
import model
import hxio
import tools
import argparse
import math
import cProfile


##### File/Directory Setup ######
outputdir="./output5"
datadir="./"
workbench_file=datadir + "single_frag.csv"

####### Input Parameters #######
offset=0        # offset between fragment start/end values and FASTA sequence
sigma0=5.0        # Estimate for slope of experimental error (slope of sd vs. %D plot)
num_exp_bins=10   # Number of bins
annealing_steps=200
inseq="YSMKSKNVVPLYDLL"

sequence="AGSSHHHHHHSSGLVPRGSNSLALSLTADQMVSALLDAEPPILYSEYDPTRPFSEASMMGLLTNLADRELVHMINWAKRVPGFVDLTLHDQVHLLESAWLEILMIGLVWRSMEHPGKLLFAPNLLLDRNQGKSVEGMVEIFDMLLATSSRFRMMNLQGEEFVCLKSIILLNSGVYTFLSSTLKSLEEKDHIHRVLDKITDTLIHLMAKAGLTLQQQHQRLAQLLLILSHIRHMSNKGMEHLYSMKSKNVVPLYDLLLEMLDAHRLHAPTSRGGASVEETDQSHLATAGSTSSHSLQK"


data  = {10   : [21.33928571428864, 15.702380952384537,  6.877976190478665],
         30   : [21.101190476192265,12.139880952382162, 23.738095238098282],
         60   : [29.78869047619391, 27.208333333334846, 17.913690476193878],
         300  : [32.130952380956586,20.264880952382892, 31.663690476193555],
         900  : [38.83035714286129, 30.360119047621147, 49.33035714286134],
         3600 : [51.97321428571825, 37.898809523810655, 41.9226190476218], }


data  = {10   : [6.31378],
         30   : [10.83015],
         60   : [13.55446],
         300  : [22.86113],
         900  : [31.88971],
         3600 : [42.51679] }

print "Hello"
###############################
###   System Setup:
###############################


# Initialize model  (name, FASTA sequence, offset)
sys = system.System(output_dir=outputdir, noclobber=False)
mol = sys.add_macromolecule(sequence, "Single", initialize_apo=False)

datasets = hxio.import_HDXWorkbench(workbench_file,       # Workbench input file
                          macromolecule=mol,
                          sequence=None,      # FASTA sequence string
                          sigma0=sigma0)

sys.get_output().write_datasets()

for state in mol.get_states():
    output_model = model.ResidueGridModel(state, grid_size=num_exp_bins)
    state.set_output_model(output_model)
    sys.output.initialize_output_model_file(state, output_model.pf_grids)

sampler = sampling.MCSampler(sys, sigma_sample_level="timepoint")
# First, run a short minimization step
sampler.run(50, 0.0001)


for dataset in datasets:
  for pep in dataset.get_peptides():
    for tp in pep.get_timepoints():
        #try:
        i = tp.get_replicates()[0]
        rep_score = -1*math.log(state.scoring_function.replicate_score(tp.get_model_deuteration()/pep.num_observable_amides*100, tp.get_replicates()[0].deut, tp.get_sigma()))
        print pep.sequence, tp.time, tp.get_model_deuteration()/pep.num_observable_amides*100, tp.get_replicates()[0].deut, tp.get_score(), rep_score, "|", tp.get_sigma(), -1*math.log(state.scoring_function.experimental_sigma_prior(tp.get_sigma(), sigma0))
        #except:
        #    pass


sampler.run(annealing_steps, 3)
#sampler.residue_sampler.set_adjacency(True, 4)
sampler.run(annealing_steps, 2)
#sampler.residue_sampler.set_adjacency(True, 3)
sampler.run(annealing_steps, 0.6)
# This temperature tends to sit around 15% MC acceptance rate, which seems to be good.
sampler.run(2000, 0.6, write=True)

