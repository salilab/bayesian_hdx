## \example hdx/ERa/ERa_088074_modeling.py
import scoring
import sampling
import system
import model
import data
import tools
import argparse
import math
import cProfile
import analysis
import plots


##### File/Directory Setup ######
outputdir="./sim_output_enum"
output_dir_sample="./sim_output_sample"

####### Input Parameters #######
offset=0        # offset between fragment start/end values and FASTA sequence
sigma0=3      # Estimate 4or slope of experimental error (slope of sd vs. %D plot)
num_exp_bins=31   # Number of bins
inseq="AAAAAAA"

dat  = { 10   : [21.33928571428864, 15.702380952384537,  6.877976190478665],
         30   : [21.101190476192265, 12.139880952382162, 23.738095238098282],
         60   : [29.78869047619391, 27.208333333334846, 17.913690476193878],
         300  : [32.130952380956586, 20.264880952382892, 31.663690476193555],
         900  : [38.83035714286129, 30.360119047621147, 49.33035714286134],
         3600 : [51.97321428571825, 37.898809523810655, 41.9226190476218], }

'''
dat  = { 10   : [3.903],
         30   : [11.258],
         60   : [30.113],
         300  : [76.145],
         900  : [97.221],
         3600 : [99.9999] }
'''
###############################
###   System Setup:
###############################

# create dataset
dataset = data.Dataset("Sim", data.Conditions(), inseq, offset=0)
pep = dataset.create_peptide(inseq, 1)
for d in dat.keys():
    tp = pep.add_timepoint(d)
    for r in dat[d]:
        tp.add_replicate(r, 1)


# Initialize model  (name, FASTA sequence, offset)
sys = system.System(output_dir=outputdir, noclobber=False)
mol = sys.add_macromolecule(inseq, "Single_pep")

state = mol.get_apo_state()
state.add_dataset(dataset)
om = model.ResidueGridModel(state, grid_size=num_exp_bins)
state.set_output_model(om)
sys.output.initialize_output_model_file(state, om.pf_grids)

output_models = [om]
states = [state]

for i in range(3):
    s = mol.add_state(name="Apo"+str(i+1))
    om = model.ResidueGridModel(state, grid_size=num_exp_bins)
    s.set_output_model(om)
    s.add_dataset(dataset)
    states.append(s)
    output_models.append(om)
    sys.output.initialize_output_model_file(state, om.pf_grids)

#sampler = sampling.EnumerationSampler(sys)
#sampler.run(write=True)

#pof = analysis.ParseOutputFile(outputdir + "/models_scores_sigmas-Apo.dat", state)

sys.output.change_output_directory(output_dir_sample)
#sys.output.initialize_output_model_file(state, output_model.pf_grids)


sampler = sampling.MCSampler(sys)
sampler.run(10000, 2.0, write=True)

pof = analysis.ParseOutputFile(output_dir_sample + "/models_scores_sigmas-Apo.dat", states[0])
pof1 = analysis.ParseOutputFile(output_dir_sample + "/models_scores_sigmas-Apo1.dat", states[1])
pof2 = analysis.ParseOutputFile(output_dir_sample + "/models_scores_sigmas-Apo2.dat", states[2])
pof3 = analysis.ParseOutputFile(output_dir_sample + "/models_scores_sigmas-Apo3.dat", states[3])


plots.plot_residue_protection_factors([pof, pof1, pof2, pof3], num_best_models=1000, sort_sectors=True, show=True)

#plots.plot_po_model_scores(pof)
#plots.plot_po_model_scores(pof2)

#for i in range(2,10):
#    pof2.cluster_models_kmeans(nmodels=1000, nclust=i)

'''
exit()

for pep in dataset.get_peptides():
    for tp in pep.get_timepoints():
        #try:
        i = tp.get_replicates()[0]
        rep_score = -1*math.log(state.scoring_function.replicate_score(tp.get_model_deuteration()/pep.num_observable_amides*100, tp.get_replicates()[0].deut, tp.get_sigma()))
        print pep.sequence, tp.time, tp.get_model_deuteration()/pep.num_observable_amides*100, tp.get_replicates()[0].deut, tp.get_score(), rep_score, "|", tp.get_sigma(), -1*math.log(state.scoring_function.experimental_sigma_prior(tp.get_sigma(), sigma0))
        #except:
        #    pass

sampler.run(20000, 2, write=True)

for pep in dataset.get_peptides():
    for tp in pep.get_timepoints():
        #try:
        i = tp.get_replicates()[0]
        rep_score = -1*math.log(state.scoring_function.replicate_score(tp.get_model_deuteration()/pep.num_observable_amides*100, tp.get_replicates()[0].deut, tp.get_sigma()))
        print pep.sequence, tp.time, tp.get_model_deuteration()/pep.num_observable_amides*100, tp.get_replicates()[0].deut, tp.get_score(), rep_score, "|", tp.get_sigma(), -1*math.log(state.scoring_function.experimental_sigma_prior(tp.get_sigma(), sigma0))
        #except:
        #    pass
'''