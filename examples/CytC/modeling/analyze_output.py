# Usage:  /package_home/setup.py python analyze_output.py
import scoring
import system
import analysis
import plots


##########################################
###    File/Directory Setup 
outputdir = "./ph_6.5" # output directory for the simulation results. 
outputdir2 = "./ph_7.4" # output directory for the simulation results.
resrange = (40,105) # Residue range is a tuple in pdb numbering (starts at 1).
num_best_models = 1000


res_pfs = {}
infile = "../data/nmr_milne_1998"
with open(infile) as f:
    for line in f:
        if line[0]=='#' or line[0]=='':
            continue
        else:
            res_pfs[int(line.split()[0])] = float(line.split()[1])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###   Analysis.  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Joining multiple runs and evaluating sampling precision (convergence)

#files = ["../../../test/input/cytC_output/prod1/models_scores_sigmas-CytC_pH_6.5.dat",
#        "../../../test/input/cytC_output/prod2/models_scores_sigmas-CytC_pH_6.5.dat"]
import glob
files6 = glob.glob("../../../test/input/cytC_output/production*/models_scores_sigmas-CytC_pH_6.5.dat")
files7 = glob.glob("../../../test/input/cytC_output/production*/models_scores_sigmas-CytC_pH_7.4.dat")

oa7 = analysis.OutputAnalysis(files7)

conv7 = oa7.get_convergence(100)

# First, calculate the distance matrix (num_models)
distmat7 = conv7.get_distance_matrix(num_models=20)

#mind = distmat.min()
#maxd = distmat.max()

# Second, get the cutoffs list (grid step size)
# The list of cutoffs will be range(min, max, n)
# Make this smaller if max pairwise distance (printed to screen or, distmat.max()) is very small, like, under 5
cutoff_list = conv7.get_cutoffs_list(1.0)

# Third, cluster at each cutoff list and return the list of convergence metrics
pvals, cvs, percents = conv7.get_clusters(cutoff_list)

sampling_precision,pval_converged,cramersv_converged,percent_converged = conv.get_sampling_precision(cutoff_list, pvals, cvs, percents)


# here, sampling_precision is a float in our RMSD-like distance metric.
# It can be equal to the sampling precision, or slightly lower if the convergence metrics
# are very good.

# The output here is a list of clusters, as POF objects, containing the models from that cluster.
# We use these for downstream analysis/plotting 
pofs = conv7.cluster_at_threshold_and_return_pofs(sampling_precision)


# This command plots the protection factor distribution curves. It will do it for up to 5? POF files. After that, you
# run out of colors and it gets too busy anyways.
# - first input is a list of POF objects.
# - num_best_models is self-explanatory for POF. Set to "all" for all of them (should do that if you've already clustered)
# - resrange is the range to plot
# - true_vals are a list of numerical (log) protection values per residue. They will show up as horizontal green lines. Default is None.
# - sort_sectors sorts residues in the sectors by increasing Pf value by model. 
# - outputdir is self explanatory. 
plots.plot_residue_protection_factors(pofs, num_best_models=num_best_models, 
    resrange=resrange, true_vals = res_pfs, sort_sectors=True, outputdir="./")

# Plot the incorporation curves.
#
# Use imagemagick montage to put them all together
# Name

exit()

plots.plot_incorporation_curve_fits(pofs[0], num_best_models, write_plots=True, output_directory=outputdir)

plots.plot_incorporation_curve_fits(pofs[1], num_best_models, write_plots=True, output_directory=outputdir+"cl1")

############################
# Difference between two states
############################
# Get pH_7.4 POFS
files = ["../../../test/input/cytC_output/prod1/models_scores_sigmas-CytC_pH_7.4.dat",
        "../../../test/input/cytC_output/prod2/models_scores_sigmas-CytC_pH_7.4.dat"]

oa7 = analysis.OutputAnalysis(files)
conv7 = oa7.get_convergence(20)
distmat = conv7.get_distance_matrix(num_models=20)
cutoff_list7 = conv7.get_cutoffs_list(1.0)
pvals, cvs, percents = conv7.get_clusters(cutoff_list)
sampling_precision7,pval_converged,cramersv_converged,percent_converged = conv7.get_sampling_precision(cutoff_list, pvals, cvs, percents)
pofs7 = conv7.cluster_at_threshold_and_return_pofs(sampling_precision)

#
# POFs for pH 6.5 are in list pofs
# POFs for pH 7.4 are in list pofs7

dhdx = analysis.DeltaHDX(pofs[0], pofs7[0])
diff, Z, mean1, mean2, sd1, sd2 = dhdx.calculate_dhdx()
dhdx.write_dhdx_file()



