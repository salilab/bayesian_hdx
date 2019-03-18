import sys
sys.path.append( "../../pyext/src" )
import plots
import analysis

outputdir="./test_simulated_data_50k/"

pof = analysis.ParseOutputFile(outputdir + "/models_scores_sigmas-Apo.dat", "Apo")
pof2 = analysis.ParseOutputFile(outputdir + "/models_scores_sigmas-Apo2.dat", "Apo2")
pof.generate_datasets()
pof2.generate_datasets()
#pof.calculate_random_sample_convergence()
#pof2.calculate_random_sample_convergence()

conv = analysis.Convergence(pof, pof2, 500)

print(conv.total_score_pvalue_and_cohensd())

ranges = [0.01, 0.1, 0.2, 0.3, 0.4]
#print(conv.get_clusters(ranges))

#exit()

#print(conv.residue_pvalue_and_cohensd())
plots.plot_incorporation_curve_fits(pof, 500, outputdir+"/incorporation_plots/")
plots.plot_incorporation_curve_fits(pof2, 500, outputdir+"/incorporation_plots2/")
plots.plot_po_model_scores(pof, False, outputdir+"/apo_total_score.png", 500)
plots.plot_po_model_scores(pof2, False, outputdir+"/apo2_total_score.png", 500)
plots.plot_residue_protection_factors([pof, pof2], num_best_models=500, sort_sectors=True)
