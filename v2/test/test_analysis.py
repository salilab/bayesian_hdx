'''
Test objects and functions in analysis.py
'''
from __future__ import print_function
import unittest
import os
import analysis
import numpy

input_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'input'))

print(input_dir)

class TestParseOutputFile(unittest.TestCase):
    '''
    ParseOutputFile reads an output file created from a sampler
    '''

    def test_create_pof(self):
        output_file = input_dir + "/cytC_output/prod1/models_scores_sigmas-CytC_pH_6.5.dat"

        pof=analysis.ParseOutputFile(output_file)
        self.assertEqual(pof.output_file, output_file)
        self.assertEqual(len(pof.get_datasets()),1)
        self.assertEqual(pof.molecule_name,"CytC")
        self.assertEqual(pof.grid_size,50)

        pof.clear_models()
        self.assertEqual(len(pof.models), 0)

    def test_get_models_datasets(self):
        output_file = input_dir + "/cytC_output/prod1/models_scores_sigmas-CytC_pH_6.5.dat"

        pof=analysis.ParseOutputFile(output_file) 
        self.assertEqual(len(pof.get_all_models()), 5100)     
'''
class TestOutputAnalysis(unittest.TestCase):

    def create_output_analysis(self):
        files = [input_dir+"/cytC_output/prod"+str(i)+"/models_scores_sigmas-CytC_pH_6.5.dat" for i in range(1,7)]
        oa = analysis.OutputAnalysis(files)

        return oa

    def test_concatenate_output_files(self):
        oa = self.create_output_analysis()
        self.assertEqual(len(oa.output_files), 6)
        self.assertEqual(len(oa.pof1.models), 5100*3)
        self.assertEqual(len(oa.pof2.models), 5100*3)
        
    def test_get_best_scoring_models(self):
        oa = self.create_output_analysis() 
        self.assertEqual(len(oa.get_all_scores()), 5100*6)        
        self.assertEqual(len(oa.get_best_scoring_models(1000)), 1000)     

class TestConvergence(unittest.TestCase):

    def create_convergence(self):
        files = [input_dir+"/cytC_output/prod"+str(i)+"/models_scores_sigmas-CytC_pH_6.5.dat" for i in range(1,7)]
        oa = analysis.OutputAnalysis(files)
        return oa.get_convergence(100)

    def test_pvalue_cohend(self):
        conv = self.create_convergence()
        #print(conv.total_score_pvalue_and_cohensd(1000))
        self.assertIsInstance(conv.total_score_pvalue_and_cohensd(100), tuple)

    def test_sampling_convergence(self):
        conv = self.create_convergence()

        n = 100
        # First, calculate the distance matrix (num_models)
        distmat = conv.get_distance_matrix(num_models=n)
        self.assertEqual(distmat.shape[0], n*2)

        mind = distmat.min()
        maxd = distmat.max()

        # Second, get the cutoffs list (num_cutoffs)
        cutoff_list = conv.get_cutoffs_list(1)
        self.assertEqual(mind, cutoff_list[0])

        # Third, cluster at each cutoff list and return the list of convergence metrics
        pvals, cvs, percents = conv.get_clusters(cutoff_list)

        self.assertEqual(len(pvals), len(cutoff_list))
        self.assertEqual(pvals[-1], 1.0)
        self.assertEqual(cvs[-1], 0.0)
        self.assertEqual(percents[-1], 100.0)

        sampling_precision,pval_converged,cramersv_converged,percent_converged = conv.get_sampling_precision(cutoff_list, pvals, cvs, percents)

        pofs = conv.cluster_at_threshold_and_return_pofs(sampling_precision)

        num_mods = 0
        for p in pofs:
            num_mods += len(p.models)

        self.assertEqual(num_mods, n*2)
        self.assertGreater(len(pofs[0].models),len(pofs[1].models))
'''

class TestDeltaHDX(unittest.TestCase):

    def test_DHDX_same_file(self):

        output_file = input_dir + "/cytC_output/prod1/models_scores_sigmas-CytC_pH_6.5.dat"
        pof1 = analysis.ParseOutputFile(output_file)
        pof2 = analysis.ParseOutputFile(output_file)

        dhdx = analysis.DeltaHDX(pof1, pof2)

        diff, Z,mean1, mean2, sd1, sd2  = dhdx.calculate_dhdx()

        for res in pof1.observed_residues:
            self.assertEqual(diff[res-1],0)
            self.assertEqual(Z[res-1],0)

    def test_DHDX_different_files(self):

        output_file1 = input_dir + "/cytC_output/prod1/models_scores_sigmas-CytC_pH_6.5.dat"
        output_file2 = input_dir + "/cytC_output/prod1/models_scores_sigmas-CytC_pH_7.4.dat"
        pof1 = analysis.ParseOutputFile(output_file1)
        pof2 = analysis.ParseOutputFile(output_file2)

        dhdx = analysis.DeltaHDX(pof1, pof2)

        diff, Z, mean1, mean2, sd1, sd2 = dhdx.calculate_dhdx()

        print(Z.shape)

        dhdx.write_dhdx_file()

        


if __name__ == '__main__':
    unittest.main()