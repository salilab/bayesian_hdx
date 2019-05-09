'''
Test objects and functions in analysis.py
'''
from __future__ import print_function
import unittest
import os
import analysis

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

    def test_get_models_datasets(self):
        output_file = input_dir + "/cytC_output/prod1/models_scores_sigmas-CytC_pH_6.5.dat"

        pof=analysis.ParseOutputFile(output_file) 
        self.assertEqual(len(pof.get_all_models()), 5100)     

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



if __name__ == '__main__':
    unittest.main()