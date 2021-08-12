'''
Test objects and functions in analysis.py
'''
from __future__ import print_function
import unittest
import os

input_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'input'))

print(input_dir)


class TestParseOutputFile(unittest.TestCase):
    '''
    ParseOutputFile reads an output file created from a sampler
    '''

    output_file = input_dir

    def test_create_pof(self):
        pass


class TestOutputAnalysis(unittest.TestCase):

    def test_concatenate_output_files(self):
        pass
