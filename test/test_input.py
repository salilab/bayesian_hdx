'''
Test the input of various file types
'''
from __future__ import print_function
import utils
import unittest
import os

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

import hxio
import data

input_path = os.path.join(TOPDIR, "test", "input")

class TestHelperFunctions(unittest.TestCase):

    def test_read_fasta(self):
        infile = os.path.join(input_path, "test.fasta")
        seq = ("TEST_PROTEIN1", "ANIMAGINARYPEPTIDE")


        i = hxio.read_fasta(infile)

        (h, s) = next(i)

        self.assertEqual(h, seq[0])
        self.assertEqual(s, seq[1])


        in2file = os.path.join(input_path, "test4.fasta")

        i = hxio.read_fasta(in2file)

        seqs = [("TEST_PROTEIN1", "THECATINTHEHAT"),
            ("TEST_PROTEIN2", "GREENEGGSANDHAM"),
            ("TEST_PROTEIN3", "THEKINGSSTILTS"),
            ("TEST_PROTEIN4", "YERTLETHETERTLE")]


        for s in range(len(seqs)):
            tup = next(i)
            self.assertEqual(tup, seqs[s])



class TestImportFiles(unittest.TestCase):

    def test_import_workbench(self):
        infile = os.path.join(input_path, "Workbench_VDR_VD3_01.csv")
        datasets = hxio.import_HDXWorkbench(infile)
        self.assertEqual(len(datasets), 2)

    def test_import_columns(self):

        infile = os.path.join(input_path, "HXColumns_test_small.csv")
        fastafile = os.path.join(input_path, "test.fasta")

        sequence = next(hxio.read_fasta(fastafile))[1]
        #print("JDHS", sequence.next())
        dataset = hxio.import_HXcolumns(infile, sequence)

        self.assertEqual(len(dataset.get_peptides()),1)
        self.assertEqual(len(dataset.get_peptides()[0].get_timepoints()), 9)
        self.assertEqual(len(dataset.get_peptides()[0].get_timepoints()[0].get_replicates()), 2)

        total_observations = 0
        for pep in dataset.get_peptides():
            for tp in pep.get_timepoints():
                total_observations+=len(tp.get_replicates())

        self.assertEqual(total_observations, 18)



if __name__ == '__main__':
    unittest.main()
