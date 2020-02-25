'''
Test the input of various file types
'''
from __future__ import print_function
import hxio
import data
import unittest
import os

input_path = os.path.dirname(os.path.realpath(__file__))+"/input/"

class TestHelperFunctions(unittest.TestCase):

    def test_read_fasta(self):
        infile = input_path + "test.fasta"
        seq = ("TEST_PROTEIN1", "ANIMAGINARYPEPTIDE")


        i = hxio.read_fasta(infile)

        (h, s) = i.next()

        self.assertEqual(h, seq[0])
        self.assertEqual(s, seq[1])


        in2file = input_path + "test4.fasta"

        i = hxio.read_fasta(in2file)

        seqs = [("TEST_PROTEIN1", "THECATINTHEHAT"),
            ("TEST_PROTEIN2", "GREENEGGSANDHAM"),
            ("TEST_PROTEIN3", "THEKINGSSTILTS"),
            ("TEST_PROTEIN4", "YERTLETHETERTLE")]


        for s in range(len(seqs)):
            tup = i.next()
            self.assertEqual(tup, seqs[s])



class TestImportFiles(unittest.TestCase):

    def test_import_workbench(self):
        infile = input_path + "Workbench_VDR_VD3_01.csv"
        datasets = hxio.import_HDXWorkbench(infile)
        self.assertEqual(len(datasets), 2)

    def test_import_columns(self):

        infile = input_path + "HXColumns_test_small.csv"
        fastafile = input_path + "test.fasta"

        sequence = hxio.read_fasta(fastafile).next()[1]
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
