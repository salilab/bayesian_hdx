#!/usr/bin/env python

import os
import sys
import argparse
import sampling
import system
import model
import hxio

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Bayesian Analysis of HDX-MS Data.\nThis script process a CVS datafile exported by HDXWorkbench')

    parser.add_argument('-w', help='HDXWorkbench CSV file', required=True)
    parser.add_argument('-o','--outputdir', help='Output directory.', required=True)
    parser.add_argument('--init', help='How to initialize - either "random" or "enumerate". '
                                       'Enumerate is slower but sampling will converge faster. '
                                       'Default: enumerate',
                        default="enumerate",
                        required=False)
    parser.add_argument('--num_exp_bins', help='Number of log(kex) values for sampling. 20 is generally sufficient. '
                                               'Default: 20',
                        default=20,
                        required=False)
    parser.add_argument('--offset', help='Offset between fragment start/end values and FASTA sequence. '
                                         'Default: 0',
                        default=0, required=False)
    parser.add_argument('--sigma0', help='Estimate for experimental error in %%D Units. '
                                         'Default: 5',
                        default=5, required=False)
    parser.add_argument('--annealing_steps', help='Steps per temperature in annealing - 100-200 sufficient. '
                                                  'Default: 20',
                        default=20, required=False)
    parser.add_argument('--nsteps', help='Equilibrium steps. 5000 to 10000. '
                                         'Default: 1000',
                        default=1000, required=False)

    parser.add_argument('--saturation', help='Deuterium saturation in experiment. '
                                         'Default: 1.0',
                        default=1.0, required=False)

    args = parser.parse_args()
    if not os.path.exists(args.w):
        print('HDXWorkbench CSV file should exist')
        sys.exit(-1)
    if not os.path.exists(args.outputdir):
        print('{0} directory was created'.format(args.outputdir))
        os.mkdir(args.outputdir)

    inseq = ""
    percentD = True  # Is the data in percent D (True) or Deuterium units? - Always percentD for Workbench.

    print("Reading sequence from file {}".format(args.w))
    with open(args.w) as fin:
        for line in fin:
            words = line.split(",")
            if words[0] == "Offset":
                args.offset = int(words[1].strip())
                print('Offset defined in file: {}'.format(args.offset))
            elif words[0] == "Experiment Protein Sequence":
                inseq = words[1].strip()
            elif words[0] == "Deuterium solution concentration":
                args.saturation = float(words[1].strip())
                print('Deuterium saturation defined in file: {}'.format(args.saturation))
            elif words[0] == "Experiment name":
                name = words[1].strip().replace(" ", "_")

    if inseq == "":
        print('HDX Workbench file does not contain FASTA sequence. '
              'Please manually add the sequence to the command line using the flag -s or --inseq')
        sys.exit(-1)
    print('Using sequence from file:\n{}'.format(inseq))

    # Initialize model
    sys = system.System(output_dir=args.outputdir, noclobber=False)
    mol = sys.add_macromolecule(inseq, "VDR", initialize_apo=False)

    # Import data
    datasets = hxio.import_HDXWorkbench(args.w,  # Workbench input file
                                        macromolecule=mol,
                                        sequence=None,  # FASTA sequence string
                                        error_estimate=args.sigma0)  # The initial estimate for experimental SD in % deuterium untis.

    sigs = []
    for d in datasets:
        sigs += d.get_all_tp_avg_sigmas()

    # Add data to molecule states and initialize models
    for s in range(len(mol.get_states())):
        state = mol.get_states()[s]
        state.add_dataset(datasets[s])
        output_model = state.set_output_model(model.ResidueGridModel(state, grid_size=args.num_exp_bins))
        sys.output.initialize_output_model_file(state, output_model.pf_grids)

    sys.get_output().write_datasets()

    sampler = sampling.MCSampler(sys, sigma_sample_level="timepoint")

    # First, run a short minimization step
    sampler.run(50, 0.0001, write=True)

    '''
    for dataset in datasets:
      for pep in dataset.get_peptides():
        for tp in pep.get_timepoints():
            #try:
            i = tp.get_replicates()[0]
            rep_score = -1*math.log(state.scoring_function.replicate_score(tp.get_model_deuteration()/pep.num_observable_amides*100, tp.get_replicates()[0].deut, tp.get_sigma()))
            print pep.sequence, tp.time, tp.get_model_deuteration()/pep.num_observable_amides*100, tp.get_replicates()[0].deut, tp.get_score(), rep_score, "|", tp.get_sigma(), -1*math.log(state.scoring_function.experimental_sigma_prior(tp.get_sigma(), sigma0))
            #except:
            #    pass
    '''
    # Slowly cool system
    sampler.run(args.annealing_steps, 3)

    # Optional modification of the number of grid points that the residue can shift (default is 5)
    # (False) allows for any bin in the sampling space to be chosen at random
    # sampler.residue_sampler.set_adjacency(True, 4)

    sampler.run(args.annealing_steps, 2)
    # sampler.residue_sampler.set_adjacency(True, 3)
    sampler.run(args.annealing_steps, 1)

    # This temperature tends to sit around 15% MC acceptance rate, which seems to be good.
    sampler.run(args.nsteps, 1, write=True)
