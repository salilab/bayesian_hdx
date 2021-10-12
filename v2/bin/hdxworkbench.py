#!/usr/bin/env python

import os
import sys
import argparse
import sampling
import system
import model
import hxio
import analysis
from operator import itemgetter


def get_best_scoring_models(o, minsize=100):
    new_pof = analysis.deepcopy(o.pof1)
    pof_all = analysis.concatenate_pofs(new_pof, o.pof2)
    smt = pof_all.get_models(return_pf=False)
    smt.sort(key=itemgetter(0))
    sum2=0.0
    sum1=0.0
    finali=minsize #HB store the final i; default is the minimum number of points to use, obviously!
    minstderr2=1.0E+34 #HB 23-May-2018. Can never be negative, so a dumb value to initialize for debugging purposes.
    print('Processing {} models'.format(len(smt)))
    for i in range(0,len(smt),1):
        sum1+=smt[i][0]
        sum2+=smt[i][0]*smt[i][0]
        avg=sum1/(i+1)
        var=(sum2/(i+1)) - (avg * avg)
        stderr2=var/(i+1) #square of stderr, monotonic with stderr and faster to calculate
        if (i >= minsize) and (stderr2 < minstderr2):
            # print('Adjusting stderr: {}'.format(stderr2))
            finali=i+1 #HB this is the i'th entry in smt, which is the current minimum in the stderr
            minstderr2=stderr2 #HB 23-May-2018
    return (finali, minstderr2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Bayesian Analysis of HDX-MS Data.\nThis script process a CVS datafile exported by HDXWorkbench')

    parser.add_argument('-w', help='HDXWorkbench CSV file', required=True)
    parser.add_argument('--mol_name', help='Molecule name',
                        required=True)
    parser.add_argument('-o','--outputdir', help='Output directory.', required=True)
    parser.add_argument('--control', help='Control sample name', required=True)
    parser.add_argument('--ligand', help='Ligand sample name', required=True)
    parser.add_argument('--init', help='How to initialize - either "random" or "enumerate". '
                                       'Enumerate is slower but sampling will converge faster. '
                                       'Default: enumerate',
                        default="enumerate",
                        required=False)
    parser.add_argument('--num_exp_bins', help='Number of log(kex) values for sampling. 20 is generally sufficient. '
                                               'Default: 20',
                        default=20,
                        type=int,
                        required=False)
    parser.add_argument('--offset', help='Offset between fragment start/end values and FASTA sequence. '
                                         'Default: 0',
                        default=0,
                        type=float,
                        required=False)
    parser.add_argument('--sigma0', help='Estimate for experimental error in %%D Units. '
                                         'Default: 5',
                        default=5,
                        type=float,
                        required=False)
    parser.add_argument('--annealing_steps', help='Steps per temperature in annealing - 100-200 sufficient. '
                                                  'Default: 20',
                        default=20,
                        type=int,
                        required=False)
    parser.add_argument('--nsteps', help='Equilibrium steps. 5000 to 10000. '
                                         'Default: 1000',
                        default=1000,
                        type=int,
                        required=False)

    parser.add_argument('--saturation', help='Deuterium saturation in experiment. '
                                         'Default: 1.0',
                        default=1.0,
                        type=float,
                        required=False)

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
    mol = sys.add_macromolecule(inseq, args.mol_name, initialize_apo=False)

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
    sampler.run(100, 0.0001, write=True)
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

    files = [f for dr, ds, files in os.walk(args.outputdir) for f in files if f.endswith('.dat')]
    control_files = []
    ligand_files = []
    for f in files:
        if args.control in f:
            control_files.append(os.path.join(args.outputdir, f))
        elif args.ligand in f:
            ligand_files.append(os.path.join(args.outputdir, f))
    print('Using control files: {}'.format(control_files))
    print('Using ligand files: {}'.format(ligand_files))
    oa = analysis.OutputAnalysis(control_files)
    oa1 = analysis.OutputAnalysis(ligand_files)

    num_models, minstderr = get_best_scoring_models(oa, 100)
    num_models1, minstderr1 = get_best_scoring_models(oa1, 100)
    print('Control best num of models {} with stderr: {}'.format(num_models, minstderr))
    print('Ligand best num of models {} with stderr: {}'.format(num_models1, minstderr1))

    conv = oa.get_convergence(num_models)
    conv1 = oa1.get_convergence(num_models1)

    distmat = conv.get_distance_matrix(num_models=num_models)
    distmat1 = conv1.get_distance_matrix(num_models=num_models1)

    cutoff_list = conv.get_cutoffs_list(1.0)
    cutoff_list1 = conv1.get_cutoffs_list(1.0)

    pvals, cvs, percents = conv.get_clusters(cutoff_list)
    pvals1, cvs1, percents1 = conv1.get_clusters(cutoff_list1)

    sampling_precision,pval_converged,cramersv_converged,percent_converged = conv.get_sampling_precision(cutoff_list, pvals, cvs, percents)
    sampling_precision1,pval_converged1,cramersv_converged1,percent_converged1 = conv1.get_sampling_precision(cutoff_list1, pvals1, cvs1, percents1)

    pofs = conv.cluster_at_threshold_and_return_pofs(sampling_precision)
    pofs1 = conv1.cluster_at_threshold_and_return_pofs(sampling_precision1)

    dhdx = analysis.DeltaHDX(pofs[0], pofs1[0])
    diff, Z, mean1, mean2, sd1, sd2 = dhdx.calculate_dhdx()
    dhdx.write_dhdx_file(prefix='{}/'.format(args.outputdir))


