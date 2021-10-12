"""
   Useful tools for HDX and MS data analysis

   Can be used as a standalone library outside of IMP
   and the RResolvedHX
"""

from __future__ import print_function
import numpy
import numpy.random
import math


AA_list = "ARNDCQEGHILKMFPSTWYV"


def get_residue_neighbor_effects(AA, pDcorr, T):
    # For each residue, a tuple containing:
    # 0:Milne acid lambda
    # 1:Milne acid rho
    # 2:Milne base lambda
    # 3:Milne base rho

    R = 1.987

    # Calculate Temp dependent pKa's
    pK_D = -1 * numpy.log10(
        10**(-1*4.48)*numpy.exp(-1*1000*((1.0/T-1.0/278)/R)))
    pK_E = -1 * numpy.log10(
        10**(-1*4.93)*numpy.exp(-1*1083*((1.0/T-1.0/278)/R)))
    pK_H = -1 * numpy.log10(
        10**(-1*7.42)*numpy.exp(-1*7500*((1.0/T-1.0/278)/R)))

    eff_dict = {"A": (0.0, 0.0, 0.0, 0.0),
                "R": (-0.59, -0.32, 0.07671225, 0.22),
                "N": (-0.58, -0.13, 0.49, 0.32),
                "C": (-0.54, -0.46, 0.62, 0.55),
                "Q": (-0.47, -0.27, 0.06, 0.20),
                "G": (-0.22, 0.21817605, 0.26725157, 0.17),
                "I": (-0.91, -0.59, -0.73, -0.23),
                "L": (-0.57, -0.13, -0.57625273, -0.21),
                "K": (-0.56, -0.29, -0.04, 0.12),
                "M": (-0.64, -0.28, -0.00895484, 0.11),
                "F": (-0.52, -0.43, -0.23585946, 0.06313159),
                "P": ("", -0.19477347, "", -0.24),
                "S": (-0.43799228, -0.38851893, 0.37, 0.29955029),
                "T": (-0.79, -0.46807313, -0.06625798, 0.20),
                "W": (-0.40, -0.44, -0.41, -0.11),
                "Y": (-0.41, -0.37, -0.27, 0.05),
                "V": (-0.73902227, -0.30, -0.70193448, -0.14),
                "NT": ("", -1.32, "", 1.62)}

    # Ionizable AA data from
    # Susumu Mori, Peter C.M. van Zijl, and David Shortle
    # PROTEINS: Structure, Function, and Genetics 28:325-332 (1997)

    if AA == "D":
        ne0 = numpy.log10(
            10**(-0.9-pDcorr)/(10**(-pK_D)+10**(-pDcorr))
            + 10**(0.9-pK_D)/(10**(-pK_D)+10**(-pDcorr)))
        ne1 = numpy.log10(
            10**(-0.12-pDcorr)/(10**(-pK_D)+10**(-pDcorr))
            + 10**(0.58-pK_D)/(10**(-pK_D)+10**(-pDcorr)))
        ne2 = numpy.log10(
            10**(0.69-pDcorr)/(10**(-pK_D)+10**(-pDcorr))
            + 10**(0.1-pK_D)/(10**(-pK_D)+10**(-pDcorr)))
        ne3 = numpy.log10(
            10**(0.6-pDcorr)/(10**(-pK_D)+10**(-pDcorr))
            + 10**(-0.18-pK_D)/(10**(-pK_D)+10**(-pDcorr)))
    elif AA == "E":
        ne0 = numpy.log10(
            10**(-0.6-pDcorr)/(10**(-pK_E)+10**(-pDcorr))
            + 10**(-0.9-pK_E)/(10**(-pK_E)+10**(-pDcorr)))
        ne1 = numpy.log10(
            10**(-0.27-pDcorr)/(10**(-pK_E)+10**(-pDcorr))
            + 10**(0.31-pK_E)/(10**(-pK_E)+10**(-pDcorr)))
        ne2 = numpy.log10(
            10**(0.24-pDcorr)/(10**(-pK_E)+10**(-pDcorr))
            + 10**(-0.11-pK_E)/(10**(-pK_E)+10**(-pDcorr)))
        ne3 = numpy.log10(
            10**(0.39-pDcorr)/(10**(-pK_E)+10**(-pDcorr))
            + 10**(-0.15-pK_E)/(10**(-pK_E)+10**(-pDcorr)))
    elif AA == "H":
        ne0 = numpy.log10(
            10**(-0.8-pDcorr)/(10**(-pK_H)+10**(-pDcorr))
            + 10**(0-pK_H)/(10**(-pK_H)+10**(-pDcorr)))
        ne1 = numpy.log10(
            10**(-0.51-pDcorr)/(10**(-pK_H)+10**(-pDcorr))
            + 10**(0-pK_H)/(10**(-pK_H)+10**(-pDcorr)))
        ne2 = numpy.log10(
            10**(0.8-pDcorr)/(10**(-pK_H)+10**(-pDcorr))
            + 10**(-0.1-pK_H)/(10**(-pK_H)+10**(-pDcorr)))
        ne3 = numpy.log10(
            10**(0.83-pDcorr)/(10**(-pK_H)+10**(-pDcorr))
            + 10**(0.14-pK_H)/(10**(-pK_H)+10**(-pDcorr)))
    elif AA == "CT":
        ne0 = numpy.log10(
            10**(0.05-pDcorr)/(10**(-pK_E)+10**(-pDcorr))
            + 10**(0.96-pK_E)/(10**(-pK_E)+10**(-pDcorr)))
        ne1 = ""
        ne2 = -1.8
        ne3 = ""
    else:
        (ne0, ne1, ne2, ne3) = eff_dict.get(AA)

    # print(AA, pDcorr, (ne0, ne1, ne2, ne3))

    return (ne0, ne1, ne2, ne3)


ResidueChemicalContent = {
    # Tuple containing number of atoms of
    #  (Carbon, Hydrogen, Nitrogen, Oxygen, Sulfur)
    #  for a free amino aicd
    "A": (3, 5, 1, 1, 0),
    "R": (6, 12, 4, 1, 0),
    "N": (4, 6, 2, 2, 0),
    "D": (4, 5, 1, 3, 0),
    "C": (3, 5, 1, 1, 1),
    "Q": (5, 8, 2, 2, 0),
    "E": (5, 7, 1, 3.0),
    "G": (2, 3, 1, 1, 0),
    "H": (6, 7, 3, 1, 0),
    "I": (6, 11, 1, 1, 0),
    "L": (6, 11, 1, 1, 0),
    "K": (6, 12, 2, 1, 0),
    "M": (5, 9, 1, 1, 1),
    "F": (9, 9, 1, 1, 0),
    "P": (5, 7, 1, 1, 0),
    "S": (3, 5, 1, 2, 0),
    "T": (4, 7, 1, 2, 0),
    "W": (11, 10, 2, 1, 0),
    "Y": (9, 9, 1, 2, 0),
    "V": (5, 9, 1, 1, 0),
    "CT": (0, 1, 0, 1, 0),
    "NT": (0, 1, 0, 0, 0)
    # Eventually add AA modifications
}

ElementMasses = {
    # List of tuples containing the mass and abundance of atoms
    # in biological peptides (C, H, N, O, S)
    #
    # Data from
    # https://www.ncsu.edu/chemistry/msf/pdf/IsotopicMass_NaturalAbundance.pdf
    #
    # Original references:
    # The isotopic mass data is from:
    #   G. Audi, A. H. Wapstra Nucl. Phys A. 1993, 565, 1-65
    #   G. Audi, A. H. Wapstra Nucl. Phys A. 1995, 595, 409-480.
    # The percent natural abundance data is from the 1997 report of the
    # IUPAC Subcommittee for Isotopic Abundance Measurements by
    # K.J.R. Rosman, P.D.P. Taylor Pure Appl. Chem. 1999, 71, 1593-1607.

    "C": [(12.000000, 98.93), (13.003355, 1.07)],
    "H": [(1.007825, 99.9885), (2.14101, 0.0115)],
    "N": [(14.0030764, 99.632), (15.000109, 0.368)],
    "O": [(15.994915, 99.757), (16.999132, 0.038), (17.999160, 0.205)],
    "S": [(31.972071, 94.93), (32.97158, 0.76),
          (33.967867, 4.29), (35.967081, 0.02)]
}


def calc_intrinsic_rate(Laa, Raa, pH, T, La2="A", Ra2="A", log=False):
    ''' Calculates random coil hydrogen exchange rate for amide cooresponding
    to side chain Raa
    @param Laa - amino acid letter code N-terminal to amide
    @param Raa - amino acid letter of amide
    @param pH - The pH of the experiment
    @param T - Temperature of the experiment

    Equation and constants derived from Bai, Englander (1980)
    '''

    if Raa == "P" or Raa == "CT" or Raa == "NT" or Laa == "NT":
        return 0
    # Constants
    pKD20c = 15.05
    ka = 0.694782306
    kb = 187003075.7
    kw = 0.000527046
    R = 1.987
    EaA = 14000
    EaB = 17000
    EaW = 19000
    # the pD is different than the pH by +0.4 units
    pDcorr = pH+0.4

    inv_dTR = (1./T-1./293)/R

    FTa = numpy.exp(-1*EaA*inv_dTR)
    FTb = numpy.exp(-1*EaB*inv_dTR)
    FTw = numpy.exp(-1*EaW*inv_dTR)
    Dplus = 10**(-1*pDcorr)
    ODminus = 10**(pDcorr-pKD20c)

    # Get residue-specific effect factors

    L_ne = get_residue_neighbor_effects(Laa, pDcorr, T)
    R_ne = get_residue_neighbor_effects(Raa, pDcorr, T)

    Fa = L_ne[1]+R_ne[0]
    Fb = L_ne[3]+R_ne[2]

    if La2 == "NT":
        Fa += get_residue_neighbor_effects(La2, pDcorr, T)[1]
        Fb += get_residue_neighbor_effects(La2, pDcorr, T)[3]
    if Ra2 == "CT":
        Fa += get_residue_neighbor_effects(Ra2, pDcorr, T)[0]
        Fb += get_residue_neighbor_effects(Ra2, pDcorr, T)[2]

    Fa = 10**(Fa)
    Fb = 10**(Fb)

    krc = Fa*Dplus*ka*FTa + Fb*ODminus*kb*FTb + Fb*kw*FTw

    return krc


def calculate_number_of_observable_amides(inseq, n_fastamides=2):
    # number of observable amides is equal to peptide length - 2,
    # minus remaining prolines
    num_prolines = inseq.count('P', n_fastamides) \
        + inseq.count('p', n_fastamides)
    return len(inseq) - num_prolines - n_fastamides


def get_sequence_intrinsic_rates(seq, pH, T, log=False):
    i_rates = numpy.zeros(len(seq))
    i_rates[0] = calc_intrinsic_rate("NT", seq[0], pH, T)
    i_rates[1] = calc_intrinsic_rate(seq[0], seq[1], pH, T, La2="NT")
    for n in range(2, len(seq)-1):
        # print(n, seq[n],seq[n+1])
        L = seq[n-1]
        R = seq[n]
        i_rates[n] = calc_intrinsic_rate(L, R, pH, T)

    i_rates[-1] = calc_intrinsic_rate(seq[-2], seq[-1], pH, T, Ra2="CT")
    if log:
        # Suppress divide by zero error.
        with numpy.errstate(divide='ignore'):
            i_rates = numpy.log10(i_rates)

        # print("LOG", seq, i_rates)
        return i_rates
    else:
        return i_rates


def get_pdb_sequence(pdb_file, chain=None, offset=0):
    # Returns macromolecular sequence as a Pandas sequence
    seq = pd.Series()
    with open(pdb_file) as f:
        for rawl in f:
            # print rawl.split()
            line = rawl.split()
            if len(line) < 4:
                continue
            elif line[2] == "CA" and (chain is None or chain == line[4]):
                resnum = int(line[5])+offset
                # print resnum, line, len(self.get_fasta_sequence())
                if len(seq) == 0:
                    seq = pd.Series(ThreeToOne.get(line[3], '.'),
                                    index=[resnum])
                else:
                    seq = seq.append(pd.Series(ThreeToOne.get(line[3], '.'),
                                               index=[resnum]))
    # seq.index = seq.index + offset
    return seq


def calculate_peptide_mass(string):
    ''' Calculates the mass of a peptide string using
    only the most abundant isotopes.
    @param string - character string of one-letter amino acids

    Returns a float
    '''
    mass = 0
    for i in string:
        if i not in AA_list:
            raise Exception(
                "tools.calculate_peptide_mass : unknown amino acid " + i + "")

        els = ResidueChemicalContent(i)

        ellist = ["C", "H", "O", "N", "S"]

        for i in range(len(els)):
            mass += els[i] * ElementMasses(ellist[i])[0][0]

    return mass


def calculate_simple_deuterium_incorporation(rate, time):
    # Calculates the deuterium incorporation for a log(kex)
    # at time = t (seconds) assuming full saturation
    return 1 - math.exp(-1*10**rate*time)


def get_residue_deuteration_at_each_timepoint(dataset, protection_factors):
    # Returns the deuterium incorporation for each residue at each timepoint in
    # the given dataset and given protection factors
    timepoints = dataset.get_set_of_timepoints()
    deuterations_by_time = {}

    for tp in timepoints:
        deuterations = []
        for i in range(len(protection_factors)):
            if (math.isnan(dataset.intrinsic[i]) or i == numpy.inf
                    or i == -1 * numpy.inf):
                deuterations.append(0)
            else:

                log_kex = dataset.intrinsic[i] - protection_factors[i]
                # Deuterium incorporation is scaled by the amount of
                # deuterium in solution
                deut = calculate_simple_deuterium_incorporation(
                    log_kex, tp.time) * dataset.conditions.saturation
                deuterations.append(deut)
            # print(log_kex, t.time, deut)
        tp.set_deuteration(deut)
        deuterations_by_time[tp.time] = deuterations
        # print(deuterations_by_time)

    return deuterations_by_time


def calculate_incorporation(intrinsic, protection_factors, timepoints):
    if len(intrinsic) != len(protection_factors):
        raise Exception(
            "Intrinsic exchange factors and protection factor list not "
            "the same length")

    incorporations = {}
    for n in range(len(intrinsic)):
        res_incorps = {}
        for tp in timepoints:
            log_kex = intrinsic[n] - protection_factors[n]
            res_incorps[tp] = calculate_simple_deuterium_incorporation(
                log_kex, tp)

        incorporations[n+1] = res_incorps

    return incorporations


def get_peptide_deuteration(peptide, protection_factors):
    deuts = {}
    for tp in peptide.get_timepoints():
        deuts[tp.time] = get_timepoint_deuteration(
            peptide, tp.time, protection_factors)
    return deuts


def get_timepoint_deuteration(peptide, time, protection_factors):
    for r in peptide.get_observable_residue_numbers():
        pf = protection_factors[r-1]
        kr = peptide.dataset.intrinsic[r-1]
        deut += calculate_simple_deuterium_incorporation(kr - pf, time) \
            * peptide.dataset.conditions.saturation
    return deut


def get_residue_peptide_deuteration_at_each_timepoint(sequence, peptides,
                                                      protection_factors):
    # Returns the deuterium incorporation for each residue at each timepoint in
    # the given dataset and given protection factors
    # deuterations_by_time = {}

    for pep in peptides:
        dataset = pep.dataset
        # print(pep.sequence, pep.get_observable_residue_numbers())
        for tp in pep.get_timepoints():
            total_deut = 0
            for i in pep.get_observable_residue_numbers():

                if dataset.intrinsic[i-1] * 0 == 0:
                    log_kex = dataset.intrinsic[i-1] - protection_factors[i-1]
                    # Deuterium incorporation is scaled by the amount
                    # of deuterium in solution
                    deut = calculate_simple_deuterium_incorporation(
                        log_kex, tp.time) * dataset.conditions.saturation
                    total_deut += deut
            # print(pep.sequence, tp.time, total_deut)
            tp.set_deuteration(total_deut)
            # deuterations_by_time[tp.time] = deuterations
            # print(deuterations_by_time)

        # return deuterations_by_time


def calc_peptide_isotopic_distribution(string, threshold=0.1):
    ''' Calculates the natrual isotopic distribution of a peptide
    using the multinomial expansion
    @param string - character string of one-letter amino acids
    @param threshold - returned isotopes must be at least this percent of total

    Returns a list of tuples (mass, percent)

    Reference: Yergey (1983)
    '''
    for i in string:
        if i not in AA_list:
            raise Exception(
                "tools.calculate_peptide_mass : unknown amino acid " + i + "")


def calculate_deut(rate, time):
    return 1-math.exp(-10**rate*time)


def simulate_peptide_data(seq, start_res, exch_rates, timepoints,
                          replicates=3, obs_error=5, percentD=True):
    '''
    For a list of exchange rates, calculate simulated data including
    gaussian noise at each timepoint.
    @param exch_rates - list of exchange rates for a set of contiguous residues
    @param timepoints - list of integer timepoints (in seconds)
    @param error - Standard deviation (in total D units)

    Returns a list of lists of 2D incorporation values.
    '''

    peptide = system.Peptide(seq, start_res, start_res+len(seq))

    # Add data to the fragment
    for t in timepoints:
        tp = frag.add_timepoint(t)
        deut = 0
        for r in range(replicates):
            for i in exch_rates[2:]:
                # calc exch
                d = calculate_deut(i, tp.time)
                deut += d + normal(0, d*obs_error/100)

        tp.add_replicate(deut)

    return peptide


def simulate_state_data(state, rates, replicates, timepoints, pH=7, temp=283):
    return 0


def calculate_mass(string, percent_cutoff=100):
    ''' Calculates the average mass of a peptide string
    @param string - AA sequence.  N and C-terms are implicitly added
    @param percent_cutoff - below this percentage, don't report the
           isotope (100 = only the major isotope)
    '''
    return 0


def subsequence_consistency(sequence, subsequence, start_residue):
    # Simple utility that compares a subsequence to a whole sequence
    # Returns True if 100% consistent. False otherwise
    return sequence[
        start_residue-1:start_residue-1+len(subsequence)] == subsequence
