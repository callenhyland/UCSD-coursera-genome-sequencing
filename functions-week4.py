# -*- coding: utf-8 -*-
"""
Functions for Week 4 of UCSD Genome Sequencing Course
Callen Hyland, August 2018
"""

from itertools import combinations

def cyclic_subpeptides(pep_seq):
    """Generates list of all sub-peptides from a circular sequence
    Input: peptide amino acid sequence
    Returns: list of subsequences
    """
    n = len(pep_seq)
    cyclo_seq = pep_seq + pep_seq
    sub_peptides = []

    for size in range(1,n):
        for start in range(n):
            sub_peptides.append(cyclo_seq[start:start+size])

    sub_peptides.append(pep_seq)
    return(sub_peptides)
    

def linear_subpeptides(pep_seq):
    """Generates list of all sub-peptides from a linear sequence
    Input: peptide amino acid sequence
    Returns: list of subsequences
    """
    n = len(pep_seq)
    sub_peptides = []
    for size in range(1,n):
        for start in range(n-size+1):
            sub_peptides.append(pep_seq[start:start+size])

    sub_peptides.append(pep_seq)
    return(sub_peptides)


def spectrum(peptides):
    """generate list of peptide weights from list of peptides
    """
    spec = [0]
    for pep in peptides:
        spec.append(sum(pep))
    return(sorted(spec))


def score_peptides(experimental_spectrum, theoretical_spectrum):
    """Computes score for experimenal spectrum based on overlap with theoretical spectrum
    Inputs:
        experimental_spectrum : list of masses in experimental spectrum
        theoretical_spectrum : list of masses in theoretical spectrum
    Returns: integer score
    """
    score = 0
    for i in experimental_spectrum:
        if i in theoretical_spectrum:
            score += 1
            theoretical_spectrum.remove(i)
    return(score)


def aa_masses(filename):
    """Read amino acid masses from file and return as dictionary
    """
    mass_table = {}
    with open (filename, "r") as myfile:
        for line in myfile:
            aa_mass = line.rstrip().split(" ")
            mass_table[aa_mass[0]] = int(aa_mass[1])
    return(mass_table)



def aa2masses(peptide, masses):
    """Computes list of masses from peptide string and amino acid mass dictionary
    """
    return([masses[i] for i in peptide])


def trim(peptide_leaderboard, exp_spectrum, trim_length):
    """Trims the list of peptides to only the top scoring peptides
    Inputs:
        peptide_leaderboard : list of peptides
        exp_spectrum : experimental mass spectrum
        trim_length : number of top-scoring peptides to keep
    Returns: trimmed list of peptides
    """
    scores = []
    for lb in peptide_leaderboard:
        score = score_peptides(exp_spectrum, spectrum(linear_subpeptides(lb)))
        scores.append(score)
    if len(scores) > trim_length:
        cutoff = sorted(scores, reverse = True)[trim_length-1]
        num_keep = len([i for i in scores if i >= cutoff])
        lb_sorted = [x for s,x in sorted(zip(scores, peptide_leaderboard), reverse = True)]
        return(lb_sorted[0:num_keep])
    else:
        return(peptide_leaderboard)


def cyclic_leaderboard_sequence(exp_spec, trim_length, mass_set):
    """Finds sequence of cyclic peptide from mass spectrum
    Inputs:
        exp_spec : experimental mass spectrum
        trim_length : number of top sequences to retain at each step
        mass_set : set of amino acid masses to use for sequencing
    Returns: Peptide sequence
    """
    parent_mass = max(exp_spec)
    peptides = [[]]
    leader_peptide = []

    while len(peptides) != 0:
        new_peptides = []
        
        # expand peptide list
        for pep in peptides:
            for mass in mass_set:
                new_pep = []
                new_pep.extend(pep)
                new_pep.append(mass)
                new_peptides.append(new_pep)
        
        # test each peptide to see if it's the leader peptide or greater than parent mass
        filtered_peptides = []
        for pep in new_peptides:

            if sum(pep) == parent_mass:
                new_spec = spectrum(cyclic_subpeptides(pep))
                leader_spec = spectrum(cyclic_subpeptides(leader_peptide))
                if score_peptides(exp_spec, new_spec) > score_peptides(exp_spec, leader_spec):
                    leader_peptide = pep
            if sum(pep) <= parent_mass:
                filtered_peptides.append(pep)

        # trim the list to just the top N peptides (and ties)
        peptides = trim(filtered_peptides, exp_spec, trim_length)
    
    return(leader_peptide)


def convolve(spectrum):
    """For every pair of masses in the spectrum, compute the difference
    Input: list of masses in spectrum
    Returns: list of differences
    """
    pairs = list(combinations(spectrum, 2))
    conv = [abs(pair[1]-pair[0]) for pair in pairs]
    return(conv)


def choose_alphabet(spectrum, M):
    """From a spectrum, generates a list of M amino acid masses appearing most frequently
    """
    conv = [i for i in convolve(spectrum) if (i > 56) and (i < 201)]
    counts = [conv.count(i) for i in set(conv)]
    conv_sorted = [x for s,x in sorted(zip(counts, set(conv)), reverse = True)]
    if len(conv_sorted) > M:
        return(conv_sorted[0:M])
    else:
        return(conv_sorted)

