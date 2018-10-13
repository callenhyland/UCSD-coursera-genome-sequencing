# -*- coding: utf-8 -*-
"""
Functions for Week 3 of UCSD Genome Sequencing Course
Callen Hyland, August 2018
"""


def cyclic_sequence(test_spec, masses):
    """Calculates sequence of peptide from a mass spectrum
    Inputs:
        test_spec : mass spectrum of circular peptide
        masses : table of integer amino acid masses
    Returns: list of candidate peptides
    """
    final_peptides = []
    parent_mass = max(test_spec)
    peptides = [[]]

    while len(peptides) != 0:
        new_peptides = []
        for pep in peptides:
            for mass in masses:

                new_pep = []
                new_pep.extend(pep)
                new_pep.append(mass)

                if sum(new_pep) == parent_mass:
                    new_pep_spec = sorted([sum(i) for i in cyclic_subpeptides(new_pep)])
                    if ((len(new_pep_spec) == len(test_spec)) and (all(i in test_spec for i in new_pep_spec))):
                        final_peptides.append(new_pep)
                        break

                consistent = True

                for sub in linear_subpeptides(new_pep):
                    sum_sub = sum(sub)
                    if sum_sub not in test_spec:
                        consistent = False
                        break

                if consistent:
                    new_peptides.append(new_pep)

        peptides = new_peptides

    return(final_peptides)


def reverse_complement(dna_sequence):
    """Finds reverse complement of a DNA sequence
    Input: DNA sequence, all caps, ATGC only
    Returns: reverse complement of sequence
    """
    bp_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
    comp = ''.join([bp_dict[b] for b in dna_sequence])
    reverse_comp = comp[::-1]
    return(reverse_comp)


def transcribe(dna_sequence):
    """Transcribes DNA sequence into RNA
    Input: DNA sequence, all caps, ATGC only
    Returns: RNA sequence, AUGC
    """
    bp_dict = {"T":"U", "A":"A", "G":"G", "C":"C"}
    return(''.join([bp_dict[b] for b in dna_sequence]))
    
    
def translate(rna_sequence, codons, a = 0):
    """Translates RNA sequence into an amino acid sequence
    Inputs:
        rna_sequence : RNA sequence, AUGC all caps
        codons: dictionary with key = codon, value = amino acid
        a : translation frame, starting position
    Returns: amino acid sequence
    """
    translated = ""
    while a < len(rna_sequence)-2:
        triplet = rna_sequence[a:a+3]
        aa = codons[triplet]
        translated = translated + aa
        a = a+3
    return(translated)


def peptide_encoding(dna, peptide, codons):
    """Finds sub-sequences encoding a peptide within both forward and reverse
    DNA strands. Calls find_peptide function.
    Inputs:
        dna : DNA sequence, AUGC all caps
        peptide : amino acid sequence of peptide
        codons: dictionary with key = codon, value = amino acid
    Returns: list of substrings encoding peptide
    """
    rc_dna = reverse_complement(dna)
    sub_dna = find_peptide(dna, peptide, codons)
    rev_subseqs = find_peptide(rc_dna, peptide, codons)
    for i in rev_subseqs:
        sub_dna.append(reverse_complement(i))
    return(sub_dna)


def find_peptide(dna, peptide, codons):
    """Finds sub-sequences encoding a peptide within a DNA sequence
    Inputs:
        dna : DNA sequence, AUGC all caps
        peptide : amino acid sequence of peptide
        codons: dictionary with key = codon, value = amino acid
    Returns: list of substrings encoding peptide
    """
    substrings = []
    pep_len = len(peptide)
    rna = transcribe(dna)
    
    for i in range(3):
        protein = translate(rna, codons, a=i)
        #print(protein)
        ind = 0
        while ind >= 0:
            ind = protein.find(peptide, ind)
            #print(ind)
            if ind == -1:
                break
            substrings.append(dna[(ind*3+i):(ind*3+pep_len*3+i)])
            ind += 1
    return(substrings)


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


# generate list of all sub-peptides from a linear sequence

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



def calc_mass(subpep_seq, mass_dict):
    """calculates mass of a peptide
    Inputs:
        subpep_seq : amino acid sequence of peptide
        mass_dict : dictionary of amino acid masses
    Returns: mass of peptide
    """
    mass = 0
    for i in subpep_seq:
        mass = mass + mass_dict[i]
    return(mass)


def calc_spectrum(pep_seq, mass_dict):
    """Calculates the mass spectrum of a circular peptide
    Inputs:
        pep_seq : amino acid sequence of  circular peptide
        mass_dict : dictionary of amino acid masses
    Returns: mass spectrum of circular peptide
    """
    subpep_list = subpeptides(pep_seq)
    mass_list = [calc_mass(pep, mass_dict) for pep in subpep_list]
    mass_list.append(0)
    mass_list.sort()
    return(mass_list)
