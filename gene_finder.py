# -*- coding: utf-8 -*-
"""
Gene Finder Project

@author: Alexander Li

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if(nucleotide == 'A'):
        return 'T'
    elif(nucleotide == 'G'):
        return 'C'
    elif(nucleotide == 'C'):
        return 'G'
    elif(nucleotide == 'T'):
        return 'A'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reversedna = ""
    count = 0
    while count < len(dna):
        nuc = get_complement(dna[count])
        reversedna += nuc
        count += 1
    return reversedna[::-1]


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    findStop = 0
    i = 0
    testCodon = ""
    while i <= len(dna):
        testCodon = dna[i:i+3]
        if(testCodon == 'TAG' or testCodon == 'TAA' or testCodon == 'TGA'):
            findStop = 1
            return dna[:i]
            break
        i += 3
    if(findStop == 0):
        return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    ORFs = []
    i = 0
    for x in range(len(dna)):
        if(dna[i:i+3] == "ATG"):
            ORFs.append(rest_of_ORF(dna[i:]))
            i = i + len(str(rest_of_ORF(dna[i:])))
        else:
            i = i+3
    return ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    ORFs = []
    for i in range(3):
        ORFs.extend(find_all_ORFs_oneframe(dna[i:]))
    return ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    ORFs = []
    forwardDna = dna
    backwardsDna = get_reverse_complement(dna)
    ORFs.extend(find_all_ORFs(forwardDna))
    ORFs.extend(find_all_ORFs(backwardsDna))
    return ORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    ORFboth = find_all_ORFs_both_strands(dna)
    maxORF = -1
    for i in range(len(ORFboth)):
        if(len(ORFboth[maxORF]) < len(ORFboth[i])):
            maxORF = i
            i += 1
    return ORFboth[maxORF]


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    randomDNA = dna
    maxlengthDNA = 0
    for i in range(num_trials):
        if(maxlengthDNA < len(longest_ORF(randomDNA))):
            maxlengthDNA = len(longest_ORF(randomDNA))
        shuffle_string(randomDNA)
    return maxlengthDNA


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    amino = []

    if len(dna) % 3 == 1:
        dna = dna[:-1]  # takes the last dna codon out if not a multiple of 3
    if len(dna) % 3 == 2:
        dna = dna[:-2]  # same thing but with 2 at the end
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        amino.append(aa_table[codon])
    return ''.join(amino)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    aasequence = []
    limit = longest_ORF_noncoding(dna, 1000)  # sets reasonable limit for dna
    allOrfs = find_all_ORFs_both_strands(dna)
    for oneOrf in allOrfs:
        if len(oneOrf) > limit:
            aasequence.append(coding_strand_to_AA(oneOrf))
    return aasequence


if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(coding_strand_to_AA, globals())
