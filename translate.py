#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """
    rna_sequence = rna_sequence.upper()
    rna_list = list(rna_sequence)
    x = 0
    codons = ""
    codon_list = list(codons)
    while x + 2 < len(rna_list):
        triplet = rna_list[x] + rna_list[x+1] + rna_list[x+2]
        codon_list.append(triplet)
        x = x + 3
    new = ""
    for elem in codon_list:
        if genetic_code[elem] == "*":
            break
        else:
            new = new + genetic_code[elem]
    return new

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """
    aa_list = []
    reading_frame_1_string = rna_sequence
    reading_frame_2_string = rna_sequence[1:]
    reading_frame_3_string = rna_sequence[2:]
    reading_frames_tuple = (reading_frame_1_string, reading_frame_2_string, reading_frame_3_string)
    for x in range(3):
        rf = reading_frames_tuple[x]
        rf_upper = rf.upper()
        codon_list = []
        y = 0
        while y + 2 < len(rf_upper):
            triplet = rf_upper[y] + rf_upper[y+1] + rf_upper[y+2]
            codon_list.append(triplet)
            y = y + 3
        for number in range(len(codon_list)):
            if genetic_code[codon_list[number]] == "M":
                z = number
                aa_seq = ""
                while z < len(codon_list) and genetic_code[codon_list[z]] != "*":
                    aa_seq = aa_seq + genetic_code[codon_list[z]]
                    z = z + 1
                aa_list.append(aa_seq)
    return aa_list

def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    """
    sequence = sequence.upper()
    original_list = list(sequence)
    original_list.reverse()
    reverse = "".join(original_list)
    return reverse

def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    """
    sequence = sequence.upper()
    complements = {"A": "U", "U": "A", "C": "G", "G": "C"}
    sequence_list = list(sequence)
    new_seq = ""
    for elem in sequence_list:
        new_seq += complements[elem]
    return new_seq

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    """
    sequence = sequence.upper()
    backwards = get_reverse(sequence)
    RC = get_complement(backwards)
    return RC

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
    FWD = rna_sequence
    RC = reverse_and_complement(rna_sequence)
    FWD_translations = get_all_translations(FWD, genetic_code)
    RC_translations = get_all_translations(RC, genetic_code)
    all_translations = FWD_translations + RC_translations
    if all_translations == []:
        longest_peptide = ""
    else:
        longest_peptide = max(all_translations, key=len)
    return longest_peptide

if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
