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

#OLD DOES NOT WORK
#make your sequence uppercase    
#   rna_sequence = rna_sequence.upper()
#split your sequence into a list where each nucleotide is a separate element
#   rna_list = list(rna_sequence)
#create a variable x = 0, then create the object rf1, and turn it into an empty list (rf1_list)
#   x = 0
#   rf1 = ""
#   rf1_list = list(rf1)
#create a while loop that will return a list of nucleotide triplets in reading frame 1 (rf1_list)
#   while x + 2 < len(rna_list):
#       triplet = rna_list[x] + rna_list[x+1] + rna_list[x+2]
#       rf1_list.append(triplet)
#       x = x + 3
#do the same thing for reading frame 2
#   rf2 = ""
#   rf2_list = list(rf2)
#   y = 1
#   while y + 2 < len(rna_list):
#       triplet = rna_list[y] + rna_list[y+1] + rna_list[y+2]
#       rf2_list.append(triplet)
#       y = y + 3
#and reading frame 3
#   rf3 = ""
#   rf3_list = list(rf3)
#   z = 2
#   while z + 2 < len(rna_list):
#       triplet = rna_list[z] + rna_list[z+1] + rna_list[z+2]
#       rf3_list.append(triplet)
#       z = z + 3
#everything up to here works
#for every codon in your reading frame 1 list, cycle through. if the codon codes for a stop codon, break the loop without adding a stop codon to your aa sequence.
#otherwise, add the appropriate amino acid for your sequence to the string "new1"
#   new1 = ""
#   for elem in rf1_list:
#       new1 = new1 + genetic_code[elem]
#take the string new1 and turn it into a tuple, where each amino acid is an item in the tuple.
#   frame1_tup = tuple(new1)
#create an empty variable k, then turn it into a list called frame1_aa_list
#   k = ""
#   frame1_aa_list = list(k)
#create a for loop that iterates the same number of times as the number of amino acids in your tuple
#if the index of frame1_tup at that particular number is "M", create a new variable d
#d is a string that contains the joined characters from M onwards, until the end of your amino acid tuple
#then, append d to your amino acid list for reading frame 1, and then make d an empty string again
#if the index of frame1_tup is NOT M, return and start over for the next num
#   for num in range(len(frame1_tup)):
#       if frame1_tup[num] == "M":
#           d = "".join(frame1_tup[num:])
#           frame1_aa_list.append(d)
#           d = ""
#   new2 = ""
#   for elem in rf2_list:
#       new2 = new2 + genetic_code[elem]
#   frame2_tup = tuple(new2)
#   m = ""
#   frame2_aa_list = list(m)
#   for num in range(len(frame2_tup)):
#       if frame2_tup[num] == "M":
#           e = "".join(frame2_tup[num:])
#           frame2_aa_list.append(e)
#           e = ""
#   new3 = ""
#   for elem in rf3_list:
#       new3 = new3 + genetic_code[elem]
#   frame3_tup = tuple(new3)
#   n = ""
#   frame3_aa_list = list(n)
#   for num in range(len(frame3_tup)):
#       if frame2_tup[num] == "M":
#           f = "".join(frame3_tup[num:])
#           frame3_aa_list.append(f)
#           f = ""
#   all_aa = frame1_aa_list + frame2_aa_list + frame3_aa_list
#   return all_aa

def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    """
    sequence.upper()
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
    pass


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
