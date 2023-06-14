# Script: Takes output from MAFFT and produces a sequence percentage identity matrix.
# Usage:  python calculate_sequence_identity.py mafft_output_file identity_matrix_file

from math import floor
import sys

import numpy as np
import pandas as pd

from Bio import AlignIO

def return_special_chars(string):
    # Mapping for all special characters that where substituted.
    sp_char_codes = {"mns": "-", "fstp": ".", "coln": ":", "uscr": "__"}
    
    for char in sp_char_codes:
        string = string.replace(char, sp_char_codes[char])

    return string

def generate_identity_matrix(alignment):
    """
    Function: produces identity matrix from sequence alignment.
    Input:    alignment, a Bio.Align.MultipleSeqAlignment object.
    Output:   pandas DataFrame with three columns: gene 1 ID, gene 2 ID and percentage identity.
    """

    gene_1 = []
    gene_2 = []
    non_gap_1 = []  # how many non-gap characters are there in the first sequence?
    non_gap_2 = []  # how many non-gap characters are there in the second sequence?
    perc_identity = []
    perc_overlap = []  # how much of the alignment does the overlap between two sequences cover?
    perc_overlap_short = []  # how much of the shorter sequence does the alignment cover?

    for sequence_1 in range(len(alignment) - 1):
        for sequence_2 in range(sequence_1 + 1, len(alignment)):
            
            total = 0
            matches = 0
            len_seq_1 = 0
            len_seq_2 = 0

            for base in range(alignment.get_alignment_length()):
                if alignment[sequence_1][base] != "-" and alignment[sequence_2][base] != "-":
                    # For identity, only counting positions where both sequences don't have a gap.
                    total += 1
                    len_seq_1 += 1
                    len_seq_2 += 1
                    if alignment[sequence_1][base] == alignment[sequence_2][base]:
                            matches += 1
                elif alignment[sequence_1][base] != "-":
                    len_seq_1 += 1
                elif alignment[sequence_2][base] != "-":
                    len_seq_2 += 1

            if total == 0:
                # If sequences don't overlap in alignment, not possible to calculate identity.
                identity = np.nan
            else:
                identity = round(100 * (matches/total), 2)
            
            # Had to remove special characters for alignment and tree generation, inserting them back.
            gene_id_1 = return_special_chars(alignment[sequence_1].id)
            gene_id_2 = return_special_chars(alignment[sequence_2].id)

            gene_1.extend([gene_id_1, gene_id_2])
            gene_2.extend([gene_id_2, gene_id_1])
            non_gap_1.extend([len_seq_1, len_seq_2])
            non_gap_2.extend([len_seq_2, len_seq_1])
            perc_identity.extend(2 * [identity])
            perc_overlap.extend(2 * [round(100 * total/alignment.get_alignment_length(), 2)])
            perc_overlap_short.extend(2 * [round(100 * total/min(len_seq_1, len_seq_2), 2)])
    
    result_df = pd.DataFrame({"gene_1": gene_1, "gene_2": gene_2, "len_seq_1": non_gap_1, "len_seq_2": non_gap_2, \
            "seq_identity": perc_identity, "seq_coverage": perc_overlap, "seq_coverage_short": perc_overlap_short})
    result_df = result_df.sort_values(by = "gene_1")
    
    return result_df

def main():
    alignment = AlignIO.read(sys.argv[1], "fasta")
    identity_matrix = generate_identity_matrix(alignment)
    identity_matrix.to_csv(path_or_buf = sys.argv[2], index = False)

if __name__ == "__main__":
    main()
