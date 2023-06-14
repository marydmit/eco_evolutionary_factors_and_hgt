# Script: Takes a FASTA sequence that contains specI cluster and gene ids.  
# If a specI cluster contains genes that come from multiple genomes, will select
# most similar sequence to those in other specI clusters and only keep that for
# multiple sequence alignment and tree generation.
# Usage: python remove_duplicate_artefacts.py input_file_name output_file_name log_file_name

import sys
import shutil

import numpy as np
import pandas as pd

from skbio import DNA, Protein
from skbio import read as skread 
from skbio.alignment import local_pairwise_align_ssw


def check_if_has_duplicates(fasta_file_name):
    """
    Function: Opens indicated fasta file, checks if any specI clusters contain
              more than one genome
    Input:    fasta_file_name - string, path to FASTA file for consideration.
    Output:   list, all specI cluster ids that contain more than one genome.
    """

    all_sequences = skread(fasta_file_name, format = "fasta")
    
    # Storing all genome ids extracted from FASTA file.
    all_genome_ids = []
    all_speci_cluster_ids = []
    
    for sequence in all_sequences:
        seq_id = sequence.metadata["id"]
        
        speci_cluster = seq_id.split("_")[0]
        all_speci_cluster_ids.append(speci_cluster)
        
        genome = ".".join(seq_id.split(".")[:2]).split("_")[-1]
        all_genome_ids.append(genome)
        
    # Now checking whether there are any specI clusters with more than one genome.
    speci_to_genome = pd.DataFrame({"speci": all_speci_cluster_ids, "genome": all_genome_ids})
    genome_counts = speci_to_genome.groupby("speci").genome.nunique()
    
    return genome_counts.loc[genome_counts > 1].index.values


def import_fasta_with_duplicates(fasta_file_name, duplicate_cluster_ids):
    """
    Function: Opens indicated fasta file, imports all specI cluster sequences, 
              separates specI clusters with more than one genome for the others.
    Input:    fasta_file_name - string, path to FASTA file for import.
              duplicate_cluster_ids - list, specI clusters with multiple genomes.
    Output:   tuple, first object is a list of all sequences from specI clusters
              containing a single genome. Second object is dictionary, where keys
              specI cluster id and values are sequences from the specI cluster.
    """

    all_sequences = skread(fasta_file_name, format = "fasta")
    
    # We will store specI clusters with a single genome and those with several 
    # genomes in separate objects.
    single_genomes = []
    
    multiple_genomes = {}
    for cluster_id in duplicate_cluster_ids:
        multiple_genomes[cluster_id] = []
    
    # Checking for each sequence whether it is part of a specI cluster with
    # multiple genomes.
    for sequence in all_sequences:
        speci_cluster = sequence.metadata["id"].split("_")[0]
        if speci_cluster not in duplicate_cluster_ids:
            single_genomes.append(sequence)
        else:
            multiple_genomes[speci_cluster].append(sequence)
    
    return (single_genomes, multiple_genomes)


def choose_best_scoring_genome_seqs(sequences_to_choose_from, sequences_to_compare_to):
    """
    Function: Selects a sequence out of a list based on the highest average
              alignment score when compared to provided sequence list.
    Input:    sequences_to_choose_from - list of Sequence objects.
              sequences_to_compare_to - list of Sequence objects.
    Output:   list, the selected sequence and all sequences that come from
              the same genome.
    """

    sequence_scores = {}
    
    for sequence_1 in sequences_to_choose_from:
        sequence_scores[sequence_1.metadata["id"]] = 0
        
        # Comparing selected sequence to all sequences from specI clusters with a single genome.
        for sequence_2 in sequences_to_compare_to:
            sequence_scores[sequence_1.metadata["id"]] += local_pairwise_align_ssw(DNA(sequence_1), DNA(sequence_2))[1]
        
        # Calculating average alignment score.
        sequence_scores[sequence_1.metadata["id"]] = sequence_scores[sequence_1.metadata["id"]]/len(sequences_to_compare_to)
        
    # Selecting sequence with highest average alignment score.
    selected_sequence = max(sequence_scores, key = sequence_scores.get)
    
    # If selected genome happens to have more than one copy of the gene, selecting those as well.
    selected_genome = ".".join(selected_sequence.split(".")[:2]).split("_")[-1]
    selected_sequences = []
    
    for sequence in sequences_to_choose_from:
        if selected_genome in sequence.metadata["id"]:
            selected_sequences.append(sequence)
    
    return selected_sequences


def choose_pseudo_random_genome_seqs(sequences_to_choose_from):
    """
    Function: Selects a sequence out of a list of provided sequences based
              on whichever gene id is the smallest.
    Input:    sequences_to_choose_from - list of Sequence objects.
    Output:   list, the selected sequence and all sequences that come from
              the same genome.
    """
    
    # Extracting all sequence IDs.
    sequence_ids = [x.metadata["id"] for x in sequences_to_choose_from]

    # Sorting sequence IDs alphabetically and choosing whichever is the first.
    sequence_ids.sort()
    selected_sequence_id = sequence_ids[0]
    
    # If selected genome happens to have more than one copy of the gene, selecting those as well.
    selected_genome = ".".join(selected_sequence_id.split(".")[:2]).split("_")[-1]
    selected_sequences = []
    
    for sequence in sequences_to_choose_from:
        if selected_genome in sequence.metadata["id"]:
            selected_sequences.append(sequence)
    
    return selected_sequences


def main():

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    log_file_name = sys.argv[3]

    log_file = open(log_file_name, "w")

    # First checking if considered gene family has duplicates.
    duplicate_cluster_list = check_if_has_duplicates(input_file)

    # If no specI clusters with duplicate artefacts detected, no filtering required.
    if len(duplicate_cluster_list) == 0:
        log_file.write("No duplicate artefacts detected, retaining all sequences.\n")
        shutil.copy(input_file, output_file)
        exit()

    # Otherwise, we import the file and perform gene selection.
    single_clusters, multi_clusters = import_fasta_with_duplicates(input_file, duplicate_cluster_list)

    # In rare cases, all specI clusters within the gene family have duplicate
    # artefacts. For these, we select the sequences pseudorandomly based on sequence id.
    if len(single_clusters) == 0:
        log_file.write("No specI clusters without duplicate artefacts, choosing sequences based on ID only.\n")
        log_file.write("Starting from {} specI clusters containing {} sequences in total...\n".format(len(multi_clusters), \
            sum([len(x) for x in multi_clusters.values()])))
        
        final_selection = []
        for speci_cluster in multi_clusters:
            final_selection.extend(choose_pseudo_random_genome_seqs(multi_clusters[speci_cluster]))
        
        log_file.write("After filtering, {} sequences remain.\n".format(len(final_selection)))
    
    # Otherwise, choosing the best scoring sequence.
    else:
        log_file.write("Sequences available for comparison: {}.\n".format(len(single_clusters)))
        log_file.write("Starting from {} specI multi-genome clusters containing {} sequences...\n".format(len(multi_clusters), \
            sum([len(x) for x in multi_clusters.values()])))

        final_selection = single_clusters.copy()
        for speci_cluster in multi_clusters:
            final_selection.extend(choose_best_scoring_genome_seqs(multi_clusters[speci_cluster], single_clusters))

        log_file.write("During filtering, {} sequences were selected.\n".format(len(final_selection)-len(single_clusters)))
        log_file.write("Sequences remaining in total: {}.\n".format(len(final_selection)))
    
    log_file.write("Writing results to file...\n")
    log_file.close()

    with open(output_file, "w") as out_file:
        # Writing all selected sequences into a new file in FASTA format.
        for sequence in final_selection:
            out_file.write(">{}\n".format(sequence.metadata["id"]))
            out_file.write(str(sequence) + "\n")


if __name__ == "__main__":
    main()