# Script: Prior to running the pipeline in nextflow, we generate FASTA files for each gene family
#         based on file containing gene family assignments.
#         IDs in each FASTA file are in format >{specI_id}__{gene_id}, name of file corresponds to family ID.
# Usage:  python generate_fasta_files_for_pipeline.py gene_family_assignment output_dir

import sys


import pandas as pd
from ete3 import Tree
from skbio import DNA, Sequence
from skbio.io import read as skread

family_assignment_table = sys.argv[1]  # (first column - specI cluster, second column - gene id, third column - gene family id).
output_dir = sys.argv[2]

species_tree_file = "/mnt/mnemo4/maria/horizontal_gene_transfer/species_tree_progenomes_v2.2_no_chimera_reps_only_mapping_to_otus.nwk"
gene_sequence_dump = "/mnt/mnemo4/maria/horizontal_gene_transfer/pangenomes_95nr_gene_sequences.fasta"

# Dictionary which contains codes for special characters.
sp_char_codes = {"-": "mns", ".": "fstp", ":": "coln", "__": "uscr"}

# Uploading mapping between specI clusters, gene ids, and gene family ids.
map_gene_to_gene_family = pd.read_csv(family_assignment_table, sep = None, header = None, engine = 'python')
map_gene_to_gene_family.columns = ["speci_cluster", "gene_id", "family_id"]
map_gene_to_gene_family.index = map_gene_to_gene_family["gene_id"]
# Removing colon symbol from family ids, as it can cause trouble in the pipeline (unlike other special characters).
map_gene_to_gene_family["family_id"] = map_gene_to_gene_family["family_id"].apply(lambda x: x.replace(":", ""))

print("Uploaded gene family assignment table. Total size:", len(map_gene_to_gene_family))

# Uploading the tree in which all specI clusters not mapping to OTUs have been removed.
species_tree = Tree(species_tree_file, format = 0)
mapped_speci_clusters = species_tree.get_leaf_names()

counter = 0
# Going through file containing all gene sequences and writing genes to file corresponding to gene family.
for sequence in skread(gene_sequence_dump, format = 'fasta'):

    gene_id = sequence.metadata["id"]
    if gene_id in map_gene_to_gene_family.index:
        speci_cluster = map_gene_to_gene_family.loc[gene_id, "speci_cluster"]
        # We only write genes from specI clusters that mapped to an 96% OTU.
        if speci_cluster in mapped_speci_clusters:
            # Replacing all special characters so they don't cause problems when running the software downstream.
            new_gene_id = gene_id
            for sp_char in sp_char_codes:
                new_gene_id = new_gene_id.replace(sp_char, sp_char_codes[sp_char])
            
            # Writing sequence to the file of the corresponding to gene family.
            with open("{}/{}.fna".format(output_dir, map_gene_to_gene_family.loc[gene_id, "family_id"]), "a+") as out_file:
                out_file.write(">{}__{}".format(speci_cluster.split("_")[-1], new_gene_id) + "\n")
                out_file.write(str(sequence) + "\n")