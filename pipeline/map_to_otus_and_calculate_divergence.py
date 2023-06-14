# Script: Takes a transfer table and adds OTU mapping and how taxonomically divergent are the species.
# Usage:  python map_to_otus_and_calculate_divergence.py transfer_table otu_mapping_table progenomes_taxonomy gtdb_taxonomy

import sys

import pandas as pd
import numpy as np
    
def extract_otu_divergence(otu_id_1, otu_id_2):
    """
    Function: Takes two OTU IDs and returns level at which they diverge.
    Input:    otu_id_1, string
              otu_id_2, string
    Output:   divergence, string indicating taxonomic level
    """
        
    # Some specI clusters were not mapped to any OTU, so have to make sure these
    # are not assigned to a level.
    if otu_id_1 == "unmapped" or otu_id_2 == "unmapped":
        return "NA"

    # Splitting OTU ids into the different levels.
    otu_1_list = [x.split(";") for x in otu_id_1.split(", ")]
    otu_2_list = [x.split(";") for x in otu_id_2.split(", ")]
    
    # Return lowest taxonomic level at which divergence is detected.
    divergence_levels = ["kingdom", "90", "96", "97", "98"]
    divergence = []

    for otu_1 in otu_1_list:
        for otu_2 in otu_2_list:
            for i in range(min(len(otu_1), len(otu_2))):
                if otu_1[i] != otu_2[i]:
                    divergence.append(i)
                    break

    if len(divergence) == 0:
        return "NA"
    else:
        return divergence_levels[max(divergence)] 


def extract_speci_divergence(speci_1, speci_2, taxonomy_table):
    
    """
    Function: Takes two speci cluster IDs, looks up taxonomy in table and returns
              level at which they diverge.
    Input:    speci_1, string
              speci_2, string
              taxonomy_table, pandas DataFrame, contains specI cluster and taxonomy assigned to it.
    Output:   string indicating taxonomic level
    """
        
    # Extracting taxonomic information for specific specI clusters. The specI cluster is indicated
    # in the first column of the table.
    try:
        taxonomy_1 = taxonomy_table.loc[taxonomy_table.iloc[:,0] == speci_1, :].iloc[0, 1:]
        taxonomy_2 = taxonomy_table.loc[taxonomy_table.iloc[:,0] == speci_2, :].iloc[0, 1:]

        if taxonomy_1.equals(taxonomy_2):
            return "no divergence"
        elif taxonomy_1.isna().all() or taxonomy_2.isna().all():
        # Some specI clusters failed GTDB quality test, so no taxonomy is available for them.
            return "NA"
        else:
            return taxonomy_1.loc[~taxonomy_1.eq(taxonomy_2)].index[0]

    except:
        return "NA"

def main():

    # Importing the files provided.
    with open(sys.argv[1], "r") as in_file:
        if in_file.readline().startswith("No transfers"):
            exit()
    transfer_df = pd.read_csv(sys.argv[1])
    
    mapping_specI_clusters_otus = pd.read_csv(sys.argv[2], sep = None, engine = "python")
    
    # In case the specI_v3 prefix is included, making sure we add it to our cluster ids.
    if mapping_specI_clusters_otus.iloc[0, 0].startswith("specI_v3_"):
        transfer_df["speci_1_to_map"] = "specI_v3_" + transfer_df["speci_1"]
        transfer_df["speci_2_to_map"] = "specI_v3_" + transfer_df["speci_2"]
    
    specI_cluster_progenomes_taxonomy = pd.read_csv(sys.argv[3], sep = None, engine = "python").loc[:, ["specI_cluster", "kingdom", "phylum", "class", "order", "family", "genus", "species"]]

    specI_cluster_gtdb_taxonomy = pd.read_csv(sys.argv[4], sep = None, engine = "python").loc[:, ["taxID.sampleID_Fr12", "kingdom", "phylum", "class", "order", "family", "genus", "species"]] 

    # First, mapping specI clusters to OTUs.
    transfer_df["otu_1"] = transfer_df["speci_1_to_map"].apply(lambda x: ", ".join(set(mapping_specI_clusters_otus.loc[mapping_specI_clusters_otus.iloc[:, 0] == x, "OTU97"].values)))
    transfer_df["otu_2"] = transfer_df["speci_2_to_map"].apply(lambda x: ", ".join(set(mapping_specI_clusters_otus.loc[mapping_specI_clusters_otus.iloc[:, 0] == x, "OTU97"])))

    # Based on the OTUs, determining what is the level of divergence between two species.
    transfer_df["otu_diff"] = transfer_df.apply(lambda x: extract_otu_divergence(x["otu_1"], x["otu_2"]), axis = 1)
    
    # In addition, determining species divergence based on the two taxonomies we have.
    transfer_df["speci_diff"] = transfer_df.apply(lambda x: extract_speci_divergence(x["speci_1_to_map"], x["speci_2_to_map"], specI_cluster_progenomes_taxonomy), axis = 1)
    transfer_df["gtdb_diff"] = transfer_df.apply(lambda x: extract_speci_divergence(x["genome_1"], x["genome_2"], specI_cluster_gtdb_taxonomy), axis = 1)
    
    # Rearranging columns in a convenient way prior to exporting.
    transfer_df_final = transfer_df.copy()
    transfer_df_final = transfer_df_final[["family_id", "node_id", "roots_present", "roots_total", \
                                           "transfer_node_1", "speci_1", "otu_1", "genome_1", "gene_1", \
                                           "transfer_node_2", "speci_2", "otu_2", "genome_2", "gene_2", \
                                           "score", "speci_diff", "gtdb_diff", "otu_diff"]]

    transfer_df_final.to_csv(sys.argv[1], index = False)

if __name__ == "__main__":
    main()