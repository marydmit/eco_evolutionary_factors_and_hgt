# Script: After running the pipeline and obtaining a list of all transfer events, we count the number
#         of gene families a transfer event has been detected in for each OTU pair. We can adjust thresholds
#         on whether only transfer nodes with certain support or gene pairs of certain identities are counted.
# Usage: python count_transfer_events_per_otu_pair.py input_file output_prefix node_support seq_coverage seq_identity_min seq_identity_max

import sys

import pandas as pd
import numpy as np

def aggregate_taxonomic_divergence(divergence, genes_transferred):
    """
    Function: Takes list of OTU pair divergences and aggregates the number of genes transferred 
              at each taxonomic level.
    Input:    divergence - list of different OTU divergences (e.g. 96, 97)
              genes_transferred - list of genes transferred corresponding to each OTU divergenece (same order).
    Output:   data aggregated across different degrees of divergence. 
    """
    unique_divergence = list(set(divergence))
    unique_genes_transferred = [0]*len(unique_divergence)
    for i in range(len(unique_genes_transferred)):
        unique_genes_transferred[i] = sum(np.array(genes_transferred)[np.array(divergence) == unique_divergence[i]])

    if len(unique_divergence) == 1:
        return unique_divergence[0]

    elif len(unique_divergence) == 2 and "NA" in unique_divergence:
        unique_divergence.remove("NA")
        return unique_divergence[0]

    else:
        combined_data = sorted(list(zip(unique_genes_transferred, unique_divergence)), reverse = True)
        combined_data = ",".join(["_".join([str(div_level[0]), div_level[1]]) for div_level in combined_data])
        return(combined_data)

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

def main():
    
    input_file = sys.argv[1]
    output_prefix = sys.argv[2]
    min_node_support = float(sys.argv[3])
    min_seq_coverage = float(sys.argv[4])
    min_seq_distance = float(sys.argv[5])
    max_seq_distance = float(sys.argv[6])

    column_types = {"family_id": str, "node_id": str, "min_support": float, \
                    "speci_1": str, "otu_1": str, "genome_1": str, "gene_1": str, \
                    "speci_2": str, "otu_2": str, "genome_2": str, "gene_2": str, \
                    "speci_diff": str, "gtdb_diff": str, "otu_diff": str, \
                    "distance": float, "seq_identity": float, "seq_coverage": float}
    
    transfer_table_1 = pd.read_csv(input_file, sep = "\t", dtype = column_types)
    
    # Based on previous analysis, would eliminate artificially inflated counts if I filter out some genomes causing these problems.
    problematic_genomes = ["1932669.SAMEA43981918", "550.SAMEA3856670", "1916230.SAMN05904579", "1618951.SAMN03319549", "1802238.SAMN04314631", "74031.SAMN03943304", "663.SAMN04088202"]
    transfer_table_1 = transfer_table_1.loc[~(transfer_table_1["genome_1"].isin(problematic_genomes)) & ~(transfer_table_1["genome_2"].isin(problematic_genomes)), :]
    problematic_clusters = ["Cluster22", "Cluster2375", "Cluster4851", "Cluster5010", "Cluster5514"]
    transfer_table_1 = transfer_table_1.loc[~(transfer_table_1["speci_1"].isin(problematic_clusters)) & ~(transfer_table_1["speci_2"].isin(problematic_clusters)), :]

    # Making sure each transfer is included in both directions.
    transfer_table_2 = transfer_table_1.loc[: , ["family_id", "node_id", "min_support", \
                                                 "speci_2", "otu_2", "genome_2", "gene_2", \
                                                 "speci_1", "otu_1", "genome_1", "gene_1", \
                                                 "speci_diff", "gtdb_diff", "otu_diff", \
                                                 "distance", "gene_distance", "seq_identity", "seq_coverage"]]
    transfer_table_2.columns = transfer_table_1.columns
    transfer_table = pd.concat([transfer_table_1, transfer_table_2])
    # Prior to summarizing, need to change the NaN in gtdb_diff and speci_dif columns to strings.
    # Otherwise, these rows will be dropped during the summarization.
    transfer_table.loc[transfer_table["speci_diff"].isna() , "speci_diff"] = "NA"
    transfer_table.loc[transfer_table["gtdb_diff"].isna() , "gtdb_diff"] = "NA"

    # Some specI clusters map to more than one OTU. Therefore, exploding table prior to aggregation.
    transfer_table["otu_1"] = transfer_table["otu_1"].apply(lambda x: x.split(", "))
    transfer_table["otu_2"] = transfer_table["otu_2"].apply(lambda x: x.split(", "))
    transfer_table = transfer_table.explode("otu_1")
    transfer_table = transfer_table.explode("otu_2")

    # Removing rows from specI clusters mapping to the same OTU or specI clusters that are unmapped.
    transfer_table = transfer_table.loc[transfer_table["otu_1"] != transfer_table["otu_2"], :]
    transfer_table = transfer_table.loc[(transfer_table["otu_1"] != "unmapped") & (transfer_table["otu_2"] != "unmapped"), :]

    # Based on defined thresholds, removing pairs based on sequence identity and node support.
    transfer_table = transfer_table.loc[transfer_table["min_support"] >= min_node_support, :]
    transfer_table = transfer_table.loc[transfer_table["seq_coverage"] >= min_seq_coverage, :]
    transfer_table = transfer_table.loc[(transfer_table["gene_distance"] >= min_seq_distance) & (transfer_table["gene_distance"] <= max_seq_distance), :]

    # Aggregating data from all families to count the number of genes transferred.
    transfers_per_otu_pair = transfer_table.groupby(["otu_1", "otu_2", "speci_diff", "gtdb_diff"]).family_id.nunique().reset_index()  # counting genes transferred.
    transfers_per_otu_pair = transfers_per_otu_pair.groupby(["otu_1", "otu_2"])

    genes_transferred_aggregated = transfers_per_otu_pair.family_id.sum()
    gtdb_divergence_aggregated = transfers_per_otu_pair.apply(lambda x: aggregate_taxonomic_divergence(list(x["gtdb_diff"].values), list(x["family_id"].values)))
    speci_divergence_aggregated = transfers_per_otu_pair.apply(lambda x: aggregate_taxonomic_divergence(list(x["speci_diff"].values), list(x["family_id"].values)))

    genes_transferred_aggregated = transfers_per_otu_pair.family_id.sum()
    transfers_per_otu_pair = pd.concat([speci_divergence_aggregated, gtdb_divergence_aggregated, genes_transferred_aggregated], axis = 1).reset_index()
    transfers_per_otu_pair.columns = ["otu_1", "otu_2", "speci_diff", "gtdb_diff", "genes_transferred"]
    
    # Adding data on OTU divergence for completeness.
    transfers_per_otu_pair["otu_diff"] = transfers_per_otu_pair.apply(lambda x: extract_otu_divergence(x["otu_1"], x["otu_2"]), axis = 1)
    transfers_per_otu_pair = transfers_per_otu_pair.loc[:, ["otu_1", "otu_2", "otu_diff", "speci_diff", "gtdb_diff", "genes_transferred"]]

    # Writing resulting table to file.
    transfers_per_otu_pair.to_csv("{}_{:.2f}_{:.2f}_support_{:.2f}.csv".format(output_prefix, min_seq_distance, max_seq_distance, min_node_support), index = False)

if __name__ == "__main__":
    main()
