# Script: Takes reconciliation files from Ranger-DTL and aggregates 
# results in a transfer matrix across different optimal roots.
# Usage: python aggregate_ranger_results.py prefix_ranger_dtl_files output_directory

from itertools import product
import glob
import random
import re
import string
import sys

from ete3 import Tree
import pandas as pd
import numpy as np

def get_info_from_ranger_file(file_prefix):
    
    in_file_trees = file_prefix
    in_file_aggregated = file_prefix + ".txt"
    
    flag_gene_tree = False
    flag_species_tree = False
    
    transfer_nodes = {}
    transfer_recipients = {}
    transfer_mappings = {}

    # Extracting structure of the species and gene trees.
    with open(in_file_trees) as in_file:

        for line in in_file:
            # Next line after the words species tree is always the species tree structure.
            if line.startswith("Species Tree:"):
                flag_species_tree = True
                continue
            if flag_species_tree:
                species_tree = Tree(line.strip(), format = 8)
                flag_species_tree = False

            # Next line after the words gene tree is always the tree structure.
            # The gene tree always comes second after species tree.
            if line.startswith("Gene Tree:"):
                flag_gene_tree = True
                continue
            if flag_gene_tree:
                gene_tree = Tree(line.strip(), format = 8)
                break
    
    # Extracting the names of all internal nodes that are labelled as Transfer
    # and their scores.
    with open(in_file_aggregated) as in_file:
        first = True
        for line in in_file:
            # The total number of files that went into the aggregation is
            # indicated in the first line.
            if first:
                total_score = int(line.split()[1])
                first = False
                
            # During reconciliations, RANGER-DTL automatically assigns internal node id
            # starting with m and with integer.
            if line.startswith("m") and line[1].isdigit():
                transfer_score = int(re.findall(r"\ Transfers = (.*)\], \[Most Frequent mapping", line)[0])
                if transfer_score != 0:
                    node_name = line.split(" = ")[0]
                    transfer_nodes[node_name] = transfer_score/total_score
                    
                    # Figuring out most likely mapping and its score.
                    likely_mapping = re.findall(r"\ mapping --> (.*) times],", line)[0].split(", ")[0]
                    if likely_mapping.startswith("n"):  # If transfer maps to internal node in species tree, take all descendant leaves.
                        likely_mapping = "LCA[{}]".format(", ".join([species for species in sorted((species_tree&likely_mapping).get_leaf_names()) if species.startswith("Cluster")]))
                    mapping_score = int(re.findall(r"\ mapping --> (.*) times],", line)[0].split(", ")[1])
                    transfer_mappings[node_name] = [likely_mapping, mapping_score/total_score]

                    
                    # Figuring out recipient and recipient score.
                    likely_recipient = re.findall(r"\ recipient --> (.*) times].", line)[0].split(", ")[0]
                    if likely_recipient.startswith("n"):  # If transfer to internal node in species tree, take all descendant leaves.
                        likely_recipient = "LCA[{}]".format(", ".join([species for species in sorted((species_tree&likely_recipient).get_leaf_names()) if species.startswith("Cluster")]))
                    recipient_score = int(re.findall(r"\ recipient --> (.*) times].", line)[0].split(", ")[1])
                    transfer_recipients[node_name] = [likely_recipient, recipient_score/total_score] 
                    
    return (gene_tree, transfer_nodes, transfer_recipients, transfer_mappings)

def descendant_names_with_lca(node_in_collapsed_tree, complete_tree):
    descendants_1 = node_in_collapsed_tree.children[0].get_leaf_names()
    descendants_2 = node_in_collapsed_tree.children[1].get_leaf_names()
    
    for descendants in [descendants_1, descendants_2]:
        for i in range(len(descendants)):
            if not descendants[i].startswith("Cluster"):
                descendants[i] = "LCA[{}]".format(", ".join(sorted((complete_tree&descendants[i]).get_leaf_names())))
    
    return(descendants_1, descendants_2)

def get_transfer_events(gene_tree, transfer_node_dict, recipient_dict, mapping_dict, record_no_transfers = False):
    
    # We will output the transfers as three lists, containing the
    # two parties participating the transfer and the transfer score.
    transfer_party_1 = []
    transfer_party_2 = []
    transfer_scores = []

    # To keep node identification consistent, we also record mappings and recipients in
    # the same order as transfer events.
    mapping = []
    mapping_scores = []
    recipient = []
    recipient_scores = []

    # In the first step, we need to collapse the gene tree in such a way
    # that we can obtain the last common ancestor for clades with
    # high-level transfers.
    for node in gene_tree.iter_descendants():
        # If node is labelled as transfer, we annotate all its descendants
        # with this transfer node.
        if node.name in transfer_node_dict:
            node_descendants = node.get_leaf_names()
            for leaf in node_descendants:
                try:
                    (gene_tree&leaf).transfers = (gene_tree&leaf).transfers + " " + node.name
                except:
                    (gene_tree&leaf).add_feature("transfers", node.name)

    # The collapse function will collapse all internal nodes that are not
    # labelled as transfers but contain leaves participating in the same 
    # transfers.
    node2transfers = gene_tree.get_cached_content(store_attr = "transfers")
    
    # For future reference, we may want to store those gene pairs which have
    # no history of horizontal gene transfers. We therefore record these prior
    # to collapsing the tree.
    if record_no_transfers:
        nodes_without_transfers = {"descendant_1":[], "descendant_2": []}
        for node in node2transfers:
            if node2transfers[node] == {None}:  # if no transfers recorded
                if node.name.startswith("m"): # if internal node
                    descendants_1 = node.children[0].get_leaf_names()
                    descendants_2 = node.children[1].get_leaf_names()
                    nodes_without_transfers["descendant_1"].extend([item for sublist in [[x] * len(descendants_2) for x in descendants_1] for item in sublist])
                    nodes_without_transfers["descendant_2"].extend(descendants_2 * len(descendants_1))
        nodes_without_transfers = pd.DataFrame(nodes_without_transfers)

    def collapsed_leaf(node):
        if len(node2transfers[node]) == 1 and node.name not in transfer_node_dict:
            return True
        else:
            return False
    collapsed_tree = Tree(gene_tree.write(is_leaf_fn = collapsed_leaf, \
                                          features = ["name"], format_root_node = True))
    
    # We will now iterate through the collapsed tree and extract
    # which genes participated in a transfer and transfer confidence score.
    for node in collapsed_tree.iter_descendants():
        if node.name in transfer_node_dict:
            all_descendants = descendant_names_with_lca(node, gene_tree)
            transfer_party_1.append(set(all_descendants[0]))
            transfer_party_2.append(set(all_descendants[1]))
            transfer_scores.append(transfer_node_dict[node.name])

            mapping.append(mapping_dict[node.name][0])
            mapping_scores.append(mapping_dict[node.name][1])
            recipient.append(recipient_dict[node.name][0])
            recipient_scores.append(recipient_dict[node.name][1])

    if not record_no_transfers:
        return(transfer_party_1, transfer_party_2, transfer_scores, mapping, mapping_scores, recipient, recipient_scores, [])
    else:
        return(transfer_party_1, transfer_party_2, transfer_scores, mapping, mapping_scores, recipient, recipient_scores, nodes_without_transfers)


def return_special_chars(string):
    # Mapping for all special characters that where substituted.
    sp_char_codes = {"mns": "-", "fstp": ".", "coln": ":", "uscr": "__"}
    
    for char in sp_char_codes:
        string = string.replace(char, sp_char_codes[char])

    return string


def main():

    # Files provided as arguments 
    ranger_dtl_file_prefix = sys.argv[1]
    output_directory = sys.argv[2]
    
    # The name of each gene family is contained within the file name.
    family_id = ranger_dtl_file_prefix.split("/")[-1]

    # Collecting all performed reconciliations.
    ranger_dtl_files = glob.glob("{}*.txt".format(ranger_dtl_file_prefix))
    if len(ranger_dtl_files) == 0:
        print("No RANGER output files detected with given prefix!")
        exit()

    # Creating lists for each piece of information that will go into the final table.
    node_id = []
    present_in_roots = []
    transfer_party_1 = []
    transfer_party_2 = []
    transfer_score = []
    
    node_id_mapping = []
    mapping = []
    mapping_score = []

    node_id_recipient = []
    recipient = []
    recipient_score = []
    
    tracked_events = []

    for in_file in ranger_dtl_files:

        # For each specific root, extracting the transfer events.
        gene_tree, transfer_node_dict, recipient_dict, mapping_dict = get_info_from_ranger_file(in_file[:-4])
        tp1, tp2, ts, mp, mps, rc, rcs, pairs_without_transfers = get_transfer_events(gene_tree, transfer_node_dict, recipient_dict, mapping_dict, True)
        
        # Checking if each transfer event is already in the list.
        # If not, creating new node id in the list.
        # If yes, simply adding the transfer score to the respective event.
        for i in range(len(tp1)):
            if (tp1[i], tp2[i]) in tracked_events:
                try:
                    loc_index = np.where((np.array(transfer_party_1) == tp1[i]) & (np.array(transfer_party_2) == tp2[i]))[0][0]
                except:
                    loc_index = np.where((np.array(transfer_party_2) == tp1[i]) & (np.array(transfer_party_1) == tp2[i]))[0][0]
                transfer_score[loc_index] += ts[i]
                present_in_roots[loc_index] += 1

                # For information on mappings, checking if mapping node in current event has been recorded preivously.
                # Otherwise, adding a new possible mapping node to the list.
                if mp[i] in mapping[loc_index]:
                    map_index = mapping[loc_index].index(mp[i])
                    mapping_score[loc_index][map_index] += mps[i]
                else:
                    node_id_mapping[loc_index].append(node_id[loc_index])
                    mapping[loc_index].append(mp[i])
                    mapping_score[loc_index].append(mps[i])

                # Repeating same procedure as for mapping for the recipients.
                if rc[i] in recipient[loc_index]:
                    rec_index = recipient[loc_index].index(rc[i])
                    recipient_score[loc_index][rec_index] += rcs[i]
                else:
                    node_id_recipient[loc_index].append(node_id[loc_index])
                    recipient[loc_index].append(rc[i])
                    recipient_score[loc_index].append(rcs[i])
            else:
                # For each new transfer node, generating random sequence as an id.
                node_id.append(''.join(random.choices(string.ascii_lowercase + string.digits, k = 8)))
                transfer_party_1.append(tp1[i])
                transfer_party_2.append(tp2[i])
                transfer_score.append(ts[i])

                node_id_mapping.append([node_id[-1]])
                mapping.append([mp[i]])
                mapping_score.append([mps[i]])

                node_id_recipient.append([node_id[-1]])
                recipient.append([rc[i]])
                recipient_score.append([rcs[i]])
                
                present_in_roots.append(1)

                # Also tracking this specific transfer pair in case it occurs in another root.
                tracked_events.append((tp1[i], tp2[i]))
    
    # Recording pairs without transfers in their own table prior to proceeding with the transfers.
    pairs_without_transfers["family_id"] = family_id
    pairs_without_transfers["speci_1"] = pairs_without_transfers["descendant_1"].apply(lambda x: x.split("__")[0])
    pairs_without_transfers["genome_1"] = pairs_without_transfers["descendant_1"].apply(lambda x: ".".join(x.split("__")[1].split("fstp")[:2]))
    pairs_without_transfers["gene_1"] = pairs_without_transfers["descendant_1"].apply(lambda x: ".".join(x.split("__")[1].split("fstp")[2:]))

    pairs_without_transfers["speci_2"] = pairs_without_transfers["descendant_2"].apply(lambda x: x.split("__")[0])
    pairs_without_transfers["genome_2"] = pairs_without_transfers["descendant_2"].apply(lambda x: ".".join(x.split("__")[1].split("fstp")[:2]))
    pairs_without_transfers["gene_2"] = pairs_without_transfers["descendant_2"].apply(lambda x: ".".join(x.split("__")[1].split("fstp")[2:]))

    pairs_without_transfers["gene_1"] = pairs_without_transfers["gene_1"].apply(lambda x: return_special_chars(x))
    pairs_without_transfers["gene_2"] = pairs_without_transfers["gene_2"].apply(lambda x: return_special_chars(x))

    pairs_without_transfers = pairs_without_transfers.iloc[:,2:]

    pairs_without_transfers.to_csv("{}/{}_no_transfers.csv".format(output_directory, family_id), index = False)

    # If no transfers detected in this gene_family, writing message to file with transfer matrix.
    # Creatingblank files for transfer mappings and recipients. Exiting script.
    if len(transfer_party_1) == 0:
        with open("{}/{}.csv".format(output_directory, family_id), "w") as out_file:
            out_file.write("No transfers detected.\n")
        open("{}/{}_mappings.csv".format(output_directory, family_id), "w").close()
        open("{}/{}_recipients.csv".format(output_directory, family_id), "w").close()
        exit()

    # Now we have some transfers that involve several options where it's unclear between which to parties
    # exactly the transfer happened. We therefore generate a cartesian product between all these options,
    # expanding the transfer table.
    full_node_ids = []
    full_roots_present = []
    full_transfer_party_1 = []
    full_descendants_1 = []
    full_transfer_party_2 = []
    full_descendants_2 = []
    full_transfer_scores = []

    for i in range(len(transfer_party_1)):
        transfer_party_product = product(transfer_party_1[i], transfer_party_2[i])

        for entry in transfer_party_product:
            descendants_1 = entry[0].replace("LCA[", "").replace("]", "").split(", ")
            descendants_2 = entry[1].replace("LCA[", "").replace("]", "").split(", ")
            total_descendants = len(descendants_1) * len(descendants_2)

            # Checking whether it's a single species or the last common ancestor of multiple
            # species.
            if len(descendants_1) > 1:
                final_transfer_party_1 = "LCA[" + ", ".join(sorted(set([x.split("__")[0] for x in descendants_1]))) + "]"
                full_transfer_party_1.extend([final_transfer_party_1] * total_descendants)
            else:
                full_transfer_party_1.extend([entry[0].split("__")[0]] * total_descendants)
            full_descendants_1.extend(descendants_1 * len(descendants_2))
            
            if len(descendants_2) > 1:
                final_transfer_party_2 = "LCA[" + ", ".join(sorted(set([x.split("__")[0] for x in descendants_2]))) + "]"
                full_transfer_party_2.extend([final_transfer_party_2] * total_descendants)
            else:
                full_transfer_party_2.extend([entry[1].split("__")[0]] * total_descendants)
            [full_descendants_2.extend([x] * len(descendants_1)) for x in descendants_2]

            full_node_ids.extend([node_id[i]] * total_descendants)
            full_roots_present.extend([present_in_roots[i]] * total_descendants)
            full_transfer_scores.extend([transfer_score[i]] * total_descendants)

    # Creating dataframe in which we collect all information on transfers.
    transfer_df = pd.DataFrame({"node_id": full_node_ids, "roots_present": full_roots_present, \
                                "roots_total": len(ranger_dtl_files), "transfer_node_1": full_transfer_party_1, \
                                "descendant_1": full_descendants_1, "transfer_node_2": full_transfer_party_2, \
                                "descendant_2": full_descendants_2, "score": full_transfer_scores})

    # Calculating transfer confidence score in percentages.
    transfer_df["score"] = round(100 * (transfer_df["score"] / len(ranger_dtl_files)), 2)

    # We split the long gene id into separate columns representing specI clusters, genome and gene ids.
    transfer_df["speci_1"] = transfer_df["descendant_1"].apply(lambda x: x.split("__")[0])
    transfer_df["genome_1"] = transfer_df["descendant_1"].apply(lambda x: ".".join(x.split("__")[1].split("fstp")[:2]))
    transfer_df["gene_1"] = transfer_df["descendant_1"].apply(lambda x: ".".join(x.split("__")[1].split("fstp")[2:]))
    
    transfer_df["speci_2"] = transfer_df["descendant_2"].apply(lambda x: x.split("__")[0])
    transfer_df["genome_2"] = transfer_df["descendant_2"].apply(lambda x: ".".join(x.split("__")[1].split("fstp")[:2]))
    transfer_df["gene_2"] = transfer_df["descendant_2"].apply(lambda x: ".".join(x.split("__")[1].split("fstp")[2:]))
    
    # Adding gene family id - this is useful when we'll dump all transfer events in one file.
    transfer_df["family_id"] = family_id

    # Rearranging columns in a convenient way prior to exporting.
    transfer_df_final = transfer_df.copy()
    transfer_df_final = transfer_df_final[["family_id", "node_id", "roots_present", "roots_total", \
                                           "transfer_node_1", "speci_1", "genome_1", "gene_1", \
                                           "transfer_node_2", "speci_2", "genome_2", "gene_2", "score"]]
    transfer_df_final = transfer_df_final.sort_values(by = ["roots_present", "node_id", "transfer_node_1", "transfer_node_2", "score"], ascending = False)

    # Substituting special characters back into gene names.
    transfer_df_final["gene_1"] = transfer_df_final["gene_1"].apply(lambda x: return_special_chars(x))
    transfer_df_final["gene_2"] = transfer_df_final["gene_2"].apply(lambda x: return_special_chars(x))
    
    transfer_df_final.to_csv("{}/{}.csv".format(output_directory, family_id), index = False)

    # In a separate, we store information on likely mappings.
    mapping_df = pd.DataFrame({"family_id": family_id, "node_id": list(np.concatenate(node_id_mapping)), "mapping": list(np.concatenate(mapping)), "score": list(np.concatenate(mapping_score))})

    # Dividing mapping scores to account for possible multiple roots.
    mapping_df["score"] = round(100 * (mapping_df["score"] / len(ranger_dtl_files)), 2)
    mapping_df = mapping_df.sort_values(by = ["node_id", "score"], ascending = False)
    mapping_df.to_csv("{}/{}_mappings.csv".format(output_directory, family_id), index = False)

    # Finally, we store information on likely recipients in another separate file.
    recipient_df = pd.DataFrame({"family_id": family_id, "node_id": list(np.concatenate(node_id_recipient)), "recipient": list(np.concatenate(recipient)), "score": list(np.concatenate(recipient_score))})
    
    # Dividing recipient scores to account for possible multiple roots.
    recipient_df["score"] = round(100 * (recipient_df["score"] / len(ranger_dtl_files)), 2)
    recipient_df = recipient_df.sort_values(by = ["node_id", "score"], ascending = False)
    recipient_df.to_csv("{}/{}_recipients.csv".format(output_directory, family_id), index = False)

if __name__ == "__main__":
    main()