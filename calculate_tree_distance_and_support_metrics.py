# Script: Calculates tree size, distances and support values between gene pairs.
#         ! This script was added later in the analysis to filter out gene pairs in poorly
#           resolved regions of the tree.
# Usage:  python calculate_tree_metrics.py gene_tree_file


import sys

from ete3 import Tree
from numpy import median, mean
from pandas import DataFrame

def return_special_chars(string):
    # Mapping for all special characters that where substituted.
    sp_char_codes = {"mns": "-", "fstp": ".", "coln": ":", "uscr": "__"}
    
    for char in sp_char_codes:
        string = string.replace(char, sp_char_codes[char])

    return string

input_tree_file = sys.argv[1]
gene_tree = Tree(input_tree_file)

# To retrieve minimal support value per tree, we will need to root the gene tree first
# based on optimal roots found by RANGER-DTL.
rooted_tree_file = input_tree_file.replace("gene_trees", "optimal_roots").replace("nw", "optroots")
rooted_gene_trees = []
with open(rooted_tree_file) as file_with_gene_root:
    is_tree_structure = False
    for line in file_with_gene_root:
        if is_tree_structure:
            try:
                rooted_gene_trees.append(Tree(line, format = 9))
            except:
                break
                
        if line.startswith(" ---"):
            is_tree_structure = True

# Collecting all rooted trees in a list.
rooted_gene_trees_w_support = []
for rooted_gene_tree in rooted_gene_trees:
    copy_of_gene_tree = gene_tree.copy()
    
    # Deciding on outgroup based on rooted gene tree structure.
    selected_outgroup_1 = rooted_gene_tree.children[0].get_leaf_names()
    if len(selected_outgroup_1) == 1:
        selected_outgroup_node_1 = rooted_gene_tree&selected_outgroup_1[0]
    else:
        selected_outgroup_node_1 = copy_of_gene_tree.get_common_ancestor(selected_outgroup_1)

    selected_outgroup_2 = rooted_gene_tree.children[1].get_leaf_names()
    if len(selected_outgroup_2) == 1:
        selected_outgroup_node_2 = rooted_gene_tree&selected_outgroup_2[0]
    else:
        selected_outgroup_node_2 = copy_of_gene_tree.get_common_ancestor(selected_outgroup_2)

    if len(selected_outgroup_node_1) < len(selected_outgroup_node_2):
        selected_outgroup_node = selected_outgroup_node_1
    else:
        selected_outgroup_node = selected_outgroup_node_2
    if len(selected_outgroup_1) == 1:
        selected_outgroup_node = selected_outgroup_1[0]
    if len(selected_outgroup_2) == 1:
        selected_outgroup_node = selected_outgroup_2[0]

    # Rooting the gene tree.
    try:
        copy_of_gene_tree.set_outgroup(selected_outgroup_node)
        rooting = "outgroup"
    except:
        copy_of_gene_tree.set_outgroup(copy_of_gene_tree.get_midpoint_outgroup())
        rooting = "midpoint"
    
    rooted_gene_trees_w_support.append(copy_of_gene_tree)

# Now we record branch length and support values for gene tree.
cumulative_branch_length = 0
list_of_support_values = []
minimum_support_values = {"family_id": input_tree_file.split("/")[-1][:-3], "rooting": rooting,
                          "speci": [], "gene_id": [], "minimum_support": []}
for selected_gene_tree in rooted_gene_trees_w_support:
    for node in selected_gene_tree.traverse():
        cumulative_branch_length += node.dist
        if not node.is_leaf():
            list_of_support_values.append(node.support)
        else:
            minimum_support_values["speci"].append(node.name.split("__")[0])
            minimum_support_values["gene_id"].append(return_special_chars(node.name.split("__")[1]))
            minimum_support_values["minimum_support"].append(min([ancestor.support for ancestor in node.get_ancestors()]))

overview_table = DataFrame(minimum_support_values)
overview_table = overview_table.sort_values("minimum_support").drop_duplicates(["speci", "gene_id"])
overview_table["median_support"] = median(list_of_support_values)
overview_table["mean_support"] = mean(list_of_support_values)
overview_table["total_branch_length"] = cumulative_branch_length
overview_table["number_genes"] = len(gene_tree)
overview_table.to_csv(input_tree_file.replace("gene_trees", "gene_tree_metrics").replace("nw", "csv"), index = False)
