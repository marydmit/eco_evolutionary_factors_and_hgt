# Script: Takes a gene tree and a species tree as an input, collapses species tree
# to clades containing gene of interest, returns corresponding gene and species tree
# structures.
# Usage: python prepare_input_for_ranger.py species_tree_file gene_tree_file

import sys

from ete3 import Tree
import pandas as pd

gene_tree_file = sys.argv[2]
gene_tree = Tree(gene_tree_file)
# Writing gene tree to file in structure-only format.
gene_tree.write(outfile = gene_tree_file + ".structure", format = 9)

# Recording all pairwise distances between genes in the tree.
pairwise_distances = {"gene_1": [], "gene_2": [], "distance": []}
gene_list = gene_tree.get_leaf_names()

for i in range(len(gene_list) - 1):
    for j in range(i+1, len(gene_list)):
        
        leaf_1 = gene_tree&gene_list[i]
        leaf_name_1 = gene_list[i].replace("semicolon", ";").replace("colon", ":").replace("dot", ".").replace("dash", "-")
        
        leaf_2 = gene_tree&gene_list[j]
        leaf_name_2 = gene_list[j].replace("semicolon", ";").replace("colon", ":").replace("dot", ".").replace("dash", "-")

        distance = round(leaf_1.get_distance(leaf_2), 3)

        pairwise_distances["gene_1"].extend([leaf_name_1, leaf_name_2])
        pairwise_distances["gene_2"].extend([leaf_name_2, leaf_name_1])
        pairwise_distances["distance"].extend([distance, distance])
pd.DataFrame(pairwise_distances).to_csv(gene_tree_file + ".distances", index = False)

# Determining whether there are any multifurcations in the gene tree that obscure
# horizontal gene transfer detection.
distances_gene_tree = gene_tree.get_cached_content(store_attr = "dist")
counter = 0
with open(gene_tree_file + ".multifurcation", "w") as out_file:
    for node in distances_gene_tree:
        if not node.is_leaf() and distances_gene_tree[node] == {0}:
            leaves_from_node = node.get_leaf_names()
            if len(leaves_from_node) > 2:
                for leaf in node.get_leaf_names():
                    out_file.write("{},{}\n".format(counter, leaf.replace("semicolon", ";").replace("colon", ":").replace("dot", ".").replace("dash", "-")))
                counter += 1

species_tree_file = sys.argv[1]
species_tree = Tree(species_tree_file, format = 0)
species_with_gene = set(leaf.split("_")[0] for leaf in gene_tree.get_leaf_names())
species_tree = species_tree.get_common_ancestor(list(species_with_gene))

# Traversing the species tree, naming internal nodes and assigning specI clusters that contain gene.
# In addition, we also calculate cumulative branch length of the species tree.
counter = 1
species_tree_dist_full = 0
species_tree_size_full = len(species_tree)
species_tree_max_dist = 0

for node in species_tree.traverse(strategy = "preorder"):
    species_tree_dist_full += node.dist  # branch length summation
    if node.is_leaf():
        # For each leaf, indicating whether specI cluster present in gene tree.
        if node.name in species_with_gene:
            node.add_feature("gene_tree", "yes")
        else:
            node.add_feature("gene_tree", "no")
        # In addition, obtaining the farthest leaf and recording maximum distance
        # between species in the tree.
        most_distant_node = node.get_farthest_node()[1]
        if most_distant_node > species_tree_max_dist:
            species_tree_max_dist = most_distant_node
    else:
        # For each internal node, naming it so we have it in the collapsed tree.
        node.name = "internal_" + str(counter)
        counter += 1

# Extracting specI cluster names and gene presence data in a convenient to parse format.
node_names = species_tree.get_cached_content(store_attr = "name")
node_gene_presence = species_tree.get_cached_content(store_attr = "gene_tree")

def collapsed_node(node):
    # If node is internal and all its descendants don't have the gene, collapse node.
    if len(node_names[node]) > 1 and node_gene_presence[node] == set(["no"]):
        return True
    # If we reached a leaf node, doesn't matter whether the gene is present or absent.
    elif len(node_names[node]) == 1:
        return True
    # All other internal nodes we look into (possible collapsing deeper in the tree).
    else:
        return False

# Collapsing the species tree before writing it to file in structure-only format.
collapsed_species_tree = Tree(species_tree.write(is_leaf_fn = collapsed_node))
collapsed_species_tree.write(outfile = species_tree_file + ".collapsed", format = 0)
collapsed_species_tree.write(outfile = species_tree_file + ".structure", format = 9)

# Exporting species tree metrics into separate file.
sp_tree_dict = {"size_full": [species_tree_size_full],
                "sum_dist_full": [round(species_tree_dist_full, 3)],
                "max_dist_full": [round(species_tree_max_dist, 3)],
                "size_coll": [len(collapsed_species_tree)],
                "sum_dist_coll": [round(sum([node.dist for node in collapsed_species_tree.traverse()]), 3)]}
                

sp_tree_df = pd.DataFrame(sp_tree_dict)
sp_tree_df.to_csv(path_or_buf = "species_tree_metrics.csv", index = False)