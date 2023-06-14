# Script: Takes the output from OptRoot.linux and generates a series of files 
# that is suitable for input into Ranger-DTL.linux.
# Usage: python extract_optimal_roots.py optroot_input_file optroot_output_file 

import sys

optroot_input = sys.argv[1]
optroot_output = sys.argv[2]

# Reading of the species tree (first line) from the input file for OptRoot.
with open(optroot_input) as in_file:
    species_tree = in_file.readline()

# Will store all calculated possible rootings as a set.
rootings = set()

# Opening file that was output by the OptRoot program.
with open(optroot_output, "r") as in_file:
    for line in in_file:
        if line.startswith("("):
            rootings.add(line)

# For each rooting, creating a file that can be used as input to Ranger-DTL.
counter = 1
for rooting in rootings:
    with open("species_gene_tree_{:d}.nw".format(counter), "a+") as out_file:
        # Writing species tree.
        out_file.write(species_tree)
        # Writing selected rooting of gene tree.
        out_file.write(rooting)
    counter += 1