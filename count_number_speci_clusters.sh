#!/bin/bash
# Usage: ./count_number_speci_clusters.sh $fasta_file

# Counting the number of different species in gene cluster to see whether we should
# proceed with multiple sequence alignment and building a phylogenetic tree.
numberSpecI=$(awk -F "__" '/>/{print substr($1, 2)}' $1 | sort -u | wc -l)

echo $numberSpecI
