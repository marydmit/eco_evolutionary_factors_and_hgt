#!/usr/bin/env nextflow

params.dev = false
params.number_of_inputs = 150

// Software and scripts used in the pipeline.
count_speci = params.script_dir + "count_number_speci_clusters.sh"
remove_duplicates = params.script_dir + "remove_duplicate_artefacts.py"
mafft = "/mnt/mnemo4/maria/anaconda3/bin/mafft"
calculate_identity = params.script_dir + "calculate_sequence_identity.py"
fasttree = "/mnt/mnemo4/maria/software/FastTree"
prepare_input_ranger = params.script_dir + "prepare_input_for_ranger.py"
find_optimal_root = params.ranger_dir + "OptRoot.linux"
extract_optimal_roots = params.script_dir + "extract_optimal_roots.py"
reconcile_tree = params.ranger_dir + "Ranger-DTL.linux"
aggregate_ranger = params.ranger_dir + "AggregateRanger_recipient.linux"
aggregate_results = params.script_dir + "aggregate_ranger_results.py"
look_up_divergence = params.script_dir + "map_to_otus_and_calculate_divergence.py"

// Path to directory "pipeline_input", which contains FASTA files for each individual gene family.
gene_sequences = Channel
                  .fromPath("pipeline_input/*")
                  .map { file -> tuple(file.baseName, file) }
                  .take( params.dev ? params.number_of_inputs : -1 )

// A version of the species tree which contains only structure, no branch lengths or support values.
species_tree = "/mnt/mnemo4/maria/horizontal_gene_transfer/species_tree_progenomes_v2.2_no_chimera_reps_only_mapping_to_otus.nwk"

// Table containing mapping between specI clusters and OTUs.
map_specI_clusters_to_OTUs = "/mnt/mnemo4/maria/horizontal_gene_transfer/mapping_progenomes_v2.2_speci_clusters_to_mapref_v2.2b_OTUs.tsv"

// Tables containing taxonomic information for the specI clusters used in the analysis.
specI_cluster_taxonomy = "/mnt/mnemo4/maria/horizontal_gene_transfer/progenomes_v2_updated_ncbi_taxonomy.tab"
specI_cluster_gtdb = "/mnt/mnemo4/maria/horizontal_gene_transfer/mapping_progenomes_v2.2_genomes_to_GTDB_taxonomy_release95.tsv"

/*
 * Counting how many different specI clusters are present within the gene family.
 */
process count_speci_clusters {
    
    tag "$familyID"

    input:
    set familyID, "sequences.fna" from gene_sequences

    output:
    set stdout, familyID, "sequences.fna" into gene_sequences_speci

    script:
    """  
    #!/bin/bash
    
    $count_speci sequences.fna 
    """

}

/*
 * Checking gene families with at least five different specI clusters if they contain
 * multiple genomes assigned to the same specI cluster. These would be labelled as
 * duplication events by Ranger-DTL, when in reality they are not.
 */
process filter_speci_duplicates {
    
    tag "$familyID"

    publishDir "filtered_sequences", pattern: "*.fna", mode: "copy"
    publishDir "duplication_logs", pattern: "*.log", mode: "copy"

    input:
    set numberSpecIs, familyID, "sequences.fna" from gene_sequences_speci

    output:
    set stdout, familyID, "${familyID}.fna", "${familyID}.log" into gene_sequences_filtered

    when:
    numberSpecIs.toInteger() >= params.minimum_speci

    script:
    """
    #!/bin/bash

    python $remove_duplicates sequences.fna ${familyID}.fna ${familyID}.log

    # Calculating what's the final number of sequences in family.
    grep -c ">" ${familyID}.fna
    """

}

/*
 * Running multiple sequence alignment for gene families with at least five specI clusters.
 */
process align_gene_sequences {
    
    tag "$familyID"

    publishDir "alignments", pattern: "*.fna.aln", mode: "copy"

    input:
    set clusterSize, familyID, "sequences.fna", "filter.log" from gene_sequences_filtered
    
    output:
    set familyID, "${familyID}.fna.aln" into alignments

    when:
    clusterSize.toInteger() >= params.min_cluster_size && clusterSize.toInteger() <= params.max_cluster_size

    script:
    """
    #!/bin/bash
    
    # Multiple sequence alignment.
    $mafft --auto sequences.fna >${familyID}.fna.aln
    """

}

/*
 * For downstream analyses on HGT and gene divergence, we generate a percentage identity matrix.
 */
process calculate_sequence_identity {

    tag "$familyID"

    publishDir "gene_identities", pattern: "*.csv", mode: "copy"

    input:
    set familyID, "alignment.fna" from alignments

    output:
    set familyID, "alignment.fna", "${familyID}.csv" into alignments_with_identities

    script:
    """
    #!/bin/bash

    # Running script that will generate a percentage identity matrix.
    python $calculate_identity alignment.fna ${familyID}.csv
    """
}

/*
 * Generating a maximum likelihood tree based on provided alignments.
 */
process generate_gene_tree {
    
    tag "$familyID"

    publishDir "gene_trees", mode: "copy"

    input:
    set familyID, "alignment.fna", "identity_matrix.csv" from alignments_with_identities

    output:
    set familyID, "${familyID}.nw" into gene_trees

    script:
    """
    #!/bin/bash

    # Renaming the tree file to contain family id.
    $fasttree -gtr -nt <alignment.fna >${familyID}.nw
    """

}

/*
 * Prior to running Ranger-DTL, we need to modify the format of the species and gene trees.
 * For gene tree, removing branch lengths and branch support values.
 * For species tree, retaining only clades containing our gene of interest and collapsing 
 * clades that don't contain the gene of interest.
 */
process prepare_tree_for_ranger_input {
    
    tag "$familyID"

    publishDir "species_trees", pattern: "species_tree.nw.collapsed", mode: "copy", saveAs: {filename -> "${familyID}.nw"}
    publishDir "species_sizes", pattern: "species_tree_metrics.csv", mode: "copy", saveAs: {filename -> "${familyID}.csv"}
    publishDir "multifurcations", pattern: "gene_tree.nw.multifurcation", mode: "copy", saveAs: {filename -> "${familyID}.csv"}
    publishDir "gene_distances", pattern: "gene_tree.nw.distances", mode: "copy", saveAs: {filename -> "${familyID}.csv"}
    
    input:
    set familyID, "gene_tree.nw" from gene_trees
    path "species_tree.nw" from species_tree

    output:
    set familyID, "gene_tree.nw.structure", "species_tree.nw.structure" into tree_structures
    set familyID, "species_tree.nw.collapsed", "species_tree_metrics.csv" into species_tree_metrics
    set familyID, "gene_tree.nw.multifurcation" into multifurcation_flags
    set familyID, "gene_tree.nw.distances" into gene_distances

    script:
    """
    #!/bin/bash
    python $prepare_input_ranger species_tree.nw gene_tree.nw

    # Not done by ETE3 by default, so adding white space at the end of the tree files.
    sed -i -e '\$a\\' gene_tree.nw.structure
    sed -i -e '\$a\\' species_tree.nw.structure
    sed -i -e '\$a\\' species_tree.nw.collapsed
    """

}

/*
 * Using the species tree to find the optimal root for gene tree before reconciliation.
 */
process find_gene_tree_optimal_root {
        
    tag "$familyID"

    publishDir "optimal_roots", pattern: "*.optroots", mode: "copy"

    input:
    set familyID, "gene_tree.nw", "species_tree.nw" from tree_structures
 
    output:
    set familyID, "${familyID}.optroot.input", "${familyID}.optroots" into optimal_roots

    script:
    """
    #!/bin/bash

    #  Format required for input into Ranger-DTL: two lines, first line: species tree, second line: gene tree.
    cat species_tree.nw gene_tree.nw >${familyID}.optroot.input
    
    $find_optimal_root -i ${familyID}.optroot.input -o ${familyID}.optroots
    """
    
}

/*
 * Extracting all optimal roots detected and writing them in separate files as
 * preliminary step to tree reconciliation.
 */
process extract_rooted_gene_trees {
    
    tag "$familyID"

    input:
    set familyID, "optroot_input.txt", "optroot_output.txt" from optimal_roots

    output:
    set stdout, familyID, "species_gene_tree_*.nw" into ranger_input
    
    """
    #!/bin/bash
    
    python $extract_optimal_roots optroot_input.txt optroot_output.txt

    optimal_roots=\$(grep "The total number of optimal rootings is" optroot_output.txt | awk -F": " '{ print \$2 }')
    echo \$optimal_roots
    """

}

/*
 * Running Ranger-DTL to reconcile trees for each optimal rooting, totaling
 * 500 reconciliations. If the tree has more than 50 roots (leaving us with 
 * less than 10 reconciliations per root), we don't proceed further.
 */
process run_tree_reconciliations {
    
    tag "$familyID"
    
    scratch true
    
    input:
    set rootNumber, familyID, "species_gene_tree.nw*" from ranger_input

    output:
    set familyID, "${familyID}_root_*_1", "${familyID}_root_*.txt" into ranger_results
    
    when:
    rootNumber.toInteger() < 51

    script: 
    """
    #!/bin/bash
   
    # First, determining the number of reconciliations based on optimal roots.
    reconciliations=\$(( (500+${rootNumber}-1)/${rootNumber} ))

    root=1
    # For each optimal root that was found, we carry out multiple reconciliations.
    for file in species_gene_tree.nw*
      do for i in \$(seq \$reconciliations)
        do
          $reconcile_tree -i \$file -o ${familyID}_root_\${root}_\${i}
        done
        # Aggregating all the results generated for this particular root.
        $aggregate_ranger ${familyID}_root_\${root}_ >${familyID}_root_\${root}.txt
        # Incrementing counter for the next root.
        (( root ++ ))
      done   
    """

}

/*
 * Aggregating results from Ranger-DTL into transfer matrices.
 */
process aggregate_reconciliation_results {
        
    tag "$familyID"

    publishDir "pairs_wo_transfers", pattern: "${familyID}_no_transfers.csv", mode: "copy", saveAs: {filename -> "${familyID}.csv"}
    publishDir "transfer_mappings", pattern: "${familyID}_mappings.csv", mode: "copy", saveAs: {filename -> "${familyID}.csv"}
    publishDir "transfer_recipients", pattern: "${familyID}_recipients.csv", mode: "copy", saveAs: {filename -> "${familyID}.csv"}    
    
    input:
    set familyID, "${familyID}*", "${familyID}*.txt" from ranger_results
    
    output:
    set familyID, "${familyID}.csv" into aggregated_results
    set familyID, "${familyID}_no_transfers.csv", "${familyID}_mappings.csv", "${familyID}_recipients.csv" into aggregated_additional_info
    script:
    """
    #!/bin/bash

    python $aggregate_results ${familyID} .
    """

}

/*
 * Mapping specI clusters participating in the transfer to OTUs and noting how divergent
 * are these specI clusters to prioritize transfers from distantly-related species.
 */
process add_otus_and_divergence {

    tag "$familyID"

    publishDir "transfer_matrices", mode: "copy"

    input:
    set familyID, "${familyID}.csv" from aggregated_results
    path "otu_mapping.tsv" from map_specI_clusters_to_OTUs
    path "progenomes_taxonomy.tsv" from specI_cluster_taxonomy
    path "gtdb_taxonomy.tsv" from specI_cluster_gtdb

    output:
    set familyID, "${familyID}.csv" into results_with_mapping

    script:
    """
    #!/bin/bash

    python $look_up_divergence ${familyID}.csv otu_mapping.tsv progenomes_taxonomy.tsv gtdb_taxonomy.tsv
    """
}
