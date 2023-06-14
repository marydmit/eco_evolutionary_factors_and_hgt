# Aim:  Compare rates of horizontal gene transfer in generalists and specialists.

library(plyr)
library(stringr)
library(ggplot2)
library(fitdistrplus)
library(RColorBrewer)
library(reshape2)

# ---- FUNCTIONS ----

extract_generalists_and_specialists <- function(abundance_table, environment) {
  
  selected_abundances <- abundance_table[abundance_table$pref_env == environment, ]
  
  # For the analysis, we take 200 most generalist and 200 most specialist species.
  specialists <- selected_abundances[order(selected_abundances$entr_gen), "otu_id"][1:200]
  generalists <- selected_abundances[order(selected_abundances$entr_gen, decreasing = TRUE), "otu_id"][1:200]
  
  return(list("specialists" = specialists, "generalists" = generalists))
}

determine_phylogenetic_distance_distribution <- function(otu_pairs_overview_table, phylogenetic_distance_bins) {
  # Function: takes table of OTU pairs and returns the target phylogenetic distribution in such a way that
  #           each environment will have the same number of OTU pairs to consider within each phylogenetic
  #           distance bin.
  otu_pairs_overview_table <- otu_pairs_overview_table[otu_pairs_overview_table$distance > phylogenetic_distance_bins[1] & otu_pairs_overview_table$distance <= phylogenetic_distance_bins[length(phylogenetic_distance_bins)], ]
  otu_pairs_overview_table$environment_combo <- paste(pmin(otu_pairs_overview_table$environment_1, otu_pairs_overview_table$environment_2), pmax(otu_pairs_overview_table$environment_1, otu_pairs_overview_table$environment_2))
  environment_combos <- unique(otu_pairs_overview_table$environment_combo)
  
  # We will subsample all environment groups to the one with the smallest number of observations in
  # selected phylogenetic distance bin.
  target_phylogenetic_distribution <- sapply(environment_combos, function(combo) {
    distances_per_combo <- otu_pairs_overview_table[otu_pairs_overview_table$environment_combo == combo, "distance"]
    return(hist(distances_per_combo, breaks = phylogenetic_distance_bins, plot = F)$counts)
  })
  target_phylogenetic_distribution <- apply(target_phylogenetic_distribution, 1, min)
  return(target_phylogenetic_distribution)
}

perform_phylogenetic_distance_correction_to_target <- function(otu_pairs_overview_table, phylogenetic_distance_bins, target_phylogenetic_distribution, environment_combos) {
  # Function: subsamples pairs in each environment combination in such a way that each combination
  # follows the same phylogenetic distance distribution.

  # Performing OTU pair subsampling.
  subsampled_overview_table <- lapply(environment_combos, function(combo) {
    original_table <- otu_pairs_overview_table[otu_pairs_overview_table$environment_combo == combo, ]
    sampled_table  <- lapply(seq(2, length(phylogenetic_distance_bins), 1), function(bin) {
      bin_min_distance <- phylogenetic_distance_bins[bin-1]
      bin_max_distance <- phylogenetic_distance_bins[bin]
      otu_pairs_within_bin <- original_table[(original_table$distance > bin_min_distance) & (original_table$distance <= bin_max_distance), ]
      if(nrow(otu_pairs_within_bin) == target_phylogenetic_distribution[bin-1]) {
        return(otu_pairs_within_bin)
      } else {
        sampled_pairs_within_bin <- otu_pairs_within_bin[sample(nrow(otu_pairs_within_bin), target_phylogenetic_distribution[bin-1]), ]
        return(sampled_pairs_within_bin)
      }})
    return(ldply(sampled_table, as.data.frame))})
  subsampled_overview_table <- ldply(subsampled_overview_table, as.data.frame)
  return(subsampled_overview_table)
}

summarize_transfers_per_environment_combo <- function(otu_pairs_overview_table) {
  
  summary_transfer <- tapply(otu_pairs_overview_table$transfer, otu_pairs_overview_table$environment_combo, function(x) sum(x)/length(x))
  summary_table <- data.frame(environment_1 = sapply(names(summary_transfer), function(x) strsplit(x, " ")[[1]][1]),
                              environment_2 = sapply(names(summary_transfer), function(x) strsplit(x, " ")[[1]][2]),
                              pairs_w_transfer = summary_transfer)
  
  summary_table$environment_1 <- factor(summary_table$environment_1)
  levels(summary_table$environment_1) <- sapply(levels(summary_table$environment_1), str_to_title)
  
  summary_table$environment_2 <- factor(summary_table$environment_2, levels = c("soil", "plant", "aquatic", "animal"))
  levels(summary_table$environment_2) <- sapply(levels(summary_table$environment_2), str_to_title)
  
  return(summary_table)
}

# ---- IMPORTING DATA ----

# Importing all OTU data with transfers, distances, and co-occurrence.
hgt_overview_table <- read.table("~/Documents/hgt_project/file_repository/OTU_pairs_with_phylogenetic_distance_genes_transferred_and_cooccurrence.csv",
                                 sep = ",", header = TRUE, stringsAsFactors = FALSE)
hgt_overview_table$transfer <- as.numeric(hgt_overview_table$genes_transferred != 0)

# Importing data on abundances (to select preferred environment).
otu_abundances <- read.table("~/Documents/hgt_project/file_repository/OTU_average_abundance_by_environment.csv",
                             sep = ",", header = TRUE, stringsAsFactors = F)
otu_abundances$pref_env <- apply(otu_abundances[ , -1], 1, function(x) names(x)[which(x == max(x))])
rownames(otu_abundances) <- otu_abundances$otu_id

# Assigning preferred environment based on relative abundance data.
hgt_overview_table$environment_1 <- otu_abundances[hgt_overview_table$otu97_1, "pref_env"]
hgt_overview_table$environment_2 <- otu_abundances[hgt_overview_table$otu97_2, "pref_env"]
# Generating list of environment combos to look at.
hgt_overview_table$environment_combo <- paste(pmin(hgt_overview_table$environment_1, hgt_overview_table$environment_2), pmax(hgt_overview_table$environment_1, hgt_overview_table$environment_2))
environment_combos <- unique(na.omit(hgt_overview_table)$environment_combo)

# Importing data on generalism and specialism.
otu_generalism_index <- read.table("~/Documents/hgt_project/generalism_specialism/generalism_stats_abundweighted.tsv", sep = "\t", header = TRUE,
                                   stringsAsFactors = FALSE, quote = "")
rownames(otu_generalism_index) <- otu_generalism_index$oid

# ---- ANALYSIS ----

# For each environment, determining generalists and specialists.
otu_abundances$entr_gen <- otu_generalism_index[sapply(otu_abundances$otu_id, function(x) paste0(substr(x, 1, 1), strsplit(x, ";")[[1]][4])), "entr_gen"]
otu_generalists_and_specialists <- list("animal" = extract_generalists_and_specialists(otu_abundances, "animal"),
                                        "aquatic" = extract_generalists_and_specialists(otu_abundances, "aquatic"),
                                        "plant" = extract_generalists_and_specialists(otu_abundances, "plant"),
                                        "soil" = extract_generalists_and_specialists(otu_abundances, "soil"))
list_of_generalists <- c(otu_generalists_and_specialists$animal$generalists, otu_generalists_and_specialists$aquatic$generalists, otu_generalists_and_specialists$plant$generalists, otu_generalists_and_specialists$soil$generalists)
list_of_specialists <- c(otu_generalists_and_specialists$animal$specialists, otu_generalists_and_specialists$aquatic$specialists, otu_generalists_and_specialists$plant$specialists, otu_generalists_and_specialists$soil$specialists)

# Determining phylogenetic distance distributions within all OTUs, generalists and specialists.
target_all <- determine_phylogenetic_distance_distribution(na.omit(hgt_overview_table), seq(0.1, 5.5, 0.2))

generalists_hgt_overview_table <- hgt_overview_table[hgt_overview_table$otu97_1 %in% list_of_generalists & hgt_overview_table$otu97_2 %in% list_of_generalists, ]
target_generalists <- determine_phylogenetic_distance_distribution(na.omit(generalists_hgt_overview_table), seq(0.1, 5.5, 0.2))

specialists_hgt_overview_table <- hgt_overview_table[hgt_overview_table$otu97_1 %in% list_of_specialists & hgt_overview_table$otu97_2 %in% list_of_specialists, ]
target_specialists <- determine_phylogenetic_distance_distribution(na.omit(specialists_hgt_overview_table), seq(0.1, 5.5, 0.2))

# We will subsample data from the other categories in such a way that all compared OTUs follow the same
# phylogenetic distance distribution.
target_distribution <- apply(data.frame("all" = target_all, "gen" = target_generalists, "sp" = target_specialists), 1, min)

# Bringing generalists, specialists and all species to the same phylogenetic distance distribution.
corrected_hgt_overview_same_target <- perform_phylogenetic_distance_correction_to_target(na.omit(hgt_overview_table), seq(0.1, 5.5, 0.2), target_distribution, environment_combos)
all_species_transfers_summarized_same_target <- summarize_transfers_per_environment_combo(corrected_hgt_overview_same_target)

corrected_generalists_hgt_overview_same_target <- perform_phylogenetic_distance_correction_to_target(na.omit(generalists_hgt_overview_table), seq(0.1, 5.5, 0.2), target_distribution, environment_combos)
generalists_transfers_summarized_same_target <- summarize_transfers_per_environment_combo(corrected_generalists_hgt_overview_same_target)

corrected_specialists_hgt_overview_same_target <- perform_phylogenetic_distance_correction_to_target(na.omit(specialists_hgt_overview_table), seq(0.1, 5.5, 0.2), target_distribution, environment_combos)
specialists_transfers_summarized_same_target <- summarize_transfers_per_environment_combo(corrected_specialists_hgt_overview_same_target)


# Plotting the results.
ggplot(data = all_species_transfers_summarized_same_target, aes(x = environment_1, y = environment_2, fill = pairs_w_transfer)) +
  geom_tile(color = "black") +
  geom_text(color = "black", size = 6, aes(label = format(round(pairs_w_transfer, 2), nsmall = 2))) +
  xlab("") +
  ylab("") + 
  ggtitle("All Species") +
  coord_fixed() + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_distiller(palette = "BrBG", limits = c(-0.165, 0.165), direction = 1) + # limits = c(0, 0.3)
  theme(legend.position = "bottom",
        text = element_text(size = 23),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 25,
                               barheight = 1,
                               title = "Fraction Species Pairs with Transfers",
                               title.position = "top",
                               title.hjust = 0.5))
ggplot(data = generalists_transfers_summarized_same_target, aes(x = environment_1, y = environment_2, fill = pairs_w_transfer)) +
  geom_tile(color = "black") +
  geom_text(color = "black", size = 6, aes(label = format(round(pairs_w_transfer, 2), nsmall = 2))) +
  xlab("") +
  ylab("") + 
  ggtitle("Generalists") +
  coord_fixed() + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_distiller(palette = "BrBG", limits = c(-0.165, 0.165), direction = 1) + # limits = c(0, 0.3)
  theme(legend.position = "bottom",
        text = element_text(size = 23),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 25,
                               barheight = 1,
                               title = "Fraction Species Pairs with Transfers",
                               title.position = "top",
                               title.hjust = 0.5))
ggplot(data = specialists_transfers_summarized_same_target, aes(x = environment_1, y = environment_2, fill = pairs_w_transfer)) +
  geom_tile(color = "black") +
  geom_text(color = "black", size = 6, aes(label = format(round(pairs_w_transfer, 2), nsmall = 2))) +
  xlab("") +
  ylab("") + 
  ggtitle("Specialists") +
  coord_fixed() + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_distiller(palette = "BrBG", limits = c(-0.165, 0.165), direction = 1) + # limits = c(0, 0.3)
  theme(legend.position = "bottom",
        text = element_text(size = 23),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 25,
                               barheight = 1,
                               title = "Fraction Species Pairs with Transfers",
                               title.position = "top",
                               title.hjust = 0.5))


# ---- BACKGROUND GENERATION ----

# Since we subsample very aggressively for all species to make it comparable to generalists and specialists.
# I have run 1'000 subsamplings for all species to generate a background distribution of transfer rates.
lapply(seq(1, 20), function(subset_number) {
random_otus_transfer_summary_table <- lapply(seq(1, 50), function(replicate) {
  
  corrected_hgt_overview_same_target <- perform_phylogenetic_distance_correction_to_target(na.omit(hgt_overview_table), seq(0.1, 5.5, 0.2), target_distribution, environment_combos)
  all_species_transfers_summarized_same_target <- summarize_transfers_per_environment_combo(corrected_hgt_overview_same_target)
  
  colnames(all_species_transfers_summarized_same_target) <- c("environment_1", "environment_2", paste0("frac_hgt_", (subset_number - 1) * 50 + replicate))
  return(all_species_transfers_summarized_same_target)
})

random_otus_transfer_summary_table <- join_all(random_otus_transfer_summary_table)
write.table(random_otus_transfer_summary_table, sprintf("overview_fraction_transfers_different_environments_all_otus_subsampled_%s.txt", subset_number),
            sep = "\t", row.names = FALSE, quote = FALSE)
})

list_of_tables <- list.files(pattern = "22_02_03_overview_fraction_transfers_different_environments_all_otus_subsampled")
background_estimation_table <- lapply(list_of_tables, function(x) read.table(x, header = TRUE, stringsAsFactors = FALSE))
background_estimation_table <- join_all(background_estimation_table)

# ---- STANDARD DEVIATION ----
all_species_sd <- apply(background_estimation_table[ , 3:ncol(background_estimation_table)], 2, sd)
generalist_sd <- sd(generalists_transfers_summarized_same_target$pairs_w_transfer)
specialist_sd <- sd(specialists_transfers_summarized_same_target$pairs_w_transfer)
fitted_distribution <- fitdist(all_species_sd * 10, "norm")

hist(all_species_sd, breaks = 100, xlim = c(0, 0.05))
abline(v = generalist_sd)
abline(v = specialist_sd)

generalist_score_sd <- (10*generalist_sd - fitted_distribution$estimate[1]) / fitted_distribution$estimate[2]
2 * pnorm(10*generalist_sd, mean = fitted_distribution$estimate[1], sd = fitted_distribution$estimate[2], lower.tail = TRUE)
# 1.334324e-25
specialist_score_sd <- (10*specialist_sd - fitted_distribution$estimate[1]) / fitted_distribution$estimate[2]
2 * pnorm(10*specialist_sd, mean = fitted_distribution$estimate[1], sd = fitted_distribution$estimate[2], lower.tail = FALSE)
# 2.938403e-299


# ---- RANGE ----
all_species_range <- apply(background_estimation_table[ , 3:ncol(background_estimation_table)], 2, function(x) max(x) - min(x))
generalist_range <- max(generalists_transfers_summarized_same_target$pairs_w_transfer) - min(generalists_transfers_summarized_same_target$pairs_w_transfer)
specialist_range <- max(specialists_transfers_summarized_same_target$pairs_w_transfer) - min(specialists_transfers_summarized_same_target$pairs_w_transfer)
fitted_distribution_range <- fitdist(all_species_range, "norm")

hist(all_species_range, breaks = 100, xlim = c(0.025, 0.15))
abline(v = generalist_range)
abline(v = specialist_range)

generalist_score_range <- (generalist_range - fitted_distribution_range$estimate[1]) / fitted_distribution_range$estimate[2]
2 * pnorm(generalist_range, mean = fitted_distribution_range$estimate[1], sd = fitted_distribution_range$estimate[2], lower.tail = TRUE)
# 5.295096e-11
specialist_score_range <- (specialist_range - fitted_distribution_range$estimate[1]) / fitted_distribution_range$estimate[2]
2 * pnorm(specialist_range, mean = fitted_distribution_range$estimate[1], sd = fitted_distribution_range$estimate[2], lower.tail = FALSE)
# 6.209193e-216


# ---- AVERAGE ----
all_species_mean <- apply(background_estimation_table[ , 3:ncol(background_estimation_table)], 2, mean)
generalist_mean <- mean(generalists_transfers_summarized_same_target$pairs_w_transfer)
specialist_mean <- mean(specialists_transfers_summarized_same_target$pairs_w_transfer)
fitted_distribution <- fitdist(all_species_mean * 10, "norm")

generalist_score_mean <- (10* generalist_mean - fitted_distribution$estimate[1]) / fitted_distribution$estimate[2]
2 * pnorm(10* generalist_mean, mean = fitted_distribution$estimate[1], sd = fitted_distribution$estimate[2], lower.tail = FALSE)
# 0.007630833
specialist_score_mean <- (10* specialist_mean - fitted_distribution$estimate[1]) / fitted_distribution$estimate[2]
2 * pnorm(10 * specialist_mean, mean = fitted_distribution$estimate[1], sd = fitted_distribution$estimate[2], lower.tail = FALSE)
# 1.050435e-05

p.adjust(c(1.334324e-25, 2.938403e-299, 5.295096e-11, 6.209193e-216, 0.007630833, 1.050435e-05), method = "fdr")
# 2.668648e-25 1.763042e-298  7.942644e-11 1.862758e-215  7.630833e-03  1.260522e-05
