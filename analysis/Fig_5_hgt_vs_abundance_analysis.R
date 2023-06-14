# Aim:  Determine fraction of pairs participating in HGT based on abundance.
#       Correct for phylogenetic relationship using distances from 16S rRNA gene tree.

library(ggplot2)
library(dplyr)
library(gridExtra)
library(tools)

# ---- FUNCTIONS ----

generate_ab_vs_transfers_table <- function(overview_hgt_for_analysis, abundance_table, environment = "animal") {
  # We retain only OTUs that prefer this particular environment.
  abundance_table <- abundance_table[abundance_table$pref_env == environment, c("otu_id", environment)]
  
  # Retaining only information from OTUs that both are in this specific environment.
  result_table <- overview_hgt_for_analysis[(overview_hgt_for_analysis$otu97_1 %in% abundance_table$otu_id) & (overview_hgt_for_analysis$otu97_2 %in% abundance_table$otu_id), ]

  # Adding OTU relative abundance information.
  result_table$otu97_1_relab <- abundance_table[result_table$otu97_1, environment]
  result_table$otu97_2_relab <- abundance_table[result_table$otu97_2, environment]

  return(result_table[ , c("otu97_1", "otu97_2", "distance", "genes_transferred", "otu97_1_relab", "otu97_2_relab")])
}

determine_otu_quantiles <- function(environment, number_quantiles, otu_abundance_table) {
  otu_abundance_table <- otu_abundance_table[otu_abundance_table$pref_env == environment, c("otu_id", environment)]
  colnames(otu_abundance_table) <- c("otu_id", "relab")
  
  # Calculating relative abundance values at which we set thresholds.
  ab_quantiles <- sapply(seq(0, 1, length.out = number_quantiles + 1), function(x) {
    return(quantile(log10(otu_abundance_table$relab), x))})[-1]
  
  otu_quantiles <- sapply(otu_abundance_table$relab, function(x) {
    names(ab_quantiles[ab_quantiles >= log10(x)][1])})
  names(otu_quantiles) <- otu_abundance_table$otu_id
  return(otu_quantiles)
}

plot_abundance_vs_genes_transferred <- function(env_table_full, environment_quantiles, environment, distance_bin_size, min_distance = 0, max_distance = 2.5) {
  
  # For simplification, we only look at the two extremes for OTUs.
  environment_quantiles[environment_quantiles == "100%"] <- "High"
  environment_quantiles[environment_quantiles == "20%"] <- "Low"
  environment_quantiles <- environment_quantiles[environment_quantiles %in%  c("High", "Low")]
  
  # Assigning quantiles to OTU table generated for environment and dropping OTUs from other environments.
  env_table_full$quantile_1 <- environment_quantiles[env_table_full$otu97_1]
  env_table_full$quantile_2 <- environment_quantiles[env_table_full$otu97_2]
  env_table_filt <- na.omit(env_table_full)
  env_table_filt <- env_table_filt[env_table_filt$distance >= min_distance & env_table_filt$distance < (max_distance + distance_bin_size), ]
  env_table_filt$quantile_combo <- paste(env_table_filt$quantile_1, env_table_filt$quantile_2, sep = "-")
  env_table_filt$quantile_combo[env_table_filt$quantile_combo == "Low-High"] <- "High-Low"
  
  env_table_filt$distance_bin <- min_distance + ((((env_table_filt$distance - min_distance) * 1000) %/% (distance_bin_size * 1000)) * distance_bin_size)
  print(table(env_table_filt$distance_bin))
  
  summary_env <- env_table_filt %>% group_by(quantile_combo, distance_bin)  %>%
    summarize(fraction_transferred = sum(genes_transferred != 0) / length(genes_transferred),
              observations = length(genes_transferred))
  summary_env$uncertainty <- sqrt((summary_env$fraction_transferred * (1 - summary_env$fraction_transferred)) / summary_env$observations)
  
  summary_env <- summary_env[summary_env$quantile_combo %in% c("High-High", "High-Low", "Low-Low"), ]
  if(environment == "animal") {
    palette <- c("#AD1F20", "#D62728", "#EB9393") # "#fb6a4a", "#fcbba1"
  } else if(environment == "aquatic") {
    palette <- c("#15507A", "#1F77B4", "#6FB5E6") # "#6baed6", "#c6dbef"
  } else if(environment == "plant") {
    palette <- c("#217821", "#2CA02C", "#5FD35F") # "#74c476", "#c7e9c0"
  } else if(environment == "soil") {
    palette <- c("#CC5F00", "#FF7F0E", "#FFAB61") # "#fe9929", "#fee391"
  }
  
  plot_1 <- ggplot(summary_env, aes(x = distance_bin, y = fraction_transferred)) +
    geom_ribbon(aes(ymin = fraction_transferred - uncertainty, ymax = fraction_transferred + uncertainty, fill = quantile_combo), alpha = 0.3) +
    geom_line(aes(col = quantile_combo), size = 0.5) +
    geom_point(aes(col = quantile_combo), size = 1) +
    scale_color_manual(values = palette, name = NULL) +
    scale_fill_manual(values = palette, name = NULL) +
    xlab("Phylogenetic Distance Between OTU Pair") +
    ylab("Fraction OTU Pairs with Transfer") +
    ggtitle(toTitleCase(environment)) +
    scale_x_continuous(limits = c(0, max_distance + 0.05), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_fixed(ratio = max_distance + 0.05) +
    theme_classic() +
    theme(legend.position = "bottom",
          text = element_text(size = 16),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length = unit(.3, "cm"))
  return(list("data" = summary_env, "plot" = plot_1))}

plot_abundance_distribution <- function(env_table_full, environment) {
  otus_w_abundances_1 <- env_table_full[ , c("otu97_1", "otu97_1_relab")]
  colnames(otus_w_abundances_1) <- c("otu_id", "abundance")
  
  otus_w_abundances_2 <- env_table_full[ , c("otu97_2", "otu97_2_relab")]
  colnames(otus_w_abundances_2) <- c("otu_id", "abundance")
  
  list_all_otus <- rbind(otus_w_abundances_1, otus_w_abundances_2)
  list_all_otus <- distinct(list_all_otus)
  list_all_otus$abundance <- log10(list_all_otus$abundance) 
  
  low_abundance <- quantile(list_all_otus$abundance, 0.2)
  high_abundance <- quantile(list_all_otus$abundance, 0.8)
  list_all_otus$category <- sapply(list_all_otus$abundance, function(x) {
    ifelse(x <= low_abundance, "Low", ifelse(x >= high_abundance, "High", "Medium"))})

  if(environment == "animal") {
      sel_color <- "#D62728"
    } else if(environment == "aquatic") {
      sel_color <- "#1F77B4"
    } else if(environment == "plant") {
      sel_color <- "#2CA02C"
    } else if(environment == "soil") {
      sel_color <- "#FF7F0E"
    }
  
  plot_1 <- ggplot(data = list_all_otus, aes(x = abundance)) +
    geom_histogram(breaks = sapply(seq(-91, -7), function(x) x/7), fill = sel_color, color = "white", size = 0.1, alpha = 0.25) +
    geom_vline(xintercept = low_abundance, color = sel_color, size = 1) +
    geom_vline(xintercept = high_abundance, color = sel_color, size = 1) +
    # scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-13, -1), breaks = c(-12, -10, -8, -6, -4, -2), expand = c(0, 0)) +
    xlab("log10(Average Sample Abundance)") +
    ylab("Number of OTUs") +
    ggtitle(toTitleCase(environment)) +
    theme_classic() +
    theme(text = element_text(size = 16),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length = unit(.3, "cm"),
          aspect.ratio = 1)

 return(plot_1)
}

# ---- IMPORTING DATA FOR ANALYSIS ----

hgt_overview_table <- read.table("OTU_pairs_with_phylogenetic_distance_genes_transferred_and_cooccurrence.csv",
                                 sep = ",", header = T, stringsAsFactors = F)

# Determining preferred environment based on highest relative abundance.
sample_average_ab <- read.table("OTU_average_abundance_by_environment.csv", sep = ",",
                                header = T, stringsAsFactors = F)
sample_average_ab$pref_env <- apply(sample_average_ab[ , -1], 1, function(x) names(x)[which(x == max(x))])  # preferred environment is one where abundance is highest
row.names(sample_average_ab) <- sample_average_ab$otu_id

# ---- ANALYSIS ABUNDANCE vs FRACTION OF PAIRS WITH TRANSFER ----

# Extracting relative abundances in animal for OTUs that prefer animal-associated environments.
animal_otu_table <- generate_ab_vs_transfers_table(hgt_overview_table, sample_average_ab, "animal")
# Determining whether an OTU is present at high or low abundance when compared to overall OTU
# abundance distribution in the animal environment.
quantiles_animal_otus <- determine_otu_quantiles("animal", 5, sample_average_ab)
# Generating plots based on data.
animal_plot <- plot_abundance_vs_genes_transferred(animal_otu_table, quantiles_animal_otus, "animal", 0.2, 0.1, 1.5)
animal_abundance_dist <- plot_abundance_distribution(animal_otu_table, "animal")

# Same procedure for water-associated, plant-associated and soil-associated OTUs.
aquatic_otu_table <- generate_ab_vs_transfers_table(hgt_overview_table, sample_average_ab, "aquatic")
quantiles_aquatic_otus <- determine_otu_quantiles("aquatic", 5, sample_average_ab)
aquatic_plot <- plot_abundance_vs_genes_transferred(aquatic_otu_table, quantiles_aquatic_otus, "aquatic", 0.2, 0.1, 1.5)
aquatic_abundance_dist <- plot_abundance_distribution(aquatic_otu_table, "aquatic")

plant_otu_table <- generate_ab_vs_transfers_table(hgt_overview_table, sample_average_ab, "plant")
quantiles_plant_otus <- determine_otu_quantiles("plant", 5, sample_average_ab)
plant_plot <- plot_abundance_vs_genes_transferred(plant_otu_table, quantiles_plant_otus, "plant", 0.2, 0.1, 1.5)
plant_abundance_dist <- plot_abundance_distribution(plant_otu_table, "plant")

soil_otu_table <- generate_ab_vs_transfers_table(hgt_overview_table, sample_average_ab, "soil")
quantiles_soil_otus <- determine_otu_quantiles("soil", 5, sample_average_ab)
soil_plot <- plot_abundance_vs_genes_transferred(soil_otu_table, quantiles_soil_otus, "soil", 0.2, 0.1, 1.5)
soil_abundance_dist <- plot_abundance_distribution(soil_otu_table, "soil")

# Statistically comparing transfers in High-High, High-Low and Low-Low abundance OTU pairs.
wilcoxon_tests <- c(wilcox.test(unlist(animal_plot$data[animal_plot$data$quantile_combo == "High-High", "fraction_transferred"]),unlist(animal_plot$data[animal_plot$data$quantile_combo == "High-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(animal_plot$data[animal_plot$data$quantile_combo == "High-Low", "fraction_transferred"]), unlist(animal_plot$data[animal_plot$data$quantile_combo == "Low-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(animal_plot$data[animal_plot$data$quantile_combo == "High-High", "fraction_transferred"]), unlist(animal_plot$data[animal_plot$data$quantile_combo == "Low-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(aquatic_plot$data[aquatic_plot$data$quantile_combo == "High-High", "fraction_transferred"]), unlist(aquatic_plot$data[aquatic_plot$data$quantile_combo == "High-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(aquatic_plot$data[aquatic_plot$data$quantile_combo == "High-Low", "fraction_transferred"]), unlist(aquatic_plot$data[aquatic_plot$data$quantile_combo == "Low-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(aquatic_plot$data[aquatic_plot$data$quantile_combo == "High-High", "fraction_transferred"]), unlist(aquatic_plot$data[aquatic_plot$data$quantile_combo == "Low-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(plant_plot$data[plant_plot$data$quantile_combo == "High-High", "fraction_transferred"]), unlist(plant_plot$data[plant_plot$data$quantile_combo == "High-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(plant_plot$data[plant_plot$data$quantile_combo == "High-Low", "fraction_transferred"]), unlist(plant_plot$data[plant_plot$data$quantile_combo == "Low-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(plant_plot$data[plant_plot$data$quantile_combo == "High-High", "fraction_transferred"]), unlist(plant_plot$data[plant_plot$data$quantile_combo == "Low-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(soil_plot$data[soil_plot$data$quantile_combo == "High-High", "fraction_transferred"]), unlist(soil_plot$data[soil_plot$data$quantile_combo == "High-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(soil_plot$data[soil_plot$data$quantile_combo == "High-Low", "fraction_transferred"]), unlist(soil_plot$data[soil_plot$data$quantile_combo == "Low-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value,
                    wilcox.test(unlist(soil_plot$data[soil_plot$data$quantile_combo == "High-High", "fraction_transferred"]), unlist(soil_plot$data[soil_plot$data$quantile_combo == "Low-Low", "fraction_transferred"]), paired = TRUE, alternative = "greater")$p.value)
wilcoxon_tests
# [1] 0.00390625 0.00390625 0.00390625 0.00390625 0.07421875 0.02734375 0.00390625 0.00390625 0.00390625 0.01124714 0.01801584 0.01124714  Pre-correction
wilcoxon_tests <- p.adjust(wilcoxon_tests, method = "fdr")
wilcoxon_tests
# [1] 0.006696429 0.006696429 0.006696429 0.006696429 0.074218750 0.029829545 0.006696429 0.006696429 0.006696429 0.014996181 0.021619012 0.014996181 Correction with FDR

# Plotting transfer rates at different phylogenetic distances for high-abundance and low-abundance OTUs.
grid.arrange(animal_plot$plot, aquatic_plot$plot, plant_plot$plot, soil_plot$plot, nrow = 1, ncol = 4)

# Plotting overall abundance distributions.
grid.arrange(animal_abundance_dist, aquatic_abundance_dist, plant_abundance_dist, soil_abundance_dist, nrow = 1, ncol = 4)
