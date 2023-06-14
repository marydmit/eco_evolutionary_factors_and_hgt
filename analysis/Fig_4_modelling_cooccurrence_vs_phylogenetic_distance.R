# Aim:  Model co-occurrence versus phylogenetic distance to correct fo phylogenetic signal.
#       Correct for phylogenetic signal and determine relationship between co-occurrence and HGT.

# Required packages.
library(ggplot2)
library(nlstools)
library(plyr)

# ---- FUNCTIONS USED IN ANALYSIS ----

exponential_decay_model_cooccurrence_vs_phylogeny <- function(input_table, sqrt_transformation = FALSE) {
  # Function: For selected OTU, generates model of co-occurrence predicted based
  #           on phylogenetic distance using exponential decay formula. Uses residuals
  #           from this model to calculate Spearman correlation between co-occurrence
  #           and genes transferred.
  # Input:    input_table - data frame, contains data on OTU partners with prevalence, overlapping
  #           samples and phylogenetic distance.
  #           sqrt_transformation - boolean, whether co-occurrence data should be transformed for
  #           modelling (outliers less likely to skew model).
  # Output:   result - data frame, RSE - random squared error of model, N0, lambda - fitted
  #           parameters for model, rho - correlation between co-occurrence and genes transferred.
  input_table <- na.omit(input_table)
  input_table <- input_table[input_table$minimum_overlap >= 20, ]  # for better estimates, each species should be present in at least 20 samples
  input_table <- input_table[input_table$distance <= 2.5, ]  # setting cut-off on distances (only within-kingdom transfers)
  
  if(nrow(input_table) < 30){  # modeling is not performed if too few data points per OTU
    result <- rep("too few data", 6)
    names(result) <- c("RSE", "N0", "lambda", "rho_pre", "rho_post", "bg_rho")
    return(result)
  }
  
  # Setting initial parameters for fitting. These are based on the mode of the whole tested population.
  if(sqrt_transformation) {
    input_table$fraction_overlap <- sqrt(input_table$fraction_overlap)
    N0 = 0.25
    lambda = 0.45
  } else {
    N0 = 0.05
    lambda = 1
  }
  
  result <- tryCatch({
    # Fitting parameters to exponential decay model.
    model <- nls(fraction_overlap ~ N0 * exp(-lambda * distance), data = input_table, start = list(N0 = N0, lambda = lambda), control = list(maxiter = 100))
    
    # Calculating correlation between genes transferred and non-corrected sample overlap.
    pre_cor_coef <- cor(input_table$fraction_overlap, input_table$genes_transferred, method = "spearman")
    
    # Calculating correlation between genes transferred and residuals from sample overlap.
    post_cor_coef <- cor(residuals(model), input_table$genes_transferred, method = "spearman")
    
    # Randomizing genes transferred to get estimate for background rho.
    bg_cor_coef <- cor(residuals(model), sample(input_table$genes_transferred), method = "spearman")
    
    # Returning data from model and resulting correlation.
    c(sqrt(deviance(model)/(nrow(input_table)-2)), coef(model), pre_cor_coef, post_cor_coef, bg_cor_coef)
    
  }, error = function(e){
    # If modeling failed, reporting this.
    rep("modeling failed", 6)
  })
  
  names(result) <- c("RSE", "N0", "lambda", "rho_pre", "rho_post", "bg_rho")
  
  return(result)
}

power_law_model_cooccurrence_vs_phylogeny <- function(input_table, sqrt_transformation = FALSE) {
  # Function: For selected OTU, generates model of co-occurrence predicted based
  #           on phylogenetic distance using power law formula. Uses residuals
  #           from this model to calculate Spearman correlation between co-occurrence
  #           and genes transferred.
  # Input:    input_table - data frame, contains data on OTU partners with prevalence, overlapping
  #           samples and phylogenetic distance.
  #           sqrt_transformation - boolean, whether co-occurrence data should be transformed for
  #           modelling (outliers less likely to skew model).
  # Output:   result - data frame, RSE - random squared error of model, k, a - fitted
  #           parameters for model, rho - correlation between co-occurrence and genes transferred.
  input_table <- na.omit(input_table)
  input_table <- input_table[input_table$minimum_overlap >= 20, ]  # for better estimates, each species should be present in at least 
  input_table <- input_table[input_table$distance <= 2.5, ]  # setting cut-off on distances (only within-kingdom transfers) 
  
  if(nrow(input_table) < 30){  # modeling is not performed if too few data points
    result <- rep("too few data", 6)
    names(result) <- c("RSE", "k", "a", "rho_pre", "rho_post", "bg_rho")
    return(result)
  }
  
  # Setting initial parameters for fitting. These are based on the mode of the whole tested population.
  if(sqrt_transformation) {
    input_table$fraction_overlap <- sqrt(input_table$fraction_overlap)
    k = 0.1
    a = -0.3
  } else {
    k = 0.15
    a = -0.6
  }
  
  result <- tryCatch({
    # Fitting parameters to power law model.
    model <- nls(fraction_overlap ~ k * distance^a, data = input_table, start = list(k = k, a = a), control = list(maxiter = 100))
    
    # Calculating correlation between genes transferred and non-corrected sample overlap.
    pre_cor_coef <- cor(input_table$fraction_overlap, input_table$genes_transferred, method = "spearman")
    
    # Calculating correlation between genes transferred and residuals from sample overlap.
    post_cor_coef <- cor(residuals(model), input_table$genes_transferred, method = "spearman")
    
    # Randomizing genes transferred to get estimate for background rho.
    bg_cor_coef <- cor(residuals(model), sample(input_table$genes_transferred), method = "spearman")
    
    # Returning data from model and resulting correlation.
    c(sqrt(deviance(model)/(nrow(input_table)-2)), coef(model), pre_cor_coef, post_cor_coef, bg_cor_coef)
    
  }, error = function(e){
    # If modeling failed, reporting this.
    rep("modeling failed", 6)
  })
  
  names(result) <- c("RSE", "k", "a", "rho_pre", "rho_post", "bg_rho")
  
  return(result)
}

# ---- ANALYSIS: CALCULATING CORRELATIONS POST-CORRECTION FOR PHYLOGENY ----

hgt_overview_df <- read.table("OTU_pairs_with_phylogenetic_distance_genes_transferred_and_cooccurrence.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Only considering OTU pairs that transferred at least one gene.
hgt_overview_df <- hgt_overview_df[hgt_overview_df$genes_transferred > 0, ]

# Including transfers in both directions for modeling with 1 OTU selected for reference.
hgt_overview_df_2 <- hgt_overview_df[ , c("otu97_2", "otu97_1", "distance", "genes_transferred",
                                          "fraction_overlap", "intersection", "presence_2",
                                          "presence_1", "minimum_overlap")]
colnames(hgt_overview_df_2) <- colnames(hgt_overview_df)
hgt_overview_df <- rbind(hgt_overview_df, hgt_overview_df_2)

list_of_otus <- names(table(hgt_overview_df$otu97_1)[table(hgt_overview_df$otu97_1) >= 30])
# for analysis, only using OTUs with at least 30 different data points.

# Running analysis on selected OTUs (exponential decay).
result_exponential_decay <- sapply(list_of_otus, function(x) {
  selected_data <- hgt_overview_df[hgt_overview_df$otu97_1 == x, ]
  return(exponential_decay_model_cooccurrence_vs_phylogeny(selected_data, TRUE))
})
# Table transformations for better readability.
result_exponential_decay <- t(result_exponential_decay)
result_exponential_decay <- as.data.frame(result_exponential_decay)
result_exponential_decay$otu <- row.names(result_exponential_decay)
result_exponential_decay <- result_exponential_decay[ , c(7, 1, 2, 3, 4, 5, 6)]
write.table(result_exponential_decay, file = "modelling_results_exponential_decay.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# For plotting, need to remove "too few data" and "modeling failed" statements.
result_exponential_decay[result_exponential_decay == "too few data"] <- NA
result_exponential_decay[result_exponential_decay == "modeling failed"] <- NA
# Full dataset, OTUs for which modeling is possible: 3769, correlations for 3755

# Running analysis on selected OTUs (power law).
result_power_law <- sapply(list_of_otus, function(x) {
  selected_data <- hgt_overview_df[hgt_overview_df$otu97_1 == x, ]
  return(power_law_model_cooccurrence_vs_phylogeny(selected_data, TRUE))
})
# Table transformations for better readability.
result_power_law <- t(result_power_law)
result_power_law <- as.data.frame(result_power_law)
result_power_law$otu <- row.names(result_power_law)
result_power_law <- result_power_law[ , c(7, 1, 2, 3, 4, 5, 6)]
write.table(result_power_law, file = "modelling_results_power_law.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# For plotting, need to remove "too few data" and "modeling failed" statements.
result_power_law[result_power_law == "too few data"] <- NA
result_power_law[result_power_law == "modeling failed"] <- NA
# Full dataset, OTUs for which modeling is possible: 3772, correlations for 3758


# ---- PLOTTING OVERALL DISTRIBUTION ----

# We're basing observations and 3755 (exponential decay) and 3758 (power law) OTUs for the full dataset.
exponential_decay_plot_df <- data.frame(category = c(rep("Background", nrow(result_exponential_decay)),
                                                rep("Post-Correction", nrow(result_exponential_decay)),
                                                rep("Pre-Correction", nrow(result_exponential_decay))),
                                   value = as.numeric(c(result_exponential_decay$bg_rho, result_exponential_decay$rho_post,
                                             result_exponential_decay$rho_pre)))

# Plotting distribution for exponential decay modelling.
ggplot(data = exponential_decay_plot_df, aes(x = value, fill = category)) +
  geom_histogram(position = "identity", alpha = 0.6, col = "black", bins = 55) +
  scale_fill_manual(values = c("grey75", "#d95f0e", "#fff7bc")) +
  xlab("Correlation(Fraction Samples Shared, Genes Transferred)") +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Number of OTUs") +
  ggtitle("Correction Based on Exponential Decay Equation") +
  coord_fixed(1.72/575, xlim = c(-0.86, 0.86), ylim = c(0, 575)) +
  theme_classic() +
  theme(text = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.ticks.length = unit(.3, "cm"))



# Same procedure for power law modelling.
power_law_df <- data.frame(category = c(rep("Background", nrow(result_power_law)),
                                                rep("Post-Correction", nrow(result_power_law)),
                                                rep("Pre-Correction", nrow(result_power_law))),
                                   value = as.numeric(c(result_power_law$bg_rho, result_power_law$rho_post,
                                             result_power_law$rho_pre)))

ggplot(data = power_law_df, aes(x = value, fill = category)) +
  geom_histogram(position = "identity", alpha = 0.6, col = "black", bins = 55) +
  scale_fill_manual(values = c("grey75", "#d95f0e", "#fff7bc")) +
  xlab("Correlation(Fraction Samples Shared, Genes Transferred)") +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Number of OTUs") +
  ggtitle("Correction Based on Power Law Equation") +
  coord_fixed(1.72/575, xlim = c(-0.86, 0.86), ylim = c(0, 575)) +  # full dataset
  theme_classic() +
  theme(text = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.ticks.length = unit(.3, "cm"))
dev.off()


# ---- PLOTTING AN INDIVIDUAL EXAMPLE FOR MANUSCRIPT ----

# Generating example based on the complete dataset, regardless of gene distance.
selected_otu = "B16S;90_31;96_1396;97_1644"  # Bifidobacterium pseudolongum
data_selected_otu <- hgt_overview_df[hgt_overview_df$otu97_1 == selected_otu, ]
data_selected_otu <- na.omit(data_selected_otu)
data_selected_otu <- data_selected_otu[data_selected_otu$minimum_overlap >= 20, ]
data_selected_otu$sqrt_fraction_overlap <- sqrt(data_selected_otu$fraction_overlap)

pwl_model <- nls(sqrt_fraction_overlap ~ k * distance^a, data = data_selected_otu, start = list(k = 0.1, a = -0.3), control = list(maxiter = 100))
data_selected_otu$pwl_residuals <- residuals(pwl_model)

# Alternatively, can use the exponential decay model:
# exp_model <- nls(sqrt_fraction_overlap ~ N0 * exp(-lambda * distance), data = data_selected_otu, start = list(N0 = 0.25, lambda = 0.45), control = list(maxiter = 100))
# data_selected_otu$exp_residuals <- residuals(exp_model)

pdf("individual_examples/example_for_manuscript_Bifidobacterium_pseudolongum.pdf", width = 4, height = 4)
ggplot(data = data_selected_otu, aes(x = distance, y = sqrt_fraction_overlap)) +
  geom_point(alpha = 0.7, col = "grey50") +
  annotate("text", x = 1, y = 0.7, label = expression(rho == -0.79), size = 6, col = "black", parse = TRUE, hjust = 0) +
  xlab(expression("Distance to OTU"[2]*" on 16S rRNA Gene Tree")) +
  ylab(expression(frac("Shared Samples with OTU"[2], "min(Samples"[OTU[1]]*", Samples"[OTU[2]]*")"))) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.65), breaks = c(0, 0.5, 1, 1.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.9656176)) +
  theme_classic() +
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        aspect.ratio = 1)

ggplot(data = data_selected_otu, aes(x = distance, y = sqrt_fraction_overlap)) +
  geom_line(data = data.frame(dist = (seq(4, 160, by = 1)/100), frac = coef(pwl_model)[1] * (seq(4, 160, by = 1)/100) ^ coef(pwl_model)[2]), aes(x = dist, y = frac), col = "#d95f0e", size = 1) +
  geom_point(alpha = 0.7, col = "grey50") +
  annotate("text", x = 0.075, y = 0.9, label = expression(f(x) == 0.26*x^-0.41), size = 6, col = "#d95f0e", parse = TRUE, hjust = 0) +
  xlab(expression("Distance to OTU"[2]*" on 16S rRNA Gene Tree")) +
  ylab(expression(frac("Shared Samples with OTU"[2], "min(Samples"[OTU[1]]*", Samples"[OTU[2]]*")"))) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.65), breaks = c(0, 0.5, 1, 1.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.9656176)) +
  theme_classic() +
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        aspect.ratio = 1)

ggplot(data = data_selected_otu, aes(x = distance, y = sqrt_fraction_overlap)) +
  annotate("segment", x = data_selected_otu$distance, xend = data_selected_otu$distance, y = data_selected_otu$sqrt_fraction_overlap, yend = fitted(pwl_model), col = "#d95f0e") +
  geom_line(data = data.frame(dist = (seq(4, 160, by = 1)/100), frac = coef(pwl_model)[1] * (seq(4, 160, by = 1)/100) ^ coef(pwl_model)[2]), aes(x = dist, y = frac), col = "black", size = 1) +
  geom_point(alpha = 0.7, col = "grey50") +
  annotate("text", x = 0.075, y = 0.9, label = expression(f(x) == 0.26*x^-0.41), size = 6, col = "black", parse = TRUE, hjust = 0) +
  xlab(expression("Distance to OTU"[2]*" on 16S rRNA Gene Tree")) +
  ylab(expression(frac("Shared Samples with OTU"[2], "min(Samples"[OTU[1]]*", Samples"[OTU[2]]*")"))) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.65), breaks = c(0, 0.5, 1, 1.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.9656176)) +
  theme_classic() +
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        aspect.ratio = 1)

ggplot(data = data_selected_otu, aes(x = pwl_residuals, y = genes_transferred)) +
  geom_point(alpha = 0.7, col = "grey50") +
  annotate("text", x = 0.1, y = 40, label = expression(rho == 0.18), size = 6, col = "black", parse = TRUE, hjust = 0) +
  xlab(expression("Residuals on Shared Samples with OTU"[2])) +
  ylab("Genes Transferred") +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.31, 0.31)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) +
  theme_classic() +
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        aspect.ratio = 1)

ggplot(data = data_selected_otu, aes(x = sqrt_fraction_overlap, y = genes_transferred)) +
  geom_point(alpha = 0.7, col = "grey50") +
  annotate("text", x = 0.5, y = 40, label = expression(rho == 0.73), size = 6, col = "black", parse = TRUE, hjust = 0) +
  xlab(expression(frac("Shared Samples with OTU"[2], "min(Samples"[OTU[1]]*", Samples"[OTU[2]]*")"))) +
  ylab("Genes Transferred") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.9)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) +
  theme_classic() +
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        aspect.ratio = 1)

ggplot(data = data_selected_otu, aes(x = distance, y = genes_transferred)) +
  geom_point(alpha = 0.7, col = "grey50") +
  annotate("text", x = 1, y = 40,  label = expression(rho == -0.73), size = 6, col = "black", parse = TRUE, hjust = 0) +
  xlab(expression("Distance to OTU"[2]*" on 16S rRNA Gene Tree")) +
  ylab("Genes Transferred") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.65), breaks = c(0, 0.5, 1, 1.5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65)) +
  theme_classic() +
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(.3, "cm"),
        aspect.ratio = 1)
dev.off()
