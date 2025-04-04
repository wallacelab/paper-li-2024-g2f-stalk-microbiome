#! /usr/bin/Rscript

# Run associations between environmental parameters and beta diversity by genotype/location

library(ggplot2)
library(vegan)
library(qiime2R)
library(tidyverse)
library(phyloseq)
#library(ggpubr)
library(gridExtra)

# TODO: Jason concerned that Environment Unifrac matrices didn't have same number of locations
#   (=dimentions) when first ran. Is that an error?
# TODO: Jason concerned that Mantel test generating a bunch of "permutations < minperm", implying
#   that the number of possible permutations is small.Too small a dataset to work with?
# TODO: Too many images generated. Redo so just 1 per environmental factor

# Load metadata (not sure this exactly the same as originally specified, which was 
#   G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location. Not sure what
#   "filtered selected location" means)
metadata = readRDS("3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds") %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")

# Environmental factors to iterate over
environment_factors <- c("X1.1.Soil.pH", "WDRF.Buffer.pH", "X1.1.S.Salts.mmho.cm", "Organic.Matter.LOI..", "Nitrate.N.ppm.N", "lbs.N.A",
                         "Potassium.ppm.K", "Sulfate.S.ppm.S", "Calcium.ppm.Ca", "Magnesium.ppm.Mg", "CEC.Sum.of.Cations.me.100g", "Temperature..C.", 
                         "Relative.Humidity....", "Rainfall..mm.", "Solar.Radiation..W.m2.")

# Function to read QZA file and return a matrix
read_qza_matrix <- function(path) {
  qza_data <- read_qza(path)
  return(as.matrix(qza_data$data))
}

# Function to extract distance vectors from distance matrices
calc_distance_vector <- function(DM) {
  lower_tri <- DM[lower.tri(DM)]
  return(lower_tri)
}

# Function to perform Mantel test and generate plot if significant
perform_mantel_and_plot <- function(genotype, env_factor, env_distance, unifrac_distance, prefix) {
  # Add this to make sure matrices match
  locations = intersect(rownames(env_distance), rownames(unifrac_distance))
  env_distance = env_distance[locations, locations]
  unifrac_distance = unifrac_distance[locations, locations]
  
  mantel_result <- mantel(env_distance, unifrac_distance, method = "pearson", permutations = 999)
  
  if (mantel_result$signif <= 0.05 && !is.na(mantel_result$signif)) {
    env_vector <- calc_distance_vector(env_distance)
    unifrac_vector <- calc_distance_vector(unifrac_distance)
    data_to_plot <- data.frame(env_vector, unifrac_vector)
    
    ggplot(data_to_plot, aes(x = env_vector, y = unifrac_vector)) +
      geom_point(size = 3, color = "black") +
      geom_smooth(method = "lm", color = "blue", se = FALSE) +
      xlab(env_factor) +
      ylab(paste(prefix, "unifrac distance", sep = " ")) +
      theme_minimal() +
      annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, 
               label = sprintf("p-value = %.3f\nR² = %.3f", mantel_result$signif, mantel_result$statistic), size = 4)
    
    ggsave(paste0("5_Associations/5c_", env_factor, "_", gsub("/", "X", genotype), "_", prefix, "_unifrac.png"), width = 10, height = 10)
  }
}


# Function to align matrices based on intersecting location names and perform Mantel test
perform_mantel_and_plot <- function(genotype, env_factor, env_distance, unifrac_distance, prefix) {
  # Find intersecting location names
  location_names <- rownames(unifrac_distance)
  factor_location_names <- rownames(env_distance)
  intersect_names <- intersect(location_names, factor_location_names)
  
  # Subset matrices based on intersecting names
  env_distance_aligned <- env_distance[intersect_names, intersect_names]
  unifrac_distance_aligned <- unifrac_distance[intersect_names, intersect_names]
  
  # Perform Mantel test
  mantel_result <- mantel(env_distance_aligned, unifrac_distance_aligned, method = "pearson", permutations = 999)
  
  # Plot if the association is significant
  if (mantel_result$signif <= 0.05 && !is.na(mantel_result$signif)) {
    env_vector <- calc_distance_vector(env_distance_aligned)
    unifrac_vector <- calc_distance_vector(unifrac_distance_aligned)
    data_to_plot <- data.frame(env_vector, unifrac_vector)
    
    gg <- ggplot(data_to_plot, aes(x = env_vector, y = unifrac_vector)) +
      geom_point(size = 3, color = "black") +
      geom_smooth(method = "lm", color = "blue", se = FALSE) +
      xlab(env_factor) +
      ylab(paste(prefix, "UniFrac Distance", sep = " ")) +
      theme_minimal() +
      annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
               label = sprintf("p-value = %.3f\nR² = %.3f", mantel_result$signif, mantel_result$statistic), size = 4)
    
    plot_name <- paste0("5_Associations/5c_", env_factor, "_", gsub("/", "X", genotype), "_", prefix, "_unifrac.png")
    #ggsave(plot_name, gg, width = 10, height = 10)
    return(gg)
  }
  else{
    label=paste(env_factor,"\nnot assoicated with\n",genotype)
    nullplot = ggplot() + 
             annotate("text", x = 4, y = 25, size=8, label = label, color="red") + 
             theme_void()
    return(nullplot)
  }
}

# Set up lists to hold plots
weighted_plots = list()
unweighted_plots = list()
for (env_factor in environment_factors) {
  weighted_plots[[env_factor]] = list()
  unweighted_plots[[env_factor]] = list()
}

# Iterating through genotypes and environmental factors
for (genotype in unique(metadata$Corrected_pedigree)) {
  print(genotype)
  genotype_dir_name <- gsub("\\/", "X", genotype)
  directory_name <- paste0("0_data_files/distance_matrices_for_beta_associations_by_genotype/core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-"
                           ,genotype_dir_name, "-by-location-3000")
  
  weighted_path <- paste0(directory_name, "/weighted_unifrac_distance_matrix.qza")
  unweighted_path <- paste0(directory_name, "/unweighted_unifrac_distance_matrix.qza")
  
  genotype_location_weighted_unifrac <- read_qza_matrix(weighted_path)
  genotype_location_unweighted_unifrac <- read_qza_matrix(unweighted_path)
  
  for (env_factor in environment_factors) {

    print(env_factor)
    env_factor_path <- paste0("0_data_files/distance_matrices_for_beta_associations/", env_factor, "_distance_matrix.qza")
    env_distance <- read_qza_matrix(env_factor_path)
    
    # Perform Mantel tests and plot for weighted UniFrac distances
    weightplot = perform_mantel_and_plot(genotype, env_factor, env_distance, genotype_location_weighted_unifrac, "weighted")
    weighted_plots[[env_factor]][[genotype]] = weightplot
    
    # Perform Mantel tests and plot for unweighted UniFrac distances
    unweightplot = perform_mantel_and_plot(genotype, env_factor, env_distance, genotype_location_unweighted_unifrac, "unweighted")
    unweighted_plots[[env_factor]][[genotype]] = unweightplot
  }
}

# Output graphics
for (env_factor in environment_factors) {
  weightouts = grid.arrange(grobs=weighted_plots[[env_factor]])
  ggsave(weightouts, file=paste("5_Associations/5c_", env_factor,".weighted.png", sep=""),
         width=20, height=20)
  
  unweightouts = grid.arrange(grobs=unweighted_plots[[env_factor]])
  ggsave(unweightouts, file=paste("5_Associations/5c_", env_factor,".unweighted.png", sep=""),
         width=20, height=20)
}
