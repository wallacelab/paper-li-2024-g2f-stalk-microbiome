library(ggplot2)
library(vegan)

setwd("D:/G2F_Rerun/G2F_data/G2F_data-main/")

# Function to read QZA file and return a matrix
read_qza_matrix <- function(path) {
  qza_data <- read_qza(path)
  return(as.matrix(qza_data$data))
}

# Function to perform Mantel test and generate plot if significant
perform_mantel_and_plot <- function(genotype, env_factor, env_distance, unifrac_distance, prefix) {
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
    
    ggsave(paste0("./", env_factor, "_", gsub("/", "X", genotype), "_", prefix, "_unifrac.png"), width = 10, height = 10)
  }
}

# Iterating through genotypes and environmental factors
for (genotype in unique(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$Corrected_pedigree)) {
  print(genotype)
  genotype_dir_name <- gsub("\\/", "X", genotype)
  directory_name <- paste0("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-"
                           ,genotype_dir_name, "-by-location-3000")
  
  print(genotype_dir_name)
  print(directory_name)
  
  weighted_path <- paste0(directory_name, "/weighted_unifrac_distance_matrix.qza")
  unweighted_path <- paste0(directory_name, "/unweighted_unifrac_distance_matrix.qza")
  
  genotype_location_weighted_unifrac <- read_qza_matrix(weighted_path)
  genotype_location_unweighted_unifrac <- read_qza_matrix(unweighted_path)
  
  for (env_factor in environment_factors) {
    print(env_factor)
    env_factor_path <- paste0("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/", env_factor, "_distance_matrix.qza",sep="")
    env_distance <- read_qza_matrix(env_factor_path)
    
    # Perform Mantel tests and plot for weighted UniFrac distances
    perform_mantel_and_plot(genotype, env_factor, env_distance, genotype_location_weighted_unifrac, "weighted")
    
    # Perform Mantel tests and plot for unweighted UniFrac distances
    perform_mantel_and_plot(genotype, env_factor, env_distance, genotype_location_unweighted_unifrac, "unweighted")
  }
}

# Check for significant associations
print(significant_association_weighted_unifrac_genotype_environment)
print(summary(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$Corrected_pedigree))




############################################################################################################

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
    
    plot_name <- paste0("./", env_factor, "_", gsub("/", "X", genotype), "_", prefix, "_unifrac.png")
    ggsave(plot_name, gg, width = 10, height = 10)
  }
}

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

# Iterating through genotypes and environmental factors
environment_factors <- c("X1.1.Soil.pH", "WDRF.Buffer.pH", "X1.1.S.Salts.mmho.cm", "Organic.Matter.LOI..", "Nitrate.N.ppm.N", "lbs.N.A",
                         "Potassium.ppm.K", "Sulfate.S.ppm.S", "Calcium.ppm.Ca", "Magnesium.ppm.Mg", "CEC.Sum.of.Cations.me.100g", "Temperature..C.", 
                         "Relative.Humidity....", "Rainfall..mm.", "Solar.Radiation..W.m2.")

for (genotype in unique(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$Corrected_pedigree)) {
  print(genotype)
  genotype_dir_name <- gsub("\\/", "X", genotype)
  directory_name <- paste0("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-", 
                           genotype_dir_name, "-by-location-3000")
  
  weighted_path <- paste0(directory_name, "/weighted_unifrac_distance_matrix.qza")
  unweighted_path <- paste0(directory_name, "/unweighted_unifrac_distance_matrix.qza")
  
  genotype_location_weighted_unifrac <- read_qza_matrix(weighted_path)
  genotype_location_unweighted_unifrac <- read_qza_matrix(unweighted_path)
  
  for (env_factor in environment_factors) {
    print(env_factor)
    env_factor_path <- paste0("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/", env_factor, "_distance_matrix.qza")
    env_distance <- read_qza_matrix(env_factor_path)
    
    # Perform Mantel tests and plot for weighted UniFrac distances
    perform_mantel_and_plot(genotype, env_factor, env_distance, genotype_location_weighted_unifrac, "weighted")
    
    # Perform Mantel tests and plot for unweighted UniFrac distances
    perform_mantel_and_plot(genotype, env_factor, env_distance, genotype_location_unweighted_unifrac, "unweighted")
  }
}
