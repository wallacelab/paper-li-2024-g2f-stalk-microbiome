library(dplyr)
library(ggplot2)
library(gridExtra)
library(vegan) # for mantel test
library(qiime2R)
library(ggpubr)

setwd("D:/G2F_Rerun/G2F_data/G2F_data-main/")

# Function to read QZA files and convert to a data frame
read_qza_to_df <- function(path) {
  as.data.frame(as.matrix(read_qza(path)$data))
}

# Reading environmental factor distance matrices
Relative_Humidity_distance <- read_qza_to_df("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/Relative.Humidity...._distance_matrix.qza")
Potassium_distance <- read_qza_to_df("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/Potassium.ppm.K_distance_matrix.qza")
Solar_radiation_distance <- read_qza_to_df("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/Solar.Radiation..W.m2._distance_matrix.qza")
Temperature_distance <- read_qza_to_df("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/Temperature..C._distance_matrix.qza")

# Reading UniFrac distance matrices
location_weighted_unifrac <- read_qza_to_df("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/weighted_unifrac_distance_matrix.qza")
location_unweighted_unifrac <- read_qza_to_df("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-no-NCH1-group-by-location-from-rarefied_table-18000/unweighted_unifrac_distance_matrix.qza")


# Function to intersect and align distance matrices based on common location names
align_distance_matrices <- function(unifrac_distance, env_distance) {
  intersect_names <- intersect(rownames(unifrac_distance), rownames(env_distance))
  aligned_distance <- env_distance[intersect_names, intersect_names]
  return(aligned_distance)
}

# Aligning environmental factor distance matrices with UniFrac distance matrices
Relative_Humidity_distance_aligned <- align_distance_matrices(location_weighted_unifrac, Relative_Humidity_distance)
Potassium_distance_aligned <- align_distance_matrices(location_weighted_unifrac, Potassium_distance)
Temperature_distance_aligned <- align_distance_matrices(location_weighted_unifrac, Temperature_distance)
Solar_radiation_distance_aligned <- align_distance_matrices(location_weighted_unifrac, Solar_radiation_distance)



# Function to extract distance vectors from distance matrices
calc_distance_vector <- function(DM) {
  lower_tri <- DM[lower.tri(DM)]
  return(lower_tri)
}

# Function to perform linear modeling and summarize results
perform_lm <- function(dependent_vector, independent_vector) {
  lm_model <- lm(dependent_vector ~ independent_vector)
  return(summary(lm_model))
}

# Function to create a ggplot with regression line and annotations
create_plot <- function(x, y, x_label, y_label, p_value, r_squared) {
  ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(size = 3, color = "black") +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    xlab(x_label) + ylab(y_label) +
    theme_minimal() +
    annotate("text", x = Inf, y = Inf, label = sprintf("p-value: %.3f\nR²: %.3f", p_value, r_squared), hjust = 1.1, vjust = 1.1, size = 5)
}

# Function to perform Mantel test and extract results
perform_mantel <- function(dm1, dm2) {
  mantel_test <- mantel(dm1, dm2, method = "pearson", permutations = 999)
  return(list(p_value = mantel_test$signif, r_squared = mantel_test$statistic))
}

# Extract distance vectors
Relative_Humidity_distance_vector <- calc_distance_vector(Relative_Humidity_distance_aligned)
Potassium_distance_vector <- calc_distance_vector(Potassium_distance_aligned)
Solar_radiation_distance_vector <- calc_distance_vector(Solar_radiation_distance_aligned)
Temperature_distance_vector <- calc_distance_vector(Temperature_distance_aligned)

weighted_unifrac_distance_vector <- calc_distance_vector(location_weighted_unifrac)
unweighted_unifrac_distance_vector <- calc_distance_vector(location_unweighted_unifrac)

# Combine distance vectors into a dataframe
G2F_distance_data <- data.frame(
  Temperature_distance_vector,
  Solar_radiation_distance_vector,
  Relative_Humidity_distance_vector,
  Potassium_distance_vector,
  weighted_unifrac_distance_vector,
  unweighted_unifrac_distance_vector
)

# Linear modeling
lm_results <- list(
  weighted_RH = perform_lm(weighted_unifrac_distance_vector, Relative_Humidity_distance_vector),
  weighted_K = perform_lm(weighted_unifrac_distance_vector, Potassium_distance_vector)
  # Add more as needed
)

# Mantel tests
mantel_results <- list(
  weighted_K = perform_mantel(Potassium_distance_aligned,location_weighted_unifrac),
  weighted_T = perform_mantel(Temperature_distance_aligned,location_unweighted_unifrac),
  
  # Add more as needed
)

# Plots
plots_list <- list(
  K_weighted = create_plot(
    G2F_distance_data$Potassium_distance_vector,
    G2F_distance_data$weighted_unifrac_distance_vector,
    "Distance in Relative Humidity",
    "Weighted UniFrac Distance",
    mantel_results$weighted_K$p_value,
    mantel_results$weighted_K$r_squared
  ),
  T_weighted = create_plot(
    G2F_distance_data$Temperature_distance_vector,
    G2F_distance_data$weighted_unifrac_distance_vector,
    "Distance in Relative Humidity",
    "Weighted UniFrac Distance",
    mantel_results$weighted_T$p_value,
    mantel_results$weighted_T$r_squared
  )
  # Add more plots as needed
)

plots_list



# Arrange and save the plots
G2F_environment_weighted_plot <- ggarrange(plotlist = plots_list, nrow = 1, ncol = 2)
ggsave("./G2F_environment_weighted_plot.png", G2F_environment_weighted_plot, width = 20, height = 10, device = "png")

G2F_environment_weighted_plot
