#! /usr/bin/Rscript

# Calculate GxE for beta diversity metrics

# TODO: DBRDA supposed to use distance matrix, but giving it first 2 PCs instead.
#    That seems wrong. Correct?
# TODO: Export a text table of values (with p-values); may need for later analyses
# TODO: Confirm that dbrda correctly matches distance matrix with values based on sample name

library(ggplot2)
library(phyloseq)
library(qiime2R)
library(tidyverse)
library(vegan)

# Global variables
nperms = 1000 # number of permutations for distance-based redundancy analysis

# Load necessary data
weighted_unifrac = read_qza("0_data_files/weighted_unifrac_distance_matrix.qza")
unweighted_unifrac = read_qza("0_data_files/unweighted_unifrac_distance_matrix.qza")
metadata = readRDS("3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds") %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")

########################################################################


# Function to calculate beta heritability from distance matrices and metadata
beta_heritability_cal <- function(input_qza, input_metadata) {
  # Extract distance matrix
  distance_matrix <- input_qza$data %>% as.matrix()
  
  # Subset to just the common taxa
  taxa = intersect(rownames(distance_matrix), input_metadata$SampleID)
  
  # Match the distance matrix with the provided metadata based on sample IDs
  input_metadata <- subset(input_metadata, SampleID %in% taxa)
  distance_matrix <- as.matrix(distance_matrix)[taxa, taxa]
  distance_matrix = as.dist(distance_matrix)
  
  # Run dbRDA and analyze
  dbRDA_result <- dbrda(distance_matrix ~ location + Corrected_pedigree + location:Corrected_pedigree, input_metadata)
  anova_result <- anova(dbRDA_result, by = "terms", permu = nperms)
  
  # Calculate heritability and create a dataframe for the results
  heritability_df <- heritability_calculation_beta(anova_result)
  
  return(heritability_df)
}

#########################################################################

# Function to calculate heritability based on ANOVA results
heritability_calculation_beta <- function(anova_result) {
  ssq_total <- sum(anova_result$SumOfSqs)
  heritability <- anova_result$SumOfSqs / ssq_total
  
  # Create a dataframe with heritability results
  heritability_df <- data.frame(
    Factor = rownames(anova_result),
    Heritability = heritability,
    stringsAsFactors = FALSE
  )
  
  return(heritability_df)
}

#######################################################################3

# Function to construct contribution tables from heritability data
construct_contribution_table <- function(heritability_df, type) {
  # Extract heritability values for Environment, Maize Genotype, and GXE
  environment_heritability <- heritability_df$Heritability[heritability_df$Factor == "location"]
  maize_genotype_heritability <- heritability_df$Heritability[heritability_df$Factor == "Corrected_pedigree"]
  gxe_heritability <- heritability_df$Heritability[heritability_df$Factor == "location:Corrected_pedigree"]
  
  # Construct a dataframe for contribution data
  contribution_df <- data.frame(
    test_var = c("Environment", "Maize Genotype", "GXE"),
    Value = c(environment_heritability, maize_genotype_heritability, gxe_heritability),
    Type = type,
    stringsAsFactors = FALSE
  )
  
  return(contribution_df)
}

#############################################################################

# Analyze weighted and unweighted UniFrac data
weighted_unifrac_heritability_df <- beta_heritability_cal(weighted_unifrac, metadata)
unweighted_unifrac_heritability_df <- beta_heritability_cal(unweighted_unifrac, metadata)

# Construct contribution tables
weighted_unifrac_contribution_table <- construct_contribution_table(weighted_unifrac_heritability_df, "Weighted")
unweighted_unifrac_contribution_table <- construct_contribution_table(unweighted_unifrac_heritability_df, "Unweighted")

# Combine weighted and unweighted UniFrac contribution tables
unifrac_contribution_table <- dplyr::bind_rows(weighted_unifrac_contribution_table, unweighted_unifrac_contribution_table)

# Prepare for combined alpha and beta diversity analysis
unifrac_contribution_table$category <- "Beta Diversity"
#colnames(alpha_result)[1] <- "Type"  # Ensure this aligns with the 'Type' column in the UniFrac tables

# Combine alpha and beta diversity contribution tables
#alpha_beta_diversity_table <- dplyr::bind_rows(unifrac_contribution_table, alpha_result)

# Plot the contributions of G, E, and GXE on alpha and beta diversity
betaplot = ggplot(data = unifrac_contribution_table, aes(x = Type, y = Value, fill = test_var)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Variance Explained") +
  ggtitle("GxE for Beta Diversity") + 
  theme_minimal() +
  theme(legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank())

# Save the plot
ggsave(betaplot, file="3_GxE/3c_beta_diversity_GxE.png", height = 3, width = 3, device = "png")
