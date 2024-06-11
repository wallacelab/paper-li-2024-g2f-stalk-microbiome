#! /usr/bin/Rscript

# Calculate GxE for beta diversity metrics

library(ggplot2)
library(phyloseq)
library(qiime2R)
library(tidyverse)
library(vegan)

# Global variables
nperms = 999 # number of permutations for distance-based redundancy analysis
set.seed(1) # Hoping this works for anova.cca() permutations; documentation unclear, but tests indicate should work

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
  taxamatch = match(taxa, input_metadata$SampleID)
  input_metadata <- input_metadata[taxamatch,]
  distance_matrix <- as.matrix(distance_matrix)[taxa, taxa]
  distance_matrix = as.dist(distance_matrix)
  
  # Confirm that are in the same order so dbRDA works
  if(any(input_metadata$SampleID != rownames(as.matrix(distance_matrix)))){
    stop("FAIL: Metadata and distance matrix samples do not match! Unable to run dbRDA.")
  }else{
    cat("PASS: All metadata and distance matrix names match\n")
  }
  
  # Run dbRDA and analyze
  dbRDA_result <- dbrda(distance_matrix ~ location + Corrected_pedigree + location:Corrected_pedigree, input_metadata)
  anova_result <- anova(dbRDA_result, by = "terms", permutations = nperms)
  
  # Calculate heritability and create a dataframe for the results
  heritability_df <- data.frame(anova_result, check.names=FALSE) %>%
    mutate(heritability = SumOfSqs / sum(SumOfSqs))
  
  return(heritability_df)
}

#######################################################################3

# Function to construct contribution tables from heritability data
construct_contribution_table <- function(heritability_df, type) {

  contribution_df = heritability_df %>%
    rownames_to_column("factor") %>%
    mutate(Type = type) %>%
    mutate(test_var = case_when(factor == "location" ~ "Environment",
                                factor == "Corrected_pedigree" ~ "Maize Genotype",
                                factor == "location:Corrected_pedigree" ~ "GXE",
                                factor == "Residual" ~ "Residual",
                                TRUE ~ "UNKNOWN") ) %>%
    relocate(test_var, Type) %>%
    select(-factor) %>%
    filter(test_var != "Residual")
  
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

# Plot the contributions of G, E, and GXE on alpha and beta diversity
betaplot = ggplot(data = unifrac_contribution_table, aes(x = Type, y = heritability, fill = test_var)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Variance Explained") +
  ggtitle("GxE for Beta Diversity") + 
  theme_minimal() +
  theme(legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank())

# Save the plot
ggsave(betaplot, file="3_GxE/3c_beta_diversity_GxE.png", height = 3, width = 3, device = "png")

# Save the output data 
write.csv(unifrac_contribution_table, file="3_GxE/3c_beta_diversity_GxE.csv", row.names=FALSE)
