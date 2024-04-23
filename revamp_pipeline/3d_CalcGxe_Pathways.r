#! /usr/bin/Rscript

# Look at GxE for inferred biochemical pathways

# Load libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(tidyverse)
library(phyloseq)

# TODO - There's a filtering of the pathways with an amount of 70. Is this correct? Presumably for prevalence?

# Global variables
pathway_file="0_data_files/path_abun_unstrat_descrip.tsv"

# Load metadata
metadata = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.phyloseq.rds") %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")

# Define functions for heritability calculations
heritability_calculation <- function(anova_result, sum_sq_index) {
  sapply(anova_result, function(res) {
    ssq_values <- res$`Sum Sq`
    ssq_total <- sum(ssq_values)
    ssq_values[sum_sq_index] / ssq_total
  })
}

# Read and preprocess pathway descriptions
pathway_description <- read.table(pathway_file, sep = "\t", header = TRUE, check.names = TRUE, row.names = 1)
colnames(pathway_description) <- gsub("\\.", "-", colnames(pathway_description))
pathway_description$description <- NULL  # Remove description to prevent console crash
pathway_description <- pathway_description[rowSums(pathway_description == 0) <= 70, ] #TODO - This correct?
pathways = rownames(pathway_description)

# Prepare pathway data
input_pathway_table <- t(pathway_description) %>% as.data.frame()
input_pathway_table$SampleID <- rownames(input_pathway_table)
input_pathway_table_metadata <- inner_join(input_pathway_table, metadata, by = "SampleID")

# Linear regression on each pathway
lm_results <- lapply(input_pathway_table_metadata[, pathways], function(x) lm(x ~ location + Corrected_pedigree + location:Corrected_pedigree, data = input_pathway_table_metadata))
anova_results <- lapply(lm_results, anova)

# Calculate heritabilities
location_heritabilities <- heritability_calculation(anova_results, 1)
pedigree_heritabilities <- heritability_calculation(anova_results, 2)
GxE_heritabilities <- heritability_calculation(anova_results, 3)

# Organize results into data frames
heritability_data <- data.frame(
  Location = location_heritabilities,
  Pedigree = pedigree_heritabilities,
  GxE = GxE_heritabilities,
  Taxa = pathways
  #Taxa = colnames(input_pathway_table_metadata[, 1:280]) # Fragile. Replaced with above
)

# Visualization of heritability results
plot_heritability <- function(heritability_data) {
  heritability_data %>%
    pivot_longer(-Taxa, names_to = "Type", values_to = "Heritability") %>%
    ggplot(aes(x = Type, y = Heritability, fill = Type)) +
    geom_violin(trim = TRUE) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "black")
}

# Generate plot
heritability_plot <- plot_heritability(heritability_data)
#print(heritability_plot)

# Save plot
ggsave("3_GxE/3d_MetaCyc_pathway_GXE.png", plot = heritability_plot, height = 10, width = 10, device = "png")

