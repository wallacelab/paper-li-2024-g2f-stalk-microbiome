#! /usr/bin/Rscript

# Calculate GxE for beta diversity metrics

library(ggplot2)
library(phyloseq)
library(qiime2R)
library(rbiom)
library(tidyverse) 
library(vegan)

# Global variables
nperms = 999 # number of permutations for distance-based redundancy analysis
set.seed(1) # Hoping this works for anova.cca() permutations; documentation unclear, but tests indicate should work
sig_threshold=0.05 # Threshold for calling a model term significant

# Load necessary data
asvs = readRDS("3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds")
metadata = asvs %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")

# Calculate distance matrices from scratch
## Extract needed data
asv_matrix = asvs %>% otu_table() %>% as.matrix()
mytree = asvs %>% phy_tree
## Run distance calculations
distances=list()
distances[["bray"]] = rbiom::beta.div(asv_matrix, method="bray-curtis") %>% as.matrix()
distances[["jaccard"]] = rbiom::beta.div(asv_matrix, method="jaccard") %>% as.matrix()
distances[["weighted"]] = rbiom::beta.div(asv_matrix, method="unifrac", tree=mytree, weighted=TRUE) %>% as.matrix()
distances[["unweighted"]] = rbiom::beta.div(asv_matrix, method="unifrac", tree=mytree, weighted=FALSE) %>% as.matrix()


########################################################################


# Function to calculate beta heritability from distance matrices and metadata
beta_heritability_cal <- function(distance_matrix, input_metadata) {

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
    dplyr::select(-factor) %>%
    filter(test_var != "Residual")
  
  return(contribution_df)
}

#############################################################################

# Analyze beta diversity data
herits = lapply(distances, beta_heritability_cal, input_metadata=metadata)

# Construct contribution tables
tables = lapply(names(herits), function(metric){
  construct_contribution_table(herits[[metric]], type=metric)
}) %>%
  bind_rows() %>%
  mutate(category="Beta Diversity") %>%
  relocate(category)

# Plot the contributions of G, E, and GXE on alpha and beta diversity
plotdata = tables %>%
  mutate(alpha=ifelse(`Pr(>F)` < sig_threshold, yes=1, no=0.75))
notsig = subset(plotdata, `Pr(>F)` >= sig_threshold)
betaplot = ggplot(plotdata) +
  aes(x = Type, y = heritability, fill = test_var, alpha=alpha) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Variance Explained") +
  ggtitle("GxE for Beta Diversity") + 
  theme_minimal() +
  theme(legend.position = "bottom", axis.title.x = element_blank(), legend.title = element_blank()) +
  scale_alpha(guide='none') +
  geom_text(data=notsig, label="n.s.", alpha=1, position=position_stack(vjust=0.5))

# Save the plot
ggsave(betaplot, file="3_GxE/3c_beta_diversity_GxE.jgw.png", height = 3, width = 3, device = "png")

# Save the output data 
write.csv(tables, file="3_GxE/3c_beta_diversity_GxE.jgw.csv", row.names=FALSE)
