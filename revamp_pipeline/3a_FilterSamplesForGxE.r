#! /usr/bin/Rscript

# Filter the input dataset to only those samples suitable for GxE analysis 
#  (Meaning, present in 2 reps at each location)
library(phyloseq)
library(tidyverse)
library(gridExtra)

# Filter parameters, in the order to be applied
min_reps_per_loc = 2 # How many reps a genotype needs per location to be kept
min_samples_per_loc = 10 # Minimum number of samples at a location to keep it (after above filtering)
min_locs_present = 3 # How many different locations a genotype need to be in to be kept (after above)


#############
# Load data
#############

ps.rarefied = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")

##############
# Filter for samples with min reps per location and min locationss
##############

# Extract and format metadata for manipulation
metadata=sample_data(ps.rarefied) %>%
  data.frame() %>%
  rownames_to_column("SampleID")  # Make sample ID a column and not just the row names

# Keep only genotypes where is replicated a sufficient number of times per location
replicated = metadata  %>%
  relocate(SampleID) %>%
  group_by(location, Corrected_pedigree) %>%  
  filter(n() >= min_reps_per_loc)

#keep location that has more than 10 samples 
replicated <- replicated %>% 
  group_by(location) %>% 
  filter(n() >= 10)

# TODO: THE BELOW IS INCORRECT. It's just looking at 3 pedigree instances total, not 3 per location
##keep location has at least 3 unique pedigree
replicated <- replicated %>% 
  group_by(Corrected_pedigree) %>% 
  filter(length(unique(location)) >= 3)




summary(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$location)
summary(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$Corrected_pedigree)

##check sample name 
G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$SampleID

#change . in sample name to - for later operation 
G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$SampleID <- gsub("\\.","-",G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$SampleID)

G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$SampleID

#G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$location <- as.character(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$location)

##check location and pedigree 
summary(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$location)
summary(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$Corrected_pedigree)

G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$location <- as.factor(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$location)
G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$Corrected_pedigree <- as.factor(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$Corrected_pedigree)


G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$location <- as.character(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$location)

##format to dataframe for string substitution
G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location <- as.data.frame(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location)

# #merge OH43/B37 and B37/OH43 -> Moved to earlier step (1a_)
# G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$Corrected_pedigree <- gsub(".*^OH43/B37","B37/OH43",G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$Corrected_pedigree)

#create a heatmap to visualize the distribution of pedigree across the location 
correted_pedigree_heatmap <-ggplot(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location, aes(x=location,y=Corrected_pedigree)) + 
  geom_tile()

correted_pedigree_heatmap 

