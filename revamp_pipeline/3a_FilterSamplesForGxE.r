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
filtered = metadata  %>%
  relocate(SampleID) %>%
  group_by(location, Corrected_pedigree) %>%  
  filter(n() >= min_reps_per_loc)

#keep location that has more than 10 samples 
filtered <- filtered %>% 
  group_by(location) %>% 
  filter(n() >= 10)

##keep pedigrees present in least 3 unique locations
filtered <- filtered %>% 
  group_by(Corrected_pedigree) %>% 
  filter(length(unique(location)) >= 3)

##########
# Plot diagnostics
##########

# Helper function to make plots; "tag" gets added to the title of each
make_plots = function(mydata, tag=""){
  
  # Number of locations for each genotype
  geno_counts = mydata %>%
    group_by(Corrected_pedigree) %>%
    summarize(count=length(unique(location)), .groups="drop")
  
  genoplot <- ggplot(geno_counts) +
    aes(x=Corrected_pedigree,y=count) + 
    geom_col(fill="darkgreen") +
    geom_abline(slope=0, intercept=min_locs_present, color="magenta", linetype="dashed", linewidth=1.5) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
    ggtitle(paste("Locations with each genotype", tag))
  
  # Number of genotypes in each location
  loc_counts = mydata %>%
    group_by(location) %>%
    summarize(count=length(unique(Corrected_pedigree)), .groups="drop")
  
  locplot <- ggplot(loc_counts) +
    aes(x=location,y=count) + 
    geom_col(fill="darkred") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
    ggtitle(paste("Unique genotypes per location", tag))
  
  # Overall heatmap of counts per genotype per location
  counts = mydata %>%
    group_by(location, Corrected_pedigree) %>%
    summarize(count=n(), .groups="drop")
  
  myheatmap <- ggplot(counts) +
    aes(x=location,y=Corrected_pedigree, fill=count) + 
    geom_tile() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
    ggtitle(paste("Counts of genotypes per location", tag))
  
  combined = grid.arrange(grobs=list(genoplot, locplot, myheatmap), nrow=1, widths=c(1,1,1.5))
  return(combined)

}

# Make plots
preplots = make_plots(metadata, "(pre)")
postplots = make_plots(filtered, "(post)")
allplots = grid.arrange(preplots, postplots, nrow=2)

# Save 
ggsave(allplots, file="3_GxE/3a_filter_checks.png", width=12, height=8)


############
# Write out filtered data
############

outdata = prune_samples(filtered$SampleID, ps.rarefied)
saveRDS(outdata, file="3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds")