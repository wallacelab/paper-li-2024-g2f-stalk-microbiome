#! /usr/bin/Rscript

# Identify the core microbiome of the samples

library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(phyloseq)
library(qiime2R)

# TODO - Jason disagrees with this choice of subsample. Should use the entire dataset for 
#        a location, not just the duplicated pedigrees
# TODO - How look at Unique vs "Min 10 locations" subcores?

# Global variables
min_prevalence = 0.6 # Fraction of samples a taxon has to be in to be "core"
ranks=c("Phylum", "Class", "Order", "Family", "Genus", "Species") # Taxonomic ranks to test at

# Load data 

#mydata = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")
asvs = readRDS("3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds")
metadata = asvs %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")

# Precompute collapsed taxon profiles (takes a while)
collapsed_taxa = lapply(ranks, function(myrank){
  tax_glom(asvs, taxrank = myrank)
})
names(collapsed_taxa) = ranks

###############
# Helper functions
###############

# Function to construct a name key
get_namekey = function(mydata){
  namekey = tax_table(mydata) %>%
    as.data.frame() %>%
    apply(MARGIN=1, FUN=paste, collapse=";")
  namekey = data.frame(Taxaname = namekey) %>%
    rownames_to_column("taxon")
  return(namekey)
}

# Preparse the ASV data down to just locations and genotypes
core_microbiome_preparse <- function(tax_level, mydata){

  # Grab collapsed data & name key at given taxonomic level
  collapsed = collapsed_taxa[[tax_level]]
  namekey = get_namekey(collapsed)
  mytaxa = rownames(tax_table(collapsed))
  
  # Calculate prevalence at each location
  # Start by reformatting data into tidy data with location, taxon, and count
  collapsed_counts = otu_table(collapsed) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("SampleID") %>%
    left_join(metadata, by="SampleID") %>%
    relocate(SampleID, Corrected_pedigree, location) %>%
    select(SampleID, Corrected_pedigree, location, all_of(mytaxa)) %>%
    pivot_longer(cols=all_of(mytaxa), names_to="taxon", values_to="count") %>%
    mutate(present = count >0) %>% # Add a true/false column for if taxon present
    left_join(namekey, by="taxon")
  
}

# Plot heatmap
plot_heatmap = function(mycore, myrank, xval, outfile){
  # Format taxonomy
  taxonomy = str_split_fixed(mycore$Taxaname, ';', n=Inf)
  colnames(taxonomy) = rank_names(asvs)
  mycore$plot_taxon = taxonomy[,myrank]
  
  # Set which variable is on X axis
  mycore$xval = mycore[,xval] %>% unlist()
  
  # Make plot
  myplot = ggplot(mycore) +
    aes(x=xval, y=plot_taxon) +
    geom_tile() + 
    theme(
      axis.text = element_text(size=15, face="bold"),
      axis.text.x = element_text(angle=90),
      axis.title = element_text(size=15,face="bold")
    ) +
    labs(y = myrank, x=xval)
  
  # Save plot
  ggsave(myplot, file=outfile, height=8, width=12)
}

##########
# Preparse data
##########

preparse = lapply(ranks, FUN=core_microbiome_preparse, mydata = asvs)
names(preparse) = ranks

###############
# Prevalence by location
###############

location_cores = lapply(ranks, function(myrank){
  prevalence = preparse[[myrank]] %>%
    group_by(location, taxon, Taxaname) %>%
    summarize(prevalence = sum(present) / n(), .groups="drop") %>%
    filter(prevalence >= prevalence_level)
  return(prevalence)
})
names(location_cores) = ranks

# Plot and save
locplots = lapply(ranks, function(myrank){
  mycore = location_cores[[myrank]]
  outfile=paste("4_CoreMicrobiome/4a_location_",myrank,"_core_microbiome_heatmap.jgw.png", sep="")
  plot_heatmap(mycore, myrank, xval="location", outfile=outfile)
})

# # TODO - Still to look at - Subset to just those taxa that are "core" in a single location
# unique_location_genus_core_microbiome_list <- location_genus_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
# unique_location_family_core_microbiome_list <- location_family_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
# unique_location_order_core_microbiome_list <- location_order_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
# unique_location_class_core_microbiome_list <- location_class_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
# unique_location_phylumn_core_microbiome_list <- location_phylumn_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)

# # TODO - still to look at - Subset to those taxa that are "core" in at least 10 locations
# common_location_genus_core_microbiome_list <- location_genus_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
# common_location_family_core_microbiome_list <- location_family_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
# common_location_order_core_microbiome_list <- location_order_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
# common_location_class_core_microbiome_list <- location_class_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
# common_location_phylumn_core_microbiome_list <- location_phylumn_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)


##################
# Core by Genotype
##################

genotype_cores = lapply(ranks, function(myrank){
  prevalence = preparse[[myrank]] %>%
    group_by(Corrected_pedigree, taxon, Taxaname) %>%
    summarize(prevalence = sum(present) / n(), .groups="drop") %>%
    filter(prevalence >= prevalence_level)
  return(prevalence)
})
names(genotype_cores) = ranks

# Plot and save
genoplots = lapply(ranks, function(myrank){
  mycore = genotype_cores[[myrank]]
  outfile=paste("4_CoreMicrobiome/4a_genotype_",myrank,"_core_microbiome_heatmap.jgw.png", sep="")
  plot_heatmap(mycore, myrank, xval="Corrected_pedigree", outfile=outfile)
})

