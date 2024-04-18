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
# TODO - Jason also disagrees with the workflow here. Seems QIIME was used to make separate
#        outputs at different prevalence instead of putting prevalence here. Redo so it
#        subsets a provided Phyloseq object instead?

# Global variables
min_prevalence = 0.6 # Fraction of samples a taxon has to be in to be "core"

# Load data
metadata = readRDS("3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds") %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")


#######################################################################

location_core_microbiome <- function(taxonomic_level,prevalence_level){
  
  #create an empty dataframe 
  core_taxa_list = data.frame(Taxaname=character(), location=character())
  
  #iterate through each location and get the given taxonomic-level taxa that have input prevalence 
  
  for(experiment_field in unique(metadata$location)){
    
    print(experiment_field)
    #create read path
    directory_name = paste("./0_data_files/core_metrics_by_location/core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-",experiment_field,"-1500",sep="")
    #print(directory_name)
    path_name = paste(directory_name,"/rarefied_table_",taxonomic_level,"_level_taxonomy_",prevalence_level,"_filtered.qza",sep="")
    #print(path_name)
    #read into qiime 2 qza file and take the data 
    core_taxa <- read_qza(path_name)
    core_taxa <- as.data.frame(core_taxa$data)
    
    # Only add if actually have core taxa
    if(nrow(core_taxa)>0){
      core_taxa$Taxaname <- rownames(core_taxa)
      core_taxa$location <- experiment_field
      #record down the name of core taxa and bind wih
      location_core_taxa <- core_taxa[,c("Taxaname","location")]
      core_taxa_list <- dplyr::bind_rows(core_taxa_list,location_core_taxa)
    }
  }
  
  length(unique(core_taxa_list$Taxaname))
  length(unique(core_taxa_list$location))
  
  return(core_taxa_list)
  
}

core_microbiome_list_to_matrix <- function(core_taxa_list){
  
  #build a dataframe based on core taxa list of all location 
  core_taxa_matrix <- as.data.frame(matrix(ncol = length(unique(core_taxa_list$location)), nrow = length(unique(core_taxa_list$Taxaname))))
  colnames(core_taxa_matrix) <- unique(core_taxa_list$location)
  rownames(core_taxa_matrix) <- unique(core_taxa_list$Taxaname)
  
  #make a matrix table of core taxa, if core taxa exist in a location, add 1 
  for(taxa_name in unique(core_taxa_list$Taxaname)){
    print(taxa_name)
    temp_location_taxa_combination <- subset(core_taxa_list,Taxaname == taxa_name)
    for(location in temp_location_taxa_combination$location){
      print(location)
      core_taxa_matrix[taxa_name,location] <- 1
    }
  }
  
  # if core taxa do not present in a location, change na to 0 
  core_taxa_matrix[is.na(core_taxa_matrix)] <- 0
  
  return(core_taxa_matrix)
  
}
#core microbiome 

##############################################################################################

location_species_core_microbiome_list <- location_core_microbiome(7, min_prevalence)
location_genus_core_microbiome_list <- location_core_microbiome(6, min_prevalence)
location_family_core_microbiome_list <- location_core_microbiome(5, min_prevalence)
location_order_core_microbiome_list <- location_core_microbiome(4, min_prevalence)
location_class_core_microbiome_list <- location_core_microbiome(3, min_prevalence)
location_phylumn_core_microbiome_list <- location_core_microbiome(2, min_prevalence)

# Subset to just those taxa that are "core" in a single location
unique_location_genus_core_microbiome_list <- location_genus_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
unique_location_family_core_microbiome_list <- location_family_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
unique_location_order_core_microbiome_list <- location_order_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
unique_location_class_core_microbiome_list <- location_class_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
unique_location_phylumn_core_microbiome_list <- location_phylumn_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)

# Subset to those taxa that are "core" in at least 10 locations
common_location_genus_core_microbiome_list <- location_genus_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
common_location_family_core_microbiome_list <- location_family_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
common_location_order_core_microbiome_list <- location_order_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
common_location_class_core_microbiome_list <- location_class_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
common_location_phylumn_core_microbiome_list <- location_phylumn_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)


# TODO - Start here. Can these be changed to simple pivot_wider() commands? Never used later anyway
location_genus_core_microbiome_matrix <- core_microbiome_list_to_matrix(location_genus_core_microbiome_list)
#location_genus_core_microbiome_matrix

location_family_core_microbiome_matrix <- core_microbiome_list_to_matrix(location_family_core_microbiome_list)
#location_family_core_microbiome_matrix


location_family_core_microbiome_list[c('Domain','Phylumn','Class','Order','Family')]  <- str_split_fixed(location_family_core_microbiome_list$Taxaname, ';',5)


location_family_core_microbiome_heatmap <-ggplot(location_family_core_microbiome_list, aes(x=location,y=Family)) + 
  geom_tile() +  theme(
    axis.text.y = element_text(size=15, face="bold"),
    axis.text.x = element_text(size=15,angle=90, face="bold"),
    axis.title.x = element_text(size=15,face="bold"),
    axis.title.y = element_text(size=15, face="bold")
    
  )

location_family_core_microbiome_heatmap



ggsave("4_CoreMicrobiome/4a_location_family_core_microbiome_heatmap.png",
       plot = location_family_core_microbiome_heatmap, height=8, width=12, device="png")



##########################################################################################



##############################################################################################
# Corrected_pedigree_phylumn_60_core_taxa_list <- as.data.frame(matrix(ncol = 2, nrow = 0))
# #add colname 
# colnames(Corrected_pedigree_phylumn_60_core_taxa_list) <- c("Taxaname","Corrected_pedigree")
# #change column type so that can be used later 
# Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname <- as.character(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname)
# Corrected_pedigree_phylumn_60_core_taxa_list$Corrected_pedigree <- as.character(Corrected_pedigree_phylumn_60_core_taxa_list$Corrected_pedigree)

Corrected_pedigree_phylumn_60_core_taxa_list = data.frame(Taxaname = character(),Corrected_pedigree = character())

#iterate through each Corrected_pedigree and get family level taxa that have 60% prevalence 

for(pedigree in unique(metadata$Corrected_pedigree)){
  
  #print(pedigree)
  pedigree = gsub("\\/","X",pedigree)
  print(pedigree)
  #create read path
  directory_name = paste("0_data_files/core_metrics_by_pedigree/core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-",pedigree,"-1500",sep="")
  #print(directory_name)
  path_name = paste(directory_name,"/rarefied_table_5_level_taxonomy_0.6_filtered.qza",sep="")
  #print(path_name)
  #read into qiime 2 qza file and take the data 
  phylumn_60_filtered_core_taxa <- read_qza(path_name)
  phylumn_60_filtered_core_taxa <- as.data.frame(phylumn_60_filtered_core_taxa$data)
  phylumn_60_filtered_core_taxa$Taxaname <- rownames(phylumn_60_filtered_core_taxa)
  phylumn_60_filtered_core_taxa$Corrected_pedigree <- pedigree
  #record down the name of core taxa and bind wih
  Corrected_pedigree_core_taxa <- phylumn_60_filtered_core_taxa[,c("Taxaname","Corrected_pedigree")]
  Corrected_pedigree_phylumn_60_core_taxa_list <- dplyr::bind_rows(Corrected_pedigree_phylumn_60_core_taxa_list,Corrected_pedigree_core_taxa)
}

# length(unique(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname))
# length(unique(Corrected_pedigree_phylumn_60_core_taxa_list$Corrected_pedigree))

#build a dataframe based on core taxa list of all Corrected_pedigree 
Corrected_pedigree_phylumn_60_core_taxa_matrix <- as.data.frame(matrix(ncol = length(unique(Corrected_pedigree_phylumn_60_core_taxa_list$Corrected_pedigree)), nrow = length(unique(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname))))
colnames(Corrected_pedigree_phylumn_60_core_taxa_matrix) <- unique(Corrected_pedigree_phylumn_60_core_taxa_list$Corrected_pedigree)
rownames(Corrected_pedigree_phylumn_60_core_taxa_matrix) <- unique(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname)

#Corrected_pedigree_phylumn_60_core_taxa_matrix["d__Bacteria;p__Proteobacteria","SCH1"]

#make a matrix table of core taxa, if core taxa exist in a Corrected_pedigree, add 1 
for(taxa_name in unique(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname)){
  #print(taxa_name)
  temp_Corrected_pedigree_taxa_combination <- subset(Corrected_pedigree_phylumn_60_core_taxa_list,Taxaname == taxa_name)
  for(Corrected_pedigree in temp_Corrected_pedigree_taxa_combination$Corrected_pedigree){
    print(Corrected_pedigree)
    Corrected_pedigree_phylumn_60_core_taxa_matrix[taxa_name,Corrected_pedigree] <- 1
  }
}

# if core taxa do not present in a Corrected_pedigree, change na to 0 
Corrected_pedigree_phylumn_60_core_taxa_matrix[is.na(Corrected_pedigree_phylumn_60_core_taxa_matrix)] <- 0
Corrected_pedigree_phylumn_60_core_taxa_matrix

Corrected_pedigree_phylumn_60_core_taxa_list

Corrected_pedigree_phylumn_60_core_taxa_list[c('Domain','Phylumn','Class','Order','Family')]  <- str_split_fixed(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname, ';',5)


Corrected_pedigree_family_core_microbiome_heatmap <-ggplot(Corrected_pedigree_phylumn_60_core_taxa_list, aes(x=Corrected_pedigree,y=Family)) + 
  geom_tile() +  theme(
    axis.text.y = element_text(size=15, face="bold"),
    axis.text.x = element_text(size=15,angle=90, face="bold"),
    axis.title.x = element_text(size=15,face="bold"),
    axis.title.y = element_text(size=15, face="bold")
    
  )

Corrected_pedigree_family_core_microbiome_heatmap

ggsave("4_CoreMicrobiome/4a_genotype_family_core_microbiome_heatmap.png",
       plot = location_family_core_microbiome_heatmap, height=8, width=12, device="png")


