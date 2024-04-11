library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)



#######################################################################

location_core_microbiome <- function(taxonomic_level,prevalence_level){
  
  #create an empty dataframe 
  core_taxa_list <- as.data.frame(matrix(ncol = 2, nrow = 0))
  #add colname 
  colnames(core_taxa_list) <- c("Taxaname","location")
  #change column type so that can be used later 
  core_taxa_list$Taxaname <- as.character(core_taxa_list$Taxaname)
  core_taxa_list$location <- as.character(core_taxa_list$location)
  
  #iterate through each location and get the given taxonomic level taxa that have input prevalence 
  
  for(experiment_field in unique(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$location)){
    
    print(experiment_field)
    #create read path
    directory_name = paste("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-",experiment_field,"-1500",sep="")
    print(directory_name)
    path_name = paste(directory_name,"/rarefied_table_",taxonomic_level,"_level_taxonomy_",prevalence_level,"_filtered.qza",sep="")
    print(path_name)
    #read into qiime 2 qza file and take the data 
    core_taxa <- read_qza(path_name)
    core_taxa <- as.data.frame(core_taxa$data)
    core_taxa$Taxaname <- rownames(core_taxa)
    core_taxa$location <- experiment_field
    #record down the name of core taxa and bind wih
    location_core_taxa <- core_taxa[,c("Taxaname","location")]
    core_taxa_list <- dplyr::bind_rows(core_taxa_list,location_core_taxa)
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

location_species_core_microbiome_list <- location_core_microbiome(7,0.6)
location_genus_core_microbiome_list <- location_core_microbiome(6,0.6)
location_family_core_microbiome_list <- location_core_microbiome(5,0.6)
location_order_core_microbiome_list <- location_core_microbiome(4,0.6)
location_class_core_microbiome_list <- location_core_microbiome(3,0.6)
location_phylumn_core_microbiome_list <- location_core_microbiome(2,0.6)

unique_location_genus_core_microbiome_list <- location_genus_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
unique_location_family_core_microbiome_list <- location_family_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
unique_location_order_core_microbiome_list <- location_order_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
unique_location_class_core_microbiome_list <- location_class_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)
unique_location_phylumn_core_microbiome_list <- location_phylumn_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() <= 1)

common_location_genus_core_microbiome_list <- location_genus_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
common_location_family_core_microbiome_list <- location_family_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
common_location_order_core_microbiome_list <- location_order_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
common_location_class_core_microbiome_list <- location_class_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)
common_location_phylumn_core_microbiome_list <- location_phylumn_core_microbiome_list  %>% group_by(Taxaname) %>% filter(n() >= 10)

location_genus_core_microbiome_matrix <- core_microbiome_list_to_matrix(location_genus_core_microbiome_list)
location_genus_core_microbiome_matrix

location_family_core_microbiome_matrix <- core_microbiome_list_to_matrix(location_family_core_microbiome_list)
location_family_core_microbiome_matrix


location_family_core_microbiome_list[c('Domain','Phylumn','Class','Order','Family')]  <- str_split_fixed(location_family_core_microbiome_list$Taxaname, ';',5)


location_family_core_microbiome_heatmap <-ggplot(location_family_core_microbiome_list, aes(x=location,y=Family)) + 
  geom_tile() +  theme(
    axis.text.y = element_text(size=15, face="bold"),
    axis.text.x = element_text(size=15,angle=90, face="bold"),
    axis.title.x = element_text(size=15,face="bold"),
    axis.title.y = element_text(size=15, face="bold")
    
  )

location_family_core_microbiome_heatmap



ggsave("./location_family_core_microbiome_heatmap.png",plot = location_family_core_microbiome_heatmap,core_taxa_plot,height=8, width=12, device="png")





##########################################################################################



##############################################################################################
Corrected_pedigree_phylumn_60_core_taxa_list <- as.data.frame(matrix(ncol = 2, nrow = 0))
#add colname 
colnames(Corrected_pedigree_phylumn_60_core_taxa_list) <- c("Taxaname","Corrected_pedigree")
#change column type so that can be used later 
Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname <- as.character(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname)
Corrected_pedigree_phylumn_60_core_taxa_list$Corrected_pedigree <- as.character(Corrected_pedigree_phylumn_60_core_taxa_list$Corrected_pedigree)

#iterate through each Corrected_pedigree and get family level taxa that have 60% prevalence 

for(pedigree in unique(G2F_metadata_2019_duplicate_pedigree_ys_filtered_selected_location$Corrected_pedigree)){
  
  #print(pedigree)
  pedigree = gsub("\\/","X",pedigree)
  print(pedigree)
  #create read path
  directory_name = paste("./core-metrics-results-dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree-",pedigree,"-1500",sep="")
  print(directory_name)
  path_name = paste(directory_name,"/rarefied_table_5_level_taxonomy_0.6_filtered.qza",sep="")
  print(path_name)
  #read into qiime 2 qza file and take the data 
  phylumn_60_filtered_core_taxa <- read_qza(path_name)
  phylumn_60_filtered_core_taxa <- as.data.frame(phylumn_60_filtered_core_taxa$data)
  phylumn_60_filtered_core_taxa$Taxaname <- rownames(phylumn_60_filtered_core_taxa)
  phylumn_60_filtered_core_taxa$Corrected_pedigree <- pedigree
  #record down the name of core taxa and bind wih
  Corrected_pedigree_core_taxa <- phylumn_60_filtered_core_taxa[,c("Taxaname","Corrected_pedigree")]
  Corrected_pedigree_phylumn_60_core_taxa_list <- dplyr::bind_rows(Corrected_pedigree_phylumn_60_core_taxa_list,Corrected_pedigree_core_taxa)
}

length(unique(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname))
length(unique(Corrected_pedigree_phylumn_60_core_taxa_list$Corrected_pedigree))

#build a dataframe based on core taxa list of all Corrected_pedigree 
Corrected_pedigree_phylumn_60_core_taxa_matrix <- as.data.frame(matrix(ncol = length(unique(Corrected_pedigree_phylumn_60_core_taxa_list$Corrected_pedigree)), nrow = length(unique(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname))))
colnames(Corrected_pedigree_phylumn_60_core_taxa_matrix) <- unique(Corrected_pedigree_phylumn_60_core_taxa_list$Corrected_pedigree)
rownames(Corrected_pedigree_phylumn_60_core_taxa_matrix) <- unique(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname)

#Corrected_pedigree_phylumn_60_core_taxa_matrix["d__Bacteria;p__Proteobacteria","SCH1"]

#make a matrix table of core taxa, if core taxa exist in a Corrected_pedigree, add 1 
for(taxa_name in unique(Corrected_pedigree_phylumn_60_core_taxa_list$Taxaname)){
  print(taxa_name)
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

