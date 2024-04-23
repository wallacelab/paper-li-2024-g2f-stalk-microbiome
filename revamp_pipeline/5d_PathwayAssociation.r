#! /usr/bin/Rscript

# Run associations between imputed biochemical pathways

library(dplyr)
library(Maaslin2) # NOTE: Must be Github version, repo version threw an error making heatmaps (see https://forum.biobakery.org/t/xtfrm-error-with-maaslin2-default-example-in-r/5216/8)
library(tibble)
library(readr)
library(tidyverse)
library(phyloseq)

# TODO: In preproces_pathway_description(), it says 80% prevalence filter but check is 70. 
#       Double-check math and convert to a parameter
# TODO: Environmental variables (inc. 2 to focus on) are hardcoded. Able to do that more programmatically?

# File locations
pathway_file = "0_data_files/path_abun_unstrat_descrip.tsv"
soil_file = "0_data_files/g2f_2019_soil_data.csv"
weather_file = "0_data_files//G2F_2019_weather_average_and_sum.tsv"

# Load metadata
metadata = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.phyloseq.rds") %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")

# Environmental factors to focus on
environmental_factors <- c("X1.1.Soil.pH", "WDRF.Buffer.pH", "X1.1.S.Salts.mmho.cm", "Organic.Matter.LOI..", "Nitrate.N.ppm.N", "lbs.N.A",                   
                           "Potassium.ppm.K", "Sulfate.S.ppm.S", "Calcium.ppm.Ca", "Magnesium.ppm.Mg", "CEC.Sum.of.Cations.me.100g", "Temperature..C.", "Relative.Humidity....", "Rainfall..mm.", "Solar.Radiation..W.m2.")
target_factors = c("Relative.Humidity....", "Potassium.ppm.K")

#################
# Functions
#################

# Function to preprocess the pathway description data
preprocess_pathway_description <- function(filepath) {
  pathway_description <- read.table(filepath, sep = "\t", header = TRUE, check.names = TRUE, row.names = 1)
  colnames(pathway_description) <- gsub("\\.", "-", colnames(pathway_description))
  pathway_description$description <- NULL  # Remove description to avoid console crash
  pathway_description <- pathway_description[rowSums(pathway_description == 0) <= 70, ]  # Keep pathways with >= 80% prevalence
  return(pathway_description)
}

# Function to merge pathway data with metadata and calculate average by location
calculate_average_pathway_by_location <- function(pathway_description, metadata) {
  input_pathway_table <- t(pathway_description) %>% as.data.frame()
  input_pathway_table$SampleID <- rownames(input_pathway_table)
  
  #metadata$SampleID <- gsub("\\.", "-", rownames(metadata)) # Not needed
  merged_data <- inner_join(input_pathway_table, metadata, by = "SampleID")
  
  #average_pathway <- aggregate(. ~ location, data = merged_data, FUN = mean)
  
  #average the pathway based on location 
  average_pathway_location <- aggregate(x= merged_data[,1:280],     
                                        
                                        # Specify group indicator
                                        by = list(merged_data$location),      
                                        
                                        # Specify function (i.e. mean)
                                        FUN = mean)
  
  colnames(average_pathway_location)[1] <- "location"
  return(average_pathway_location)
}

###########################################################################

# Function to read and merge environmental data with pathway data
merge_environmental_data <- function(pathway_data, soil_file, weather_file) {
  # Load soil and weather data
  soil_data <- read.csv(soil_file, sep = ",", check.names = TRUE)
  weather_data <- read.csv(weather_file, sep = "\t", check.names = TRUE)
  
  # Merge soil data with alpha diversity data already in pathway data
  pathway_soil_merged <- left_join(pathway_data, soil_data, by = c("location" = "Location"))
  
  # Merge the above with weather data
  pathway_soil_weather_merged <- left_join(pathway_soil_merged, weather_data, by = "location")
  
  # Clean up unnecessary columns, if required
  pathway_soil_weather_merged$Texture <- NULL
  pathway_soil_weather_merged$Texture.No <- NULL
  
  return(pathway_soil_weather_merged)
}

#####################################################################################


# Function to perform Maaslin2 analysis for each environmental factor
perform_maaslin2_analysis <- function(pathway_data, metadata, environmental_factors) {
  significant_results <- data.frame()
  
  for (factor in environmental_factors) {
    save_path <- paste0("5_Associations/5d_pathways/lm_ec_pathway_Maaslin2_no_transform_", factor, "_test")
    fit_data <- Maaslin2(input_data = pathway_data, input_metadata = metadata, output = save_path,
                         transform = "None", fixed_effects = factor, min_prevalence = 0.9, max_significance = 0.1,
                         standardize = TRUE)
    
    result_path <- file.path(save_path, "significant_results.tsv")
    if (file.exists(result_path)){
      temp_results <- read.table(result_path, sep = "\t", header = TRUE) %>% 
        mutate(metadata = as.character(metadata), feature = as.character(feature), value = as.character(value))
      significant_results <- bind_rows(significant_results, temp_results)
    }
  }
  
  significant_results$feature <- gsub("\\.", "-", significant_results$feature)
  return(significant_results)
}

##############################################################################

# Function to annotate significant results with pathway descriptions
annotate_significant_results <- function(significant_results, annotation_path) {
  pathway_annotation <- read.table(annotation_path, sep = "\t", header = TRUE, check.names = TRUE) %>%
    select(pathway, description)
  
  annotated_results <- left_join(significant_results, pathway_annotation, by = c("feature" = "pathway"))
  annotated_results <- annotated_results %>% 
    filter(qval <= 0.05) %>% 
    select(Factor = metadata, Pathway = description, Association = coef, Pathway_ID = feature) %>%
    mutate(Association = ifelse(coef <= 0, "negative", "positive"))
  return(annotated_results)
}



######################
# Main execution block
######################

pathway_description_1 <- preprocess_pathway_description(pathway_file)
average_pathway_location_1 <- calculate_average_pathway_by_location(pathway_description_1, metadata)
average_pathway_location_metadata_2 <- merge_environmental_data(average_pathway_location_1, soil_file, weather_file)
#colnames(average_pathway_location_metadata_2)


#select pathway section
#average_pathway_group <- average_pathway_location_metadata_2[,2:281] # FRAGILE. Replaced with next 2 lines
pathways= rownames(pathway_description_1)
average_pathway_group <- average_pathway_location_metadata_2[,pathways]

#select environmental
pathway_location_metadata <- average_pathway_location_metadata_2[,environmental_factors]

# Set location as rownames
rownames(average_pathway_group) <- average_pathway_location_metadata_2$location
rownames(pathway_location_metadata) <- average_pathway_location_metadata_2$location

#run masslin2 analysis
significant_results <- perform_maaslin2_analysis(average_pathway_group,pathway_location_metadata, environmental_factors)
significant_result_0.05 <- subset(significant_results,qval <=0.05)
significant_result_0.01  <- subset(significant_results,qval <=0.01)


# Read and preprocess the pathway description
pathway_description <- read.table(pathway_file, sep = "\t", header = TRUE, check.names = TRUE)
pathway_description <- pathway_description[, c("pathway", "description")]

# Filter the significant results to include only specific values and add pathway descriptions
significant_result_0.05_annotation <- significant_result_0.05 %>%
  filter(value %in% target_factors) %>%
  left_join(pathway_description, by = c("feature" = "pathway")) %>%
  mutate(
    association = case_when(
      coef <= 0 ~ "negative",
      coef > 0  ~ "positive"
    ),
    Factor = case_when(
      metadata == "Potassium.ppm.K"         ~ "Potassium",
      metadata == "Relative.Humidity...."   ~ "Relative Humidity"
    )
  ) %>%
  select(Factor, Pathway = description, association, Pathway_ID = feature)

# View the updated significant results with annotations
# print(significant_result_0.05_annotation)

# Optionally, save the annotated data using write.table for compatibility
write.table(significant_result_0.05_annotation, "5_Associations/5d_pathways_combined_significant_result_rerun_final.tsv", row.names = FALSE, sep = "\t", quote = FALSE)

#########################################################################################

