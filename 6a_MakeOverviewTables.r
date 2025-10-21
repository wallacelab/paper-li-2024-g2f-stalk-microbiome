#! /usr/bin/Rscript

# Quick calculation of OTU table stats for the manuscript
library(phyloseq)
library(tidyverse)
library(argparse)

parser=ArgumentParser()
parser$add_argument("-o", "--outdir", default="6_PublicationFiles", help="Output directory")
args=parser$parse_args()

##########
# Load data
##########

initial=readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.phyloseq.rds")
rarefied=readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")
gxe = readRDS("3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds")
metadata = read.csv("0_data_files/g2f_2019_field_metadata.csv", check.names=FALSE)



#############
# Overall table stats
#############

stats=function(mydata, label="UNKNOWN"){

  # helper object for genotype counts above 10
  genocounts=sample_data(mydata)$Corrected_pedigree %>% table()
  genocounts = genocounts[genocounts>=10]
  
  outstats = c(
    "n_samples" = nsamples(mydata),
    "n_taxa_OTUs" = ntaxa(mydata),
    "total_read_counts" = sum(sample_sums(mydata)),
    "median_read_count_per_sample" = median(sample_sums(mydata)),
    "n_genotypes" = length(unique(sample_data(mydata)$Corrected_pedigree) %>% na.omit()),
    "n_genotypes_with_10_instances" = length(genocounts),
    "n_locations" = length(unique(sample_data(mydata)$location) %>% na.omit())
  ) %>% as.data.frame()
  
  # Make into data frame
  names(outstats)[1] = label
  return(outstats)
  
}

initial.counts = stats(initial, "Initial")
rarefied.counts = stats(rarefied, "Rarefied")
gxe.counts = stats(gxe, "GxE_Subset")

combined = cbind(initial.counts, rarefied.counts, gxe.counts)
write.csv(combined, file=paste(args$outdir, "/", "table.table_stats.supplemental.csv", sep=""), 
          row.names=TRUE)

##########
# Location stats
##########

# Select out field metadata
fields = metadata %>%
  select(Experiment_Code, City, 
         `Latitude_of_Field_Corner_#1 (lower left)`,
         `Longitude_of_Field_Corner_#1 (lower left)`) %>%
  rename(Location = Experiment_Code,
         Latitude=`Latitude_of_Field_Corner_#1 (lower left)`,
         Longitude=`Longitude_of_Field_Corner_#1 (lower left)`)

# Calculate counts per field
field_summary = function(mydata, label="UNKNOWN"){
  mylocs = sample_data(mydata)$location %>%
    table() %>%
    as.data.frame() %>%
    rename("Location" = ".")
  names(mylocs)[names(mylocs)=="Freq"] = label # Not sure how to do this with dplyr
  return(mylocs)
}
initial.fields = field_summary(initial, "Initial_Samples")
rarefied.fields = field_summary(rarefied, "Rarefied_Samples")
gxe.fields = field_summary(gxe, "GxE_Subset_Samples")

# Combine
outfields = fields %>%
  left_join(initial.fields, by="Location") %>%
  left_join(rarefied.fields, by="Location") %>%
  left_join(gxe.fields, by="Location")

# Replace & filter
targets = c("Initial_Samples", "Rarefied_Samples", "GxE_Subset_Samples")
outfields[,targets][is.na(outfields[,targets])] = 0 # Set missing to 0 
outfields = subset(outfields, rowSums(outfields[,targets]) > 0) # Remove rows with no samples in any set

# Output
write.csv(outfields, file=paste(args$outdir, "/", "table.field_stats.supplemental.csv", sep=""), 
          row.names=FALSE)
