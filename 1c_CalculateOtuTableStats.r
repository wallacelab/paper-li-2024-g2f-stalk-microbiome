#! /usr/bin/Rscript

# Quick calculation of OTU table stats for the manuscript
library(phyloseq)
library(tidyverse)

#############
# Load data
#############

initial=readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.phyloseq.rds")
rarefied=readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")
yellow_stripe=readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.yellow_stripe.phyloseq.rds")

stats=function(mydata, label="UNKNOWN"){
  cat(label,"has:\n")
  cat("\t", nsamples(mydata), "samples\n")
  cat("\t", ntaxa(mydata), "taxa (OTUs)\n")
  cat("\t", sum(sample_sums(mydata)), "total read counts\n")
  cat("\t", median(sample_sums(mydata)), "median read count per sample\n")
  
  # Locations and genotypes
  genocounts=sample_data(mydata)$Corrected_pedigree %>% table()
  genocounts = genocounts[genocounts>=10]
  cat("\t", length(unique(sample_data(mydata)$Corrected_pedigree) 
                   %>% na.omit()), "genotypes\n")
  cat("\t\t", length(genocounts), "of which have at least 10 instances\n")
  cat("\t", length(unique(sample_data(mydata)$location)
                   %>% na.omit()), "locations\n")
  
  
}

stats(initial, "Initial")
stats(rarefied, "Rarefied")
stats(yellow_stripe, "Yellow Stripe")


