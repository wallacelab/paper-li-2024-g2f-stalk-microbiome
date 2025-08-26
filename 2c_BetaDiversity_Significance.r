# Load necessary libraries for beta diversity analysis
library(phyloseq)
library(tidyverse)
library(vegan)

# Checking significance of beta diversity distances via PERMANOVA

# Global variables
factors_to_test=c("location", "Corrected_pedigree")
min_group_count=10  # Minimum times the above factors have to show up to be tested for
permutations=9999
n_cores=6 # For parallel processing
set.seed(1) # Random seed


#####
# Load data
#####

# Metadata
asvs = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")
metadata = sample_data(asvs) %>% data.frame()

# Distance matrices
file_pattern=".*2b_beta_diversity\\.(.+)\\.distance.csv" # For finding & capturing distance matrix files
distance_files=list.files(path="2_Diversity/", pattern=file_pattern, full.names=TRUE)
metrics = sub(distance_files, pattern=file_pattern, repl="\\1")
distances = lapply(distance_files, function(infile){
  read.csv(infile, row.names=1, check.names=FALSE) %>% as.dist()
})
names(distances)=metrics


####
# Calculate PERMANOVA significance for each metric
####

# Helper function to calculate permanova p-value given a distance matrix and metadata
calc_permanova = function(mydist, mymeta, factor="location", min_count=10, permutations=1000, n_cores=1){

  # Convert dist to matrix (if not already)
  distmat = as.matrix(mydist)
  
  # Make sure everything in right order
  if (!identical(rownames(distmat), colnames(distmat))){
    stop("Distance matrix rows and columns do not match")
  }
  mycovar = mymeta[rownames(distmat),factor] # Get metadata correct
  
  # Remove groups with too few samples
  groupcount = table(mycovar)
  too_few = names(groupcount)[groupcount < min_count]
  samples_to_remove = rownames(mymeta)[mymeta[,factor] %in% too_few]
  distmat = distmat[!rownames(distmat) %in% samples_to_remove,
                    !colnames(distmat) %in% samples_to_remove]
  
  # Double-check order and metadata
  if (!identical(rownames(distmat), colnames(distmat))){
    stop("Distance matrix rows and columns do not match")
  }
  mycovar = mymeta[rownames(distmat),factor] # Get metadata correct
  
  # Run permanova
  distmat = as.dist(distmat) # Back to distance matrix
  permanova = adonis2(distmat ~ mycovar, permutations=permutations, parallel=n_cores)
  
  # return p-value
  p=permanova["mycovar", "Pr(>F)"]
  R2=permanova["mycovar", "R2"]
  return(data.frame(p, R2, n_groups=length(unique(mycovar))))
}

# Loop over metrics and factors to test
results=list()
for(mymetric in metrics){
  mydist=distances[[mymetric]]
  for(myfactor in factors_to_test){
    perm_results = calc_permanova(mydist, metadata, factor=myfactor,
                          min_count=min_group_count, permutations=permutations, n_cores=n_cores)
    label=paste(mymetric, myfactor)
    results[[label]] = data.frame(metric=mymetric, factor=myfactor, perm_results)
  }
}
results = bind_rows(results)

# Write out results
write.csv(results, file="2_Diversity/2c_beta_diversity.permanova_results.csv")

