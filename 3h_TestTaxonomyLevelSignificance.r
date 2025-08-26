#! /usr/bin/Rscript

# Test statistical significance of taxa levels

library(dplyr)
library(gvlma)

# Load data
taxonomy = read.csv("3_GxE/3f_taxon_GXE.csv")
levels = c("Phylum","Class",  "Order",  "Family","Genus", "Species")


# Loop over the different terms
results=list()
for(myterm in unique(taxonomy$term)){

  # Fit linear model to check differences between adjacent levels
  for(i in 1:(length(levels)-1)){
    target_levels = levels[c(i, i+1)]
    mydata = taxonomy %>%
      filter(term==myterm & level %in% target_levels)
    mymodel = wilcox.test(herit ~ level, data=mydata)
    
    # Save results
    group=paste(myterm, i)
    results[[group]] = data.frame(term=myterm, level1=target_levels[1], level2=target_levels[2],
                                  test="wilcoxon", p=mymodel$p.value)
  }
}
results = bind_rows(results)


# Write out results
write.csv(results, file="3_GxE/3h_check_significance_taxa_levels.csv", row.names=FALSE)
