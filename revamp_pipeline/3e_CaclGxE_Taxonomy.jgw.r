#! /usr/bin/Rscript

# Plot the GxE breakdown of taxonomic levels
# TODO - Only 4 species make the cut. Remove from plot?

library(qiime2R)
library(phyloseq)
library(tidyverse)

# Global variables
taxa_levels=c("Phylum", "Class","Order","Family","Genus","Species")
min_prevalence=0.2  # Minimum prevalence to be included in calculations

# Load data
asvs = readRDS("3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds")
metadata = asvs %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")

###################
# Collapse & filter
###################

# Collapse taxa at various levels
collapsed = lapply(taxa_levels, function(mylevel){
  tax_glom(asvs, taxrank=mylevel)
})
names(collapsed) = taxa_levels

# Filter for just taxa with required prevalence
core_taxa = lapply(collapsed, function(mydata){
  # Calculate prevalence of each taxon
  mytable = otu_table(mydata)
  is_present = mytable > 0
  prevalence = rowSums(is_present) / ncol(mytable)
  
  # Filter and return
  tokeep = prevalence > min_prevalence
  filtered = prune_taxa(tokeep, mydata)
  return(filtered)
})


##########
# Quick check/output of prevalence
##########

prevalence = lapply(names(core_taxa), function(level){
  mydata = otu_table(core_taxa[[level]])
  presence = mydata>0
  prevalence = rowSums(presence) / ncol(mydata)
  taxstrings = tax_table(core_taxa[[level]]) %>%
    apply(MARGIN=1, FUN=paste, collapse=";")
  output = data.frame(level=level, taxon=taxstrings, prevalence=prevalence, mydata) %>%
    arrange(taxon)
}) %>% bind_rows()
  
write.csv(prevalence, file="3_GxE/3e_taxonomy_prevalence.jgw.csv", row.names=FALSE)


############
# ANOVA
############

# Helper function to calculate ANOVA results for taxa
anova_calculation_taxa <- function(input_data, input_metadata){
  
  # Reformat core taxa table
  #core_taxa_table <- as.data.frame(t(input_data))
  core_taxa_table = otu_table(input_data) %>%
    t() %>%
    as.data.frame()
  pseudo_core_taxa_table <- core_taxa_table + 1  # Add pseudocount of 1
  pseudo_core_taxa_table$SampleID <- rownames(pseudo_core_taxa_table)
  
  # Join with metadata
  input_table_core_taxa_meatadata_results <- dplyr::inner_join(input_metadata,pseudo_core_taxa_table,by="SampleID")

  # Run linear regression 
  taxa = taxa_names(input_data)
  models = lapply(taxa, function(taxon){
    #myformula = paste(taxon, "~ location + Corrected_pedigree + location:Corrected_pedigree")
    lm(input_table_core_taxa_meatadata_results[,taxon] ~ location + Corrected_pedigree + location:Corrected_pedigree,
       data = input_table_core_taxa_meatadata_results)
  })
  
  # Run ANOVAs
  anovas <- lapply(models, anova)
  names(anovas) = taxa
  return(anovas)
}

anovas = lapply(core_taxa, anova_calculation_taxa, input_metadata=metadata)


#######################
# Heritability
#######################

# Helper function to calculate heritabilities
calc_heritability <- function(anova_results){
  
  taxa = names(anova_results)
  herits = lapply(taxa, function(taxon){
    herit = anova_results[[taxon]] %>% 
      data.frame() %>%
      rename("SS" = "Sum.Sq", "pval" = "Pr..F.") %>%  # Rename
      dplyr::select(SS, pval) %>%
      rownames_to_column("term") %>%
      mutate(herit = SS / sum(SS)) %>% # Calculate heritability
      filter(term != "Residuals") %>% # Remove residuals
      mutate(term = sapply(term, switch, "location"="Environment",
                           "Corrected_pedigree" = "Maize Genotype",
                           "location:Corrected_pedigree" = "GXE"),
             taxon=taxon) %>%
      relocate(taxon)
  }) %>% bind_rows() 
  return(herits)
  
}

# Calculate heritability
herits = lapply(anovas, calc_heritability)

# Add taxon levels
for(level in names(herits)){
  # Column with level we're at
  herits[[level]] = herits[[level]] %>%
    mutate(level=level) %>%
    relocate(level)
  
  # Actual taxonomy string
  mytaxa = tax_table(core_taxa[[level]])
  taxstrings = apply(mytaxa, MARGIN=1, FUN=paste, collapse=";")
  taxkey = data.frame(taxon=names(taxstrings), taxonomy=taxstrings)
  herits[[level]] = left_join(herits[[level]], taxkey, by="taxon")
}

# Combine into a single data frame
herits = bind_rows(herits)

############
# Plot results
############

# Calculate fractions that are significant
significant = herits %>%
  group_by(level, term) %>%
  summarize(n=n(), fraction = sum(pval < 0.01)/n(), .groups="drop") %>%
  mutate(percent = round(fraction*100, digits=1) %>% paste("%", sep=""))
  #mutate(percent = round(fraction*100, digits=1) %>% paste("% sig.", sep="")


# Plot
herits$level = factor(herits$level, levels=taxa_levels)
core_taxa_plot <- ggplot(herits) +
  aes(x=level, y=herit, fill=term) +
  geom_violin(trim = TRUE) +
  ylab("variance explained") + 
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7))  + 
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey28"), 
    panel.background = element_blank(),
    axis.title.y = element_text(size=30, face="bold"),
    plot.title = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=15,angle=90, face="bold"),
    legend.position="none")  + 
  scale_fill_discrete(name = "",) +
  facet_grid(~term) +
  theme(strip.text.x = element_text(size = 30,face="bold"),
    strip.text.y = element_text(size = 30,face="bold" )) + 
  stat_summary(fun = "mean", geom = "crossbar",  width = 0.5, colour = "black") +
  geom_text(data=significant, mapping=aes(label=percent), y=max(herits$herit)* 1.02) +
  geom_text(data=significant, mapping=aes(label=n), y=max(herits$herit) * 1.05) +
  ylim(NA, max(herits$herit) * 1.1) # Adjust vertical limit

# Save
ggsave(core_taxa_plot, file="3_GxE/3e_taxon_GXE.jgw.png", height=8, width=12)
write.csv(herits, file="3_GxE/3e_taxon_GXE.jgw.csv", row.names=FALSE)



