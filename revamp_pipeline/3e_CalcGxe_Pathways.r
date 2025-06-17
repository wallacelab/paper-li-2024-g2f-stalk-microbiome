#! /usr/bin/Rscript

# Look at GxE for inferred biochemical pathways

# Load libraries
library(data.table)
library(tidyverse)
library(phyloseq)
library(MASS)

# Global variables
#pathway_file="0_data_files/path_abun_unstrat_descrip.tsv" # Original analysis; DEPRECATED
pathway_file="3_GxE/3d_picrust2_predictions/pathways_out/path_abun_unstrat.tsv.gz" # Redone with exact GxE sample dataset
max_sample_missing=70 # Pathway descriptions with 0s (=missing) in more than this many samples will be removed
n_diagnostics=12 # Number of pathways to plot diagnostic plots for transmforations

# Load metadata
metadata = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.phyloseq.rds") %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")

# Read and preprocess pathway descriptions
pathway_description <- read.table(pathway_file, sep = "\t", header = TRUE, check.names = TRUE, row.names = 1)
colnames(pathway_description) <- gsub("\\.", "-", colnames(pathway_description))
pathway_description$description <- NULL  # Remove description to prevent console crash
pathway_description <- pathway_description[rowSums(pathway_description == 0) <= max_sample_missing, ]
pathways = rownames(pathway_description)

# Prepare pathway data
input_pathway_table <- t(pathway_description) %>% as.data.frame()
input_pathway_table$SampleID <- rownames(input_pathway_table)
input_pathway_table_metadata <- inner_join(input_pathway_table, metadata, by = "SampleID")

# Linear regression on each pathway (old method)
# lm_results <- lapply(input_pathway_table_metadata[, pathways],
#                      function(x) lm(x ~ location + Corrected_pedigree + location:Corrected_pedigree, data = input_pathway_table_metadata))

# DEPRECATED BECAUSE boxcox requires all model elements still there apparently
# # Linear regression on each pathway (new method, more control)
# pathway_regression = function(mypath, mydata){
#   mydata = mydata[,c(mypath, "location", "Corrected_pedigree")] %>%
#     mutate(interaction=paste(location, Corrected_pedigree, sep="|"))
#   myformula = paste("`", mypath, "` ~ location + Corrected_pedigree + interaction", sep='')
#   mymodel = lm(myformula, data = mydata)
#   return(mymodel)
# }
# lm_results <- lapply(pathways, pathway_regression, mydata=input_pathway_table_metadata)
# names(lm_results) = pathways

# Fit linear model before and after boxcox transformation to each variable
pathways_transformed = input_pathway_table_metadata
lm_results=list()
lm_boxcox=list()
for(mypath in pathways){
  # Initial regression
  mydata = input_pathway_table_metadata[,c(mypath, "location", "Corrected_pedigree")] %>%
    mutate(interaction=paste(location, Corrected_pedigree, sep="|"))
  mydata[,mypath] = mydata[,mypath] - min(mydata[,mypath]) + 1 # Make sure is positive for boxcox
  myformula = paste("`", mypath, "` ~ location + Corrected_pedigree + interaction", sep='')
  lm_results[[mypath]] = lm(myformula, data = mydata)
  
  # Determine best boxcox transformation
  bc = boxcox(lm_results[[mypath]], lambda=seq(-2, 2, 0.05), plotit = FALSE)
  best = which.max(bc$y)
  lambda = bc$x[best]
  
  # Transform data
  x = mydata[,mypath]
  pathways_transformed[,mypath] = (x ^ lambda - 1)/x
  
  # Rerun regression
  mydata[,mypath] = pathways_transformed[,mypath]
  lm_boxcox[[mypath]] = lm(myformula, data = mydata)
  #break
}

# Output diagnostic plots of before and after transofrmation
set.seed(1)
diagnostics = sample(pathways, size=n_diagnostics, replace=FALSE) %>% sort()
diag_data = lapply(diagnostics, function(mypath){
#for(mypath in diagnostics){
  raw = data.frame(pathway=mypath, set="raw", value=input_pathway_table_metadata[,mypath])
  raw_resid = data.frame(pathway=mypath, set="raw_residuals", value=residuals(lm_results[[mypath]]))
  trans = data.frame(pathway=mypath, set="trans", value=pathways_transformed[,mypath])
  trans_resid = data.frame(pathway=mypath, set="trans_residuals", value=residuals(lm_boxcox[[mypath]]))
  return(bind_rows(raw, raw_resid, trans, trans_resid))
}) %>% bind_rows()
diplot = ggplot(diag_data) +
  aes(x=value, fill=set) +
  geom_histogram(bins=40) +
  facet_wrap(pathway ~ set, scales="free", ncol=4) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        legend.position = "none")
ggsave(diplot, file="3_GxE/3e_MetaCyc_pathway_GXE.diagnostics.png",
       height = n_diagnostics*2, width = 6)


# ANOVA of results
#####anova_results <- lapply(lm_results, anova)  # Don't use raw results
anova_results <- lapply(lm_boxcox, anova)  # Use boxcox-transformed instead

# Calculate heritabilities
herits = lapply(names(anova_results), function(mypathway){
  myanova = anova_results[[mypathway]]
  herit = myanova %>% 
    data.frame() %>%
    rename("SS" = "Sum.Sq", "pval" = "Pr..F.") %>%  # Rename
    dplyr::select(SS, pval) %>%
    rownames_to_column("term") %>%
    mutate(herit = SS / sum(SS)) %>% # Calculate heritability
    filter(term != "Residuals") %>% # Remove residuals
    mutate(term = sapply(term, switch, "location"="Environment",
                         "Corrected_pedigree" = "Maize Genotype",
                         #"location:Corrected_pedigree" = "GXE"),
                         "interaction" = "GXE"),
           pathway=mypathway) %>%
    relocate(pathway)
  return(herit)
}) %>%
  bind_rows() %>%
  mutate(fdr = p.adjust(pval, method = "fdr"))

# Calculate fractions that are significant
significant = herits %>%
  group_by(term) %>%
  summarize(fraction = sum(fdr < 0.01)/n()) %>%
  mutate(percent = round(fraction*100, digits=1) %>% paste("% sig."))

# Make violin plot  
heritability_plot = ggplot(herits) +
  aes(x = term, y = herit, fill = term) +
  geom_violin(trim = TRUE) +
  theme_minimal() +
  theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12)
  ) +
  labs(x="", y="Heritability") + 
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "black") +
  geom_text(data=significant, mapping=aes(label=percent), y=max(herits$herit))

# Save data & plots
ggsave(heritability_plot, file="3_GxE/3e_MetaCyc_pathway_GXE.png", height = 10, width = 10, device = "png")
write.csv(herits, file="3_GxE/3e_MetaCyc_pathway_GXE.csv", row.names=FALSE)
