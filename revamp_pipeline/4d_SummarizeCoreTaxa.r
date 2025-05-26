#! /usr/bin/Rscript

# Summarize the core microbiome across all samples into a table

library(dplyr)
library(tidyverse)
library(phyloseq)


# Load data 
prevalence = read.csv("4_CoreMicrobiome/4c_prevalence_summary.overall.csv")
asvs = readRDS("3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds") # To be consistent with GxE analysis

# Filter for core
core = filter(prevalence, is_core)

# Make short label
core = core %>%
  mutate(label = sub(taxon_name, pattern=".+__", repl="")) %>%
  mutate(label = paste(label, " (", prevalence*100, "%)", sep="")) %>%
  mutate(taxon_level = factor(taxon_level, levels=unique(taxon_level))) %>% # For ordering
  arrange(taxon_level, taxon_name)

# Write out (still need to do some manual formatting, unfortuntely)
write.csv(core, file="4_CoreMicrobiome/4d_short_summary.csv", row.names=FALSE)


# Determine what percentage these occupy
counts = otu_table(asvs)
totals = rowSums(counts) # Total for each ASV
taxonomy = tax_table(asvs)

levels = unique(core$taxon_level)
readcounts = lapply(levels, function(mylevel){
  # Get key for ASV table
  key = as.character(mylevel)
  taxa = taxonomy[,key] %>% sub(pattern=".+__", repl="")
  # Get names to search for
  targets = core[core$taxon_level==mylevel,"label"] 
  targets = sub(targets, pattern=" \\([0-9\\.]+%\\)", repl="")
  # Find cores
  is_core = taxa %in% targets
  core_count = sum(totals[is_core])
  noncore_count = sum(totals[!is_core])
  total = core_count + noncore_count
  # Return result
  data.frame(level=mylevel, set = c("core", "noncore"), 
             count = c(core_count, noncore_count), 
             percent = c(core_count/total, noncore_count/total))
}) %>% bind_rows() %>%
  mutate(level = factor(level, levels=levels(core$taxon_level)),
         set=factor(set, levels=c("noncore", "core"))) 


# Basic plot
ggplot(readcounts) +
  aes(x=level, y=percent, fill=set) +
  geom_col() +
  labs(x="Taxonomic level", y="Percent of Reads")



# Readcounts again, but by sample so can do a violin plot

levels = unique(core$taxon_level)
readcounts.bysample = lapply(levels, function(mylevel){
  # Get key for ASV table
  key = as.character(mylevel)
  taxa = taxonomy[,key] %>% sub(pattern=".+__", repl="")
  # Get names to search for
  targets = core[core$taxon_level==mylevel,"label"] 
  targets = sub(targets, pattern=" \\([0-9\\.]+%\\)", repl="")
  # Find cores
  is_core = taxa %in% targets
  core_count = colSums(counts[is_core])
  noncore_count = colSums(counts[!is_core])
  total = core_count + noncore_count
  # Return result
  coreset = data.frame(set="core", sample = names(core_count), count=core_count,
                       percent = core_count/total)
  noncoreset = data.frame(set="noncore", sample = names(noncore_count), count=noncore_count,
                       percent = noncore_count/total)
  # Return
  bind_rows(coreset, noncoreset) %>%
    mutate(level=mylevel) %>%
    relocate(level)

}) %>% bind_rows() %>%
  mutate(level = factor(level, levels=levels(core$taxon_level)),
         set=factor(set, levels=c("noncore", "core"))) 

# Summarize
bysample_summary = core.bysample %>%
  group_by(level) %>%
  summarize(mean_percent = mean(percent), median_percent = median(percent)) %>%
  mutate(label = paste(round(mean_percent * 100, 1), "%", sep=""))


# Plot distribution by samples
core.bysample = filter(readcounts.bysample, set=="core")
coreplot = ggplot(core.bysample) +
  aes(x=level, y=percent) +
  geom_violin(fill="lightblue") +
  labs(x="Taxonomic level", y="Percent of Reads") +
  #stat_summary(fun="mean", geom="crossbar", width=0.5) +  # Mean bar
  geom_label(data=bysample_summary, mapping=aes(x=level, y=mean_percent, label=label),
             color="black", size=3)

ggsave(coreplot, file="4_CoreMicrobiome/4d_core_percent_by_sample.png", width=4, height=3)



############
# Lineplot of decay across percentage of samples
############

# Calculate what % of taxa at each level are retained
prev.split = split(prevalence, prevalence$taxon_level)
stepdown = lapply(names(prev.split), function(mylevel){
  myprev = prev.split[[mylevel]]
  steps=(0:100) / 100
  percent = sapply(steps, function(s){
    sum(myprev$prevalence >= s) / nrow(myprev)
  })
  data.frame(level = mylevel, steps, percent)
}) %>%
  bind_rows() %>%
  mutate(level = factor(level, levels = unique(prevalence$taxon_level)))

# Plot
ggplot(stepdown) +
  aes(x=steps, y=percent, color=level) +
  geom_line()
