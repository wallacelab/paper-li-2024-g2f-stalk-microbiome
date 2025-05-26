#! /usr/bin/Rscript

# Compare old versus new GxE taxonomy calculations to see why some so different
#   Fixed problem 1 (3 Dec 2024): Roy's old Order and Class files had the incorrect names; corrected
#   Not a problem per se, but old "species" was mostly just genus regurgitated, with species of "__". Use new one instead
#     Pretty much all other "missing" in new analysis are same thing, so new analysis is more correct

# Heritabilities also highly correlated for shared taxa (97-99%), so yes, is consistent

library(tidyverse)
library(ggpubr)

# Load data
old.prevalence = read.csv("3_GxE/3e_taxonomy_prevalence.csv")
new.prevalence = read.csv("3_GxE/3e_taxonomy_prevalence.jgw.csv")
old.herit = read.csv("3_GxE/3e_taxon_GXE.csv")
new.herit = read.csv("3_GxE/3e_taxon_GXE.jgw.csv")

# Helper function to standardize taxonomy strings (mostly standardize to QIIME)
fixnames = function(s){
  s = gsub(s, pattern=" +", repl="") %>% # Remove spaces
    gsub(pattern=";NA", repl="") # Remove missing levels
}
new.prevalence$taxon = fixnames(new.prevalence$taxon)
new.prevalence$level = tolower(new.prevalence$level) # Lowercase
new.herit$taxon = fixnames(new.herit$taxonomy) # Overwrite 'taxon' column to make standard
new.herit$level = tolower(new.herit$level)


# Compare prevalence
old.short = select(old.prevalence, level, taxon, prevalence) %>%
  rename(old_prev = prevalence)
new.short = select(new.prevalence, level, taxon, prevalence) %>%
  rename(new_prev = prevalence)
prevalence = full_join(old.short, new.short, by=c("level", "taxon")) %>%
  arrange(level, taxon)

# Plot prevalences
corplot = ggplot(prevalence) +
  aes(x=old_prev, y=new_prev) +
  geom_smooth() + 
  geom_point() +
  stat_cor() +
  labs(title="Compare old and new prevalences",
       x="Old prevalence",
       y="New prevalence")
ggsave(corplot, file="3_GxE/3e1_compare_prevalence.png", width=5, height=5)


# Compare heritabilities
old.herit = old.herit %>%
  select(level, taxon, term, herit) %>%
  rename(old_herit=herit)
new.herit = new.herit %>%
  select(level, taxon, term, herit) %>%
  rename(new_herit=herit)
herits = full_join(old.herit, new.herit, by=c("level", "taxon", "term")) %>%
  arrange(level, taxon, term)

# Plot heritabilities
heritplot = ggplot(herits) +
  aes(x=old_herit, y=new_herit, color=term) +
  geom_smooth() + 
  geom_point() +
  stat_cor() +
  labs(title="Compare old and new prevalences",
       x="Old prevalence",
       y="New prevalence") +
  facet_wrap(~term, scales="free")
ggsave(heritplot, file="3_GxE/3e1_compare_herits.png", width=10, height=5)
