#! /usr/bin/Rscript

# Calculate alpha diversity on the samples
library(phyloseq)
library(tidyverse)
library(gridExtra)

# Note: this file differs by including different samples
#   For locations, it uses all samples from that location
#   For genotypes, it uses those present at least X times, not just YS filtered

min_n = 10 # Only use samples (locations/genotypes) that have at least 10 data points

#############
# Load data
#############

#ps.rarefied = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.yellow_stripe.phyloseq.rds")
ps.rarefied = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")

##############
# Helper function for filtering locations/genotypes
##############

filter_category=function(x, column, n){
  counts = table(x[,column])
  targets = names(counts)[counts >= n]
  return(subset(x, x[,column] %in% targets))
}

##############
# Alpha diversiy calculations
##############

#calculate the alpha diversity
alphaObserved = estimate_richness(ps.rarefied, measures="Observed")
alphaSimpson = estimate_richness(ps.rarefied, measures="Simpson")
alphaShannon = estimate_richness(ps.rarefied, measures="Shannon")

#merge the alpha diversity and metadata 
G2F_metadata_2019 = cbind(sample_data(ps.rarefied), alphaObserved, alphaSimpson, alphaShannon)
G2F_locs = G2F_metadata_2019 %>% filter_category(column="location", n=min_n)
G2F_genos = G2F_metadata_2019 %>% filter_category(column="Corrected_pedigree", n=min_n)

#use kruskal wallis test to test for significant difference between location and pedigree 
#  Note: unclass() is to remove the 'h1test' class so results can be put into a data frame
l1 = kruskal.test(Observed~location, data = G2F_locs) %>% unclass()
l2 = kruskal.test(Simpson~location, data = G2F_locs) %>% unclass()
l3 = kruskal.test(Shannon~location, data = G2F_locs) %>% unclass()

p1 = kruskal.test(Observed~Corrected_pedigree, data = G2F_genos) %>% unclass()
p2 = kruskal.test(Simpson~Corrected_pedigree, data = G2F_genos) %>% unclass()
p3 = kruskal.test(Shannon~Corrected_pedigree, data = G2F_genos) %>% unclass()

# Combine results into a single data frame to write out
kruskal.results = bind_rows(l1, l2, l3, p1, p2, p3) %>%
  as.data.frame() %>%
  relocate(method, data.name, statistic, parameter, p.value) %>%
  rename(dr=parameter, comparison=data.name)
write.csv(kruskal.results, file="2_Diversity/2a_alpha_diversity.kruskal_test_results.jgw.csv", row.names=FALSE)


##########################
# plot the alpha diversity 
##########################

# Subset phyloseq object to what kept
ps.locs = prune_samples(sample_data(ps.rarefied)$location %in% G2F_locs$location, ps.rarefied)
ps.genos = prune_samples(sample_data(ps.rarefied)$Corrected_pedigree %in% G2F_genos$Corrected_pedigree, ps.rarefied)


# Standard theme elements for all plots
mytheme = theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  plot.title = element_text(size=20, face="bold"),
  axis.text.y = element_text(size=20, face="bold"),
  axis.text.x = element_text(size=16, face="bold"),
  strip.text.x = element_text(size=20,face="bold"),
  panel.background = element_blank()
)

Shannon_location <- plot_richness(ps.locs, x="location", measures=c("Shannon")) +
  geom_boxplot(fill="darkred") +
  mytheme +
  labs(x="Location")

Observed_location <- plot_richness(ps.locs, x="location", measures=c("Observed")) + 
  geom_boxplot(fill="darkred") +
  mytheme +
  labs(x="Location")

Simpson_location <- plot_richness(ps.locs, x="location", measures=c("Simpson")) +
  geom_boxplot(fill="darkred") + 
  mytheme +
  labs(x="Location")

Shannon_pedigree <- plot_richness(ps.genos, x="Corrected_pedigree", measures=c("Shannon")) + 
  geom_boxplot(fill="darkgreen") + 
  mytheme +
  labs(x="Genotype")

Observed_pedigree <- plot_richness(ps.genos, x="Corrected_pedigree", measures=c("Observed")) +
  geom_boxplot(fill="darkgreen") + 
  mytheme +
  labs(x="Genotype")

Simpson_pedigree <- plot_richness(ps.genos, x="Corrected_pedigree", measures=c("Simpson")) + 
  geom_boxplot(fill="darkgreen") + 
  mytheme+
  labs(x="Genotype")

# Save plots
alphaplots = grid.arrange(grobs=list(Shannon_location, Observed_location, Simpson_location,
                                     Shannon_pedigree, Observed_pedigree, Simpson_pedigree),
                          nrow=2)
ggsave(alphaplots, file="2_Diversity/2a_alpha_diversity.plots.jgw.png", height=12, width=15)


################
# Export alpha diversity
################

alpha_output = G2F_metadata_2019 %>%
  rownames_to_column("sampleID") %>%
  select(sampleID, location, rep_number, plot_number, Corrected_pedigree, Observed, Simpson, Shannon)
write.csv(alpha_output, file="2_Diversity/2a_alpha_diversity.jgw.csv", row.names=FALSE)
