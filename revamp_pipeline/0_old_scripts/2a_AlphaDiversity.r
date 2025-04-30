#! /usr/bin/Rscript

# Calculate alpha diversity on the samples
library(phyloseq)
library(tidyverse)
library(gridExtra)

# TODO - Jason thinks the location comparisons should use the full dataset, not the YS_filtered
# TODO - Jason thinks the genotype comparisons should use ones present at least X times, not just YS filtered

#############
# Load data
#############

ps.rarefied_ys_filtered = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.yellow_stripe.phyloseq.rds")

##############
# Alpha diversiy calculations
##############

#calculate the alpha diversity
alphaObserved = estimate_richness(ps.rarefied_ys_filtered, measures="Observed")
alphaSimpson = estimate_richness(ps.rarefied_ys_filtered, measures="Simpson")
alphaShannon = estimate_richness(ps.rarefied_ys_filtered, measures="Shannon")

#merge the alpha diversity and metadata 
G2F_metadata_2019 = cbind(sample_data(ps.rarefied_ys_filtered), alphaObserved, alphaSimpson, alphaShannon)

#use kruskal wallis test to test for significant difference between location and pedigree 
#  Note: unclass() is to remove the 'h1test' class so results can be put into a data frame
l1 = kruskal.test(Observed~location, data = G2F_metadata_2019) %>% unclass()
l2 = kruskal.test(Simpson~location, data = G2F_metadata_2019) %>% unclass()
l3 = kruskal.test(Shannon~location, data = G2F_metadata_2019) %>% unclass()

p1 = kruskal.test(Observed~Corrected_pedigree, data = G2F_metadata_2019) %>% unclass()
p2 = kruskal.test(Simpson~Corrected_pedigree, data = G2F_metadata_2019) %>% unclass()
p3 = kruskal.test(Shannon~Corrected_pedigree, data = G2F_metadata_2019) %>% unclass()

# Combine results into a single data frame to write out
kruskal.results = bind_rows(l1, l2, l3, p1, p2, p3) %>%
  as.data.frame() %>%
  relocate(method, data.name, statistic, parameter, p.value) %>%
  rename(dr=parameter, comparison=data.name)
write.csv(kruskal.results, file="2_Diversity/2a_alpha_diversity.kruskal_test_results.csv", row.names=FALSE)

# # Check if more complex formulas show a pattern instead; have to use ANOVA
# # NOTE: Observed is more normal if log-transformed, Shannon as is, Simpson very non-normal both ways
# both1 = lm(Observed~location * Corrected_pedigree, data = G2F_metadata_2019) %>% anova()
# both2 = lm(Simpson~location * Corrected_pedigree, data = G2F_metadata_2019) %>% anova()
# both3 = lm(Shannon~location * Corrected_pedigree, data = G2F_metadata_2019) %>% anova()
# 
# # No interactions because lots of missing; lots of gvlma violations
# m1 = lm(Observed~location + Corrected_pedigree, data = G2F_metadata_2019) %>% gvlma()
# m2 = lm(Simpson~location + Corrected_pedigree, data = G2F_metadata_2019) %>% gvlma()
# m3 = lm(Shannon~location + Corrected_pedigree, data = G2F_metadata_2019) %>% gvlma()


##########################
# plot the alpha diversity 
##########################

# Standard theme elements for all plots
mytheme = theme(
  legend.text = element_text(color = "black", size = 20),
  legend.title = element_text(color = "black", size = 20),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  plot.title = element_text(size=20, face="bold"),
  axis.text.y = element_text(size=20, face="bold"),
  axis.text.x = element_text(size=20, face="bold"),
  strip.text.x = element_text(size=20,face="bold"),
  panel.background = element_blank()
)

Shannon_location <- plot_richness(ps.rarefied_ys_filtered, x="location", measures=c("Shannon")) +
  geom_boxplot(fill="darkred") +
  mytheme

Observed_location <- plot_richness(ps.rarefied_ys_filtered, x="location", measures=c("Observed")) + 
  geom_boxplot(fill="darkred") +
  mytheme

Simpson_location <- plot_richness(ps.rarefied_ys_filtered, x="location", measures=c("Simpson")) +
  geom_boxplot(fill="darkred") + 
  mytheme

Shannon_pedigree <- plot_richness(ps.rarefied_ys_filtered, x="Corrected_pedigree", measures=c("Shannon")) + 
  geom_boxplot(fill="darkgreen") + 
  mytheme

Observed_pedigree <- plot_richness(ps.rarefied_ys_filtered, x="Corrected_pedigree", measures=c("Observed")) +
  geom_boxplot(fill="darkgreen") + 
  mytheme

Simpson_pedigree <- plot_richness(ps.rarefied_ys_filtered, x="Corrected_pedigree", measures=c("Simpson")) + 
  geom_boxplot(fill="darkgreen") + 
  mytheme

# Save plots
alphaplots = grid.arrange(grobs=list(Shannon_location, Observed_location, Simpson_location,
                                     Shannon_pedigree, Observed_pedigree, Simpson_pedigree),
                          nrow=2)
ggsave(alphaplots, file="2_Diversity/2a_alpha_diversity.plots.png", height=12, width=15)


################
# Export alpha diversity
################

alpha_output = G2F_metadata_2019 %>%
  rownames_to_column("sampleID") %>%
  select(sampleID, location, rep_number, plot_number, Corrected_pedigree, Observed, Simpson, Shannon)
write.csv(alpha_output, file="2_Diversity/2a_alpha_diversity.csv", row.names=FALSE)


# # Check alpha after correcting for location
# ggplot(G2F_metadata_2019) + aes(x=Corrected_pedigree, y=Observed) +
#   aes(color=Corrected_pedigree) +
#   geom_point() +
#   facet_wrap(~location)
