# Load necessary libraries for beta diversity analysis
library(phyloseq)
library(tidyverse)
library(rbiom)
library(qiime2R) # Provides the read_qza function
library(ggpubr)

# Beta Diversity Analysis Using Weighted and Unweighted UniFrac Distances

#####
# Load data
#####

#asvs = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.yellow_stripe.phyloseq.rds")
asvs = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")
lockey = read.csv("1_parsed_files/1b_location_key.csv") # Location key for colors, etc

# Extract components for easier manipulation
counts = otu_table(asvs)
mytree = phy_tree(asvs)
metadata = sample_data(asvs) %>% data.frame()

####
# Calculate beta diversity using rbiom (because Phyloseq's unifract calculation as a bug: https://github.com/joey711/phyloseq/issues/956 )
####

betas = list()  #List to hold results
betas[['Weighted Unifrac']] = rbiom::beta.div(counts, method="unifrac", weighted=TRUE, tree=mytree)
betas[['Unweighted Unifrac']] = rbiom::beta.div(counts, method="unifrac", weighted=FALSE, tree=mytree)
betas[['Bray-Curtis']] = rbiom::beta.div(counts, method="bray-curtis")

# Calculate PCs
pcs = lapply(betas, cmdscale, eig=TRUE)

########
# Plot
########

# Combine with metadata
plotdata = lapply(pcs, function(p){
  # Label columns and make for axes
  indices = 1:ncol(p$points)
  percents = (p$eig / sum(p$eig) * 100) %>% round(digits=1)
  colnames(p$points) = paste("PC", indices, sep="")
  p$axis_labels = paste(colnames(p$points), " (", percents[indices], "%)", sep="")
  # Combine with metadata and reformat
  p$plotdata = bind_cols(p$points, metadata) %>%
    rownames_to_column("sample") %>%
    left_join(lockey, by="location") 
    # mutate(region = case_when(
    #   location %in% c("MOH1", "IAH2", "IAH4", "MNH1", "NEH1", "NEH2") ~ "Midwest",
    #   location %in% c("WIH1", "INH1", "MIH1", "OHH1") ~ "East Mississippi River",
    #   location %in% c("GAH1", "GAH2", "SCH1", "NCH1") ~ "South",
    #   location %in% c("NYH2", "NYH3", "DEH1") ~ "Northeast"
    # ))
  return(p)
})

# Make color keys
regionkey = lockey %>%
  dplyr::select(region, region_color) %>%
  unique() 
colors.region = setNames(regionkey$region_color, regionkey$region)
colors.location = setNames(lockey$location_color, lockey$location)

# Plot, colored by region
plots.region = lapply(names(plotdata), function(metric){
  mydata=plotdata[[metric]]
  ggplot(mydata$plotdata) +
    aes(x = PC1, y = PC2, color = region) +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = mydata$axis_labels[1], y=mydata$axis_labels[2],
         title=metric, color="Region") +
    #scale_colour_manual(values = c("#FFC107", "#D81B60", "#1E88E5", "#004D40"))
    scale_colour_manual(values = colors.region)
})

regionplot = ggarrange(plotlist=plots.region, nrow=1, common.legend=TRUE, legend="bottom")
ggsave(regionplot, file="2_Diversity/2b_beta_diversity.by_region.jgw.png", width=8, height=3)


# Plot, colored by sample location
plots.location = lapply(names(plotdata), function(metric){
  mydata=plotdata[[metric]]
  nloc = length(unique(mydata$plotdata$location)) # number of locations, for coloring
  ggplot(mydata$plotdata) +
    aes(x = PC1, y = PC2, color = location) +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = mydata$axis_labels[1], y=mydata$axis_labels[2],
         title=metric, color="Location") +
    scale_color_manual(values=colors.location )
})

locationplot = ggarrange(plotlist=plots.location, nrow=1, common.legend=TRUE, legend="bottom")
ggsave(locationplot, file="2_Diversity/2b_beta_diversity.by_location.jgw.png", width=8, height=4)


# Still by location, but split with facets by region
plots.split = lapply(names(plotdata), function(metric){
  mydata=plotdata[[metric]]
  nloc = length(unique(mydata$plotdata$location)) # number of locations, for coloring
  ggplot(mydata$plotdata) +
    aes(x = PC1, y = PC2, color = location) +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = mydata$axis_labels[1], y=mydata$axis_labels[2],
         title=metric, color="Location") +
    facet_wrap(~region, ncol=1) +
    scale_color_manual(values=colors.location )
})

splitplot = ggarrange(plotlist=plots.split, nrow=1, common.legend=TRUE, legend="bottom")
ggsave(splitplot, file="2_Diversity/2b_beta_diversity.split.jgw.png", width=8, height=10)


