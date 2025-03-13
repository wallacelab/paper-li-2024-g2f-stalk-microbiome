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

# Write out distance matrices and PCs for later
for(metric in names(betas)){
  distance_file = paste("2_Diversity/2b_beta_diversity.", metric, ".distance.csv", sep="") %>%
    gsub(pattern=" ", repl="_")
  pc_file = sub(distance_file, pattern=".distance.", repl=".pcs.", fixed=T)
  write.csv(as.matrix(betas[[metric]]), file=distance_file)
  write.csv(pcs[[metric]]$points, file=pc_file)
}



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



##############
# Plot - For publication
##############


# Overall plot
metric='Weighted Unifrac'
regions = unique(mydata$plotdata$region) %>% sort()
mydata=plotdata[[metric]]
pub.overall = ggplot(mydata$plotdata) +
  aes(x = PC1, y = PC2, color = region) +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = mydata$axis_labels[1], y=mydata$axis_labels[2],
       title=metric, color="Region") +
  scale_color_manual(values=colors.region ) +
  theme(legend.position="bottom")
  #theme(legend.position=c(0.83,0.15), 
   #    legend.background = element_rect(fill="#FFFFFFAA"))
# # Add custom legend-like annotation
# for(i in 1:length(regions)){
#   pub.overall = pub.overall + 
#     grid::textGrob(regions[i], x=1, y=0.5, 
#                    gp = gpar(color=regionkey$region_color[regionkey$region==regions[i]]))
#     #annotate("text", x = Inf, y = Inf, label = regions[i], vjust=1, hjust=0,
#      #        color=regionkey$region_color[regionkey$region==regions[i]])
# }
  

# Individual Regions
regionals = lapply(regions, function(myregion){
  ggplot(mydata$plotdata %>% filter(region==myregion)) +
    aes(x = PC1, y = PC2, color = location) +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = "PC1", y="PC2",
         title=myregion, color="Location") +
    scale_color_manual(values=colors.location ) +
    theme(legend.position = "bottom")
})
pub.regions = ggarrange(plotlist = regionals, nrow=2, ncol=2)

# Combine all together
pubplot = ggarrange(pub.overall, pub.regions, nrow=1, widths=c(1,1.5))
ggsave(pubplot, file="2_Diversity/2b_beta_diversity.publication.jgw.png", width=12, height=7)



