# Load necessary libraries for beta diversity analysis
library(tidyverse)
library(phyloseq)
library(ggpubr)

# Beta Diversity Analysis Using Weighted and Unweighted UniFrac Distances

background_color="gray96"
min_geno_count=20 # Minimum times a genotype has to appear to be plotted

#####
# Load PC data
#####

# Location key & metadata
lockey = read.csv("1_parsed_files/1b_location_key.csv") 
asvs = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")
metadata = sample_data(asvs) %>% data.frame()

# PCs
pcs = readRDS("2_Diversity/2b_beta_diversity.pcs.rds")


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
ggsave(regionplot, file="2_Diversity/2d_beta_diversity.by_region.png", width=8, height=3)


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
ggsave(locationplot, file="2_Diversity/2d_beta_diversity.by_location.png", width=8, height=4)


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
ggsave(splitplot, file="2_Diversity/2d_beta_diversity.split.png", width=8, height=10)



##############
# Plot - For publication - Weighted by location
##############


# Overall plot
metric='Weighted Unifrac'
mydata=plotdata[[metric]]
regions = unique(mydata$plotdata$region) %>% sort()
pub.overall = ggplot(mydata$plotdata) +
  aes(x = PC1, y = PC2, color = region) +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = mydata$axis_labels[1], y=mydata$axis_labels[2],
       title=metric, color="Region") +
  scale_color_manual(values=colors.region ) +
  theme(legend.position="bottom", 
        panel.background = element_rect(fill=background_color))

# Individual Regions
regionals = lapply(regions, function(myregion){
  ggplot(mydata$plotdata %>% filter(region==myregion)) +
    aes(x = PC1, y = PC2, color = location) +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = "PC1", y="PC2",
         title=myregion, color="Location") +
    scale_color_manual(values=colors.location ) +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill=background_color))
})
pub.regions = ggarrange(plotlist = regionals, nrow=2, ncol=2)

# Combine all together
pubplot = ggarrange(pub.overall, pub.regions, nrow=1, widths=c(1,1.5))
ggsave(pubplot, file="2_Diversity/2d_beta_diversity.weighted.publication.png", width=12, height=7)



##############
# Supplemental - Unweighted by location
##############

# Overall plot
metric='Unweighted Unifrac'
mydata=plotdata[[metric]]
regions = unique(mydata$plotdata$region) %>% sort()
pub.overall = ggplot(mydata$plotdata) +
  aes(x = PC1, y = PC2, color = region) +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = mydata$axis_labels[1], y=mydata$axis_labels[2],
       title=metric, color="Region") +
  scale_color_manual(values=colors.region ) +
  theme(legend.position="bottom",
        panel.background = element_rect(fill=background_color))

# Individual Regions
regionals = lapply(regions, function(myregion){
  ggplot(mydata$plotdata %>% filter(region==myregion)) +
    aes(x = PC1, y = PC2, color = location) +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = "PC1", y="PC2",
         title=myregion, color="Location") +
    scale_color_manual(values=colors.location ) +
    theme(legend.position = "bottom", 
          panel.background = element_rect(fill=background_color))
})
pub.regions = ggarrange(plotlist = regionals, nrow=2, ncol=2)

# Combine all together
pubplot = ggarrange(pub.overall, pub.regions, nrow=1, widths=c(1,1.5))
ggsave(pubplot, file="2_Diversity/2d_beta_diversity.unweighted.supplemental.png", width=12, height=7)




##############
# Supplemental - Weighted & Unweighted by genotype
##############

plot_by_geno = function(plotdata, metric, min_count=1){
  # Select and filter data for genotypes present minimum count of times
  mydata=plotdata[[metric]]
  toplot = mydata$plotdata %>%
    group_by(Corrected_pedigree) %>%
    filter(n() >= min_count)
  # Make plot
  pub.overall = ggplot(toplot) +
    aes(x = PC1, y = PC2, color = Corrected_pedigree) + 
    geom_point(size = 2, alpha = 0.8) +
    labs(x = mydata$axis_labels[1], y=mydata$axis_labels[2],
         title=metric, color="Genotype",
         subtitle=paste("Only showing genotypes with ≥", min_count,"samples")) +
    theme(legend.position="bottom",
          panel.background = element_rect(fill=background_color),
          legend.key.height=unit(0.15, "cm"),
          plot.subtitle = element_text(size=8))
  return(pub.overall)
}

# Make plots
geno.weight = plot_by_geno(plotdata, "Weighted Unifrac", min_count=min_geno_count)
geno.unweight = plot_by_geno(plotdata, "Unweighted Unifrac", min_count=min_geno_count)

# Arrange and write out
genoplots = ggarrange(geno.weight, geno.unweight, nrow=1, 
                      common.legend=TRUE, legend="right")
ggsave(genoplots, file="2_Diversity/2d_beta_diversity.genotypes.supplemental.png", width=8, height=4)


# ##############
# # Supplemental - tSNE -> Didn't give any great insight, so didn't use
# ##############
# 
# metric = 'Weighted Unifrac'
# mydist = betas[[metric]]
# 
# tsne = Rtsne(mydist, is_distance=TRUE, dims=2, perplexity=30, max_iter=2000, verbose=TRUE)
# plotdata=data.frame(t1=tsne$Y[,1], t2=tsne$Y[,2], 
#                     sample=rownames(as.matrix(mydist)))
# plotdata$location = metadata[plotdata$sample, "location"]
# plotdata = left_join(plotdata, lockey, by="location") 
# 
# tplot=ggplot(plotdata) +
#   aes(x = t1, y = t2, color = region) +
#   geom_point(size = 2, alpha = 0.8) +
#   labs(x = "tSNE-1", y="tSNE-2",
#        title=metric, color="Region") +
#   scale_color_manual(values=colors.region ) +
#   theme(legend.position="bottom")
# 
# # Individual Regions
# t.regionals = lapply(regions, function(myregion){
#   ggplot(plotdata %>% filter(region==myregion)) +
#     aes(x = t1, y = t2, color = location) +
#     geom_point(size = 2, alpha = 0.8) +
#     labs(x = "PC1", y="PC2",
#          title=myregion, color="Location") +
#     scale_color_manual(values=colors.location ) +
#     theme(legend.position = "bottom")
# })
# t.regions = ggarrange(plotlist = t.regionals, nrow=2, ncol=2)
