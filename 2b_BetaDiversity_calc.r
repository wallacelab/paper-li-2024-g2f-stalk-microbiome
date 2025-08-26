# Load necessary libraries for beta diversity analysis
library(phyloseq)
library(tidyverse)
library(rbiom)
library(qiime2R) # Provides the read_qza function
library(ggpubr)
library(Rtsne) # For tSNE plot

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

# write out pc data as RDS for easier plotting later
write_rds(pcs, file="2_Diversity/2b_beta_diversity.pcs.rds")


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
