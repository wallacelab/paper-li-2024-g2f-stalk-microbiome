#! /usr/bin/Rscript

# Look for associations between beta diversity and environmental parameters

# TODO - Make PC plots and color by the indicated variables to see if can see things visually
#    - Current setup doesn't work well. Need to separate PC plots to get separate color scales

library(dplyr)
library(ggplot2)
library(gridExtra)
library(vegan) # for mantel test
library(rbiom)
library(phyloseq)

# Variables to ignore from sampled ata
to_ignore = c("location", "pedigree", "rep_number", "plot_number", "Corrected_pedigree",
              "Sub", "Grower", "Day.before.extract", "Date.Received","Date.Reported")
ncores=7 # Parallel cores to use
plot_p_threshold=0.05 # Cutoff for whether a mantel test is significant enough to plot

# Load data - Less filtered set
asvs = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")

###############
# Preprocessing & distance calculations
###############

# Merge by location
better_mean = function(x){mean(x, na.rm=TRUE)}
asv_locs = merge_samples(asvs, group="location", fun=better_mean)
asv_locs.rarefy = rarefy_even_depth(asv_locs, rngseed=1, replace=FALSE)

# Calculate various Beta diversity distances
counts = t(otu_table(asv_locs.rarefy))  # Apparentlt flipped rows/cols
mytree = phy_tree(asv_locs.rarefy)
betas = list()  #List to hold results
betas[['Weighted Unifrac']] = rbiom::beta.div(counts, method="unifrac", weighted=TRUE, tree=mytree)
betas[['Unweighted Unifrac']] = rbiom::beta.div(counts, method="unifrac", weighted=FALSE, tree=mytree)
betas[['Bray-Curtis']] = rbiom::beta.div(counts, method="bray-curtis")

# Calculate environmental distances
targets = names(sample_data(asv_locs.rarefy)) %>%
  setdiff(to_ignore)
env_dists = lapply(targets, function(mytarget){
  # Subset environmental data
  subdata = sample_data(asv_locs.rarefy)[,mytarget]
  if(all(is.na(subdata[,mytarget]))){ # Check if any data left
    return(NULL)
  }
  subdata = subset(subdata, !is.na(subdata[,1]))

  #Get distances
  dist(subdata)
})
names(env_dists) = targets
env_dists = env_dists[!sapply(env_dists, is.null)] # Remove any null results

##################
# Mantel tests
#################

mantels=list()
for(metric in names(betas)){
  for(target in names(env_dists)){
    key = paste(metric, target, sep="|")
    dm1 = betas[[metric]] %>% as.matrix()
    dm2 = env_dists[[target]]  %>% as.matrix()
    # Standardize matrices
    locs = intersect(rownames(dm1), rownames(dm2))
    dm1 = dm1[locs, locs]
    dm2 = dm2[locs, locs]
    # Run mantel and save results
    mymantel = vegan::mantel(dm1, dm2, method = "pearson", permutations = 999, parallel=ncores)
    result = data.frame(beta=metric, env=target, size=length(locs), statistic=mymantel$statistic, p=mymantel$signif)
    mantels[[key]] = result
  }
}
mantels = bind_rows(mantels)

# Sort & write out
mantels = mantels %>% arrange(p)
write.csv(mantels, file="5_Associations/5b_mantel_tests.jgw.csv", row.names=FALSE)

################
# Plot significant results
################

# Subset to significant associations
significant = subset(mantels, mantels$p <= plot_p_threshold) %>%
  mutate(key = paste(beta, env, sep="|"),
         label = paste("r=", round(statistic, digits=3), "; p=",p, sep=""))

# Make individual plots
plots=list()
for(i in 1:nrow(significant)){
  metric=significant$beta[i]
  target=significant$env[i]
  # Recycling some code from above to standardize matrices
  dm1 = betas[[metric]] %>% as.matrix()
  dm2 = env_dists[[target]]  %>% as.matrix()
  locs = intersect(rownames(dm1), rownames(dm2))
  dm1 = dm1[locs, locs]
  dm2 = dm2[locs, locs]
  # Extract results
  results = data.frame(x=dm2[lower.tri(dm2)], y=dm1[lower.tri(dm1)])
  plots[[i]] = ggplot(results) +
    aes(x=x, y=y) +
    geom_point() + 
    geom_smooth(method="lm") +
    labs(x=target, y=metric, subtitle=significant$label[i])
}

# Combine and write out
allplots = grid.arrange(grobs=plots, ncol=3)
ggsave(allplots, file="5_Associations/5b_mantel_tests.jgw.png", width=8, height=8)


# Second try - Make longer data frame to use a single ggplot
plots=list()
get_plotdata = function(dm, locs, valname="unknown"){
  dm = dm[locs, locs] %>%
    as.data.frame() %>%
    rownames_to_column("loc1") %>%
    pivot_longer(cols=-loc1, values_to=valname) %>%
    rename(loc2=name) %>%
    filter(loc1 != loc2) # Remove self-plot
}
for(i in 1:nrow(significant)){
  metric=significant$beta[i]
  target=significant$env[i]
  # Recycling some code from above to standardize matrices
  dm.beta = betas[[metric]] %>% as.matrix()
  dm.env = env_dists[[target]]  %>% as.matrix()
  locs = intersect(rownames(dm.beta), rownames(dm.env))
  # Extract results
  beta.plotdata = get_plotdata(dm.beta, locs, "betadist") %>%
    mutate(beta=metric)
  env.plotdata = get_plotdata(dm.env, locs, "envdist") %>%
    mutate(env=target)
  plotdata = left_join(beta.plotdata, env.plotdata, by=c("loc1", "loc2"))
  plots[[i]] = plotdata
}
plotdata = bind_rows(plots)

# Plot with grid arrange
plot2 = ggplot(plotdata) + 
  aes(x=envdist, y=betadist, color=env) +
  geom_point() + 
  geom_smooth(method="lm", color="black") +
  facet_grid(beta ~ env, scales="free", switch="x") +
  theme(legend.position="none") +
  geom_text(mapping=aes(label=label), x=-Inf, y=Inf, data=significant,
            color="black", vjust=2, hjust=-0.1)

ggsave(plot2, file="5_Associations/5b_mantel_tests.jgw_alt.png", width=10, height=8)

###############
# Principal Components
###############

# Set up again for plotting
combinedata=list()
for(i in 1:nrow(significant)){
  metric=significant$beta[i]
  target=significant$env[i]
  
  # Get PCs
  pcs = cmdscale(betas[[metric]], k=2) %>%
    as.data.frame()
  names(pcs) = c("PC1", "PC2")
  
  # Get environmental data
  env = sample_data(asv_locs.rarefy)[,target] %>%
    data.frame() %>%
    rownames_to_column("location") %>%
    pivot_longer(cols=-location, names_to="env_metric")
  # Combine with environmental data
  combined = pcs %>%
    mutate(beta_metric=metric) %>%
    rownames_to_column("location") %>%
    left_join(env, by="location")
  combinedata[[i]] = combined
}
plotdata = bind_rows(combinedata)



# Plot PCs with grid arrange
pca_plot = ggplot(plotdata) + 
  aes(x=PC1, y=PC2, color=value) +
  geom_point() + 
  facet_grid(beta_metric ~ env_metric, scales="free", switch="x") # +
  #theme(legend.position="none") +
  #geom_text(mapping=aes(label=label), x=-Inf, y=Inf, data=significant,
  #          color="black", vjust=2, hjust=-0.1)

ggsave(pca_plot, file="5_Associations/5b_mantel_tests.jgw_alt_pca.png", width=10, height=8)

