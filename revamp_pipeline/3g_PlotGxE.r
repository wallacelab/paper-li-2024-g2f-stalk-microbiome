#! /usr/bin/Rscript

# Plot the various GxE metrics for publication

library(tidyverse)
library(ggpubr)

# Load data
alpha = read.csv("3_GxE/3b_alpha_diversity.variance_explained.jgw.csv")
beta = read.csv("3_GxE/3c_beta_diversity_GxE.jgw.csv")
pathways = read.csv("3_GxE/3e_MetaCyc_pathway_GXE.csv")
taxonomy = read.csv("3_GxE/3f_taxon_GXE.csv")

# Global variables
sig_threshold=0.05 # Significance threshold for alpha & beta diversity
fdr_cutoff.pathways=0.01 # FDR cutoff for GxE pathways
colorkey = c("Maize Genotype" = "#74afda", "GXE" = "#00a693", "Environment" = "#da9100")
taxa_levels=c("Phylum", "Class","Order","Family","Genus","Species")

#########
# Alpha and Beta Diversity
#########

# Standardize tables
alpha = alpha %>%
  rename(metric=alpha_diversity, pvalue=p, fraction_explained=Value) %>%
  select(category, metric, test_var, fraction_explained, pvalue) %>%
  filter(test_var != "Residuals")
beta = beta %>%
  rename(metric=Type, pvalue = Pr..F., fraction_explained=heritability)%>%
  select(category, metric, test_var, fraction_explained, pvalue) %>%
  filter(test_var != "Residuals")

# Helper function for alpha and beta diversity plots
div_barplot = function(plotdata, set="???"){
  plotdata = plotdata %>%
    mutate(alpha=ifelse(pvalue < sig_threshold, yes=1, no=0.75),
           label=ifelse(pvalue < sig_threshold, yes="", no="n.s."))
  ggplot(plotdata) +
    aes(x = metric, y = fraction_explained, fill = test_var, alpha=alpha, label=label) +
    geom_bar(position = "stack", stat = "identity") +
    ylab("Variance Explained") +
    ggtitle(paste(set, "Diversity")) + 
    theme_minimal() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), 
          legend.title = element_blank(), 
          axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_alpha(guide='none', range=c(0.3,1)) + # Use to set how transparent becomes
    scale_fill_manual(values=colorkey) + 
    geom_text(alpha=1, position=position_stack(vjust=0.5), size=2)
}

# Alpha  & Beta plots
alphaplot = div_barplot(alpha, set="Alpha")
betaplot = div_barplot(beta, set="Beta")

##########
# Pathways
##########

# Significance
path.significant = pathways %>%
  group_by(term) %>%
  summarize(n=n(),
            n_significant = sum(fdr < fdr_cutoff.pathways),
            fraction = n_significant/n) %>%
  mutate(percent = round(fraction*100, digits=1) %>% paste("%", sep=""),
         n_significant = paste("n =", n_significant))

# Make violin plot  
pathplot = ggplot(pathways) +
  aes(x = term, y = herit, fill = term) +
  geom_violin(trim = TRUE, linewidth=0.25) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "black") +
  geom_text(data=path.significant, mapping=aes(label=percent), 
            size=2.5, y=max(pathways$herit)* 1.06) +
  geom_text(data=path.significant, mapping=aes(label=n), 
            size=2, y=max(pathways$herit * 1.02)) +
  ylab("Variance Explained") +
  ggtitle(paste("Metabolism (inferred)")) + 
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor=element_blank()) +
  scale_fill_manual(values=colorkey) +
  ylim(NA, max(pathways$herit) * 1.1) # Adjust vertical limit

# geom_text(data=sig.taxa, mapping=aes(label=percent), 
#           size=2.5, y=max(taxonomy$herit)* 1.06) +
#   geom_text(data=sig.taxa, mapping=aes(label=n), 
#             size=2, y=max(taxonomy$herit) * 1.02) +



##############
# Taxonomy
##############


# Calculate fractions that are significant
sig.taxa = taxonomy %>%
  group_by(level, term) %>%
  summarize(n=n(), fraction = sum(pval < 0.01)/n(), .groups="drop") %>%
  mutate(percent = round(fraction*100, digits=1) %>% paste("%", sep=""))


# Plot
taxonomy$level = factor(taxonomy$level, levels=taxa_levels)
taxaplot <- ggplot(taxonomy) +
  aes(x=level, y=herit, fill=term) +
  geom_violin(trim = TRUE, linewidth=0.25) +
  ylab("Variance Explained") + 
  ggtitle(paste("Taxonomy")) + 
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(~term) +
  # theme(strip.text.x = element_text(size = 30,face="bold"),
  #       strip.text.y = element_text(size = 30,face="bold" )) + 
  stat_summary(fun = "mean", geom = "crossbar",  width = 0.5, colour = "black",
               show.legend=FALSE) +
  geom_text(data=sig.taxa, mapping=aes(label=percent), 
            size=2.5, y=max(taxonomy$herit)* 1.06) +
  geom_text(data=sig.taxa, mapping=aes(label=n), 
            size=2, y=max(taxonomy$herit) * 1.02) +
  ylim(NA, max(taxonomy$herit) * 1.1) + # Adjust vertical limit
  scale_fill_manual(values=colorkey)



########################
# Combine plots and save
########################

top = ggarrange(alphaplot, betaplot, pathplot, nrow=1, labels=c("A","B", "C"),
                align="h")
bottom = ggarrange(taxaplot, labels=c("D"))
combined = ggarrange(top, bottom, nrow=2, common.legend=TRUE)
ggsave(combined, file="3_GxE/3g_GxE_combined_plot.png", width = 7, height = 8)


