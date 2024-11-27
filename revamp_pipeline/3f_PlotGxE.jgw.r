#! /usr/bin/Rscript

# Plot the various GxE metrics for publication

library(tidyverse)
library(ggpubr)

# Load data
alpha = read.csv("3_GxE/3b_alpha_diversity.variance_explained.jgw.csv")
beta = read.csv("3_GxE/3c_beta_diversity_GxE.jgw.csv")
pathways = read.csv("3_GxE/3d_MetaCyc_pathway_GXE.jgw.csv")
# taxonomy = TODO

# Global variables
sig_threshold=0.05 # Significance threshold
colorkey = c("Maize Genotype" = "#74afda", "GXE" = "#00a693", "Environment" = "#da9100")

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
    theme(legend.position = "bottom", 
          axis.title.x = element_blank(), 
          legend.title = element_blank(), 
          axis.text.x = element_text(angle=90)) +
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
  summarize(fraction = sum(fdr < 0.01)/n()) %>%
  mutate(percent = round(fraction*100, digits=1) %>% paste("% sig."))

# Make violin plot  
pathplot = ggplot(pathways) +
  aes(x = term, y = herit, fill = term) +
  geom_violin(trim = TRUE, linewidth=0.25) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.2, color = "black") +
  geom_text(data=path.significant, mapping=aes(label=percent), y=max(pathways$herit),
            size=3) +
  ylab("Variance Explained") +
  ggtitle(paste("Metabolism (inferred)")) + 
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values=colorkey)

# Combine plots and save
top = ggarrange(alphaplot, betaplot, nrow=1, labels=c("A","B"))
bottom = ggarrange(pathplot, nrow=1, ncol=2, labels=c("C","D"))
combined = ggarrange(top, bottom, nrow=2, common.legend=TRUE)
ggsave(combined, file="3_GxE/3f_GxE_combined_plot.jgw.png", width = 6, height = 8)
