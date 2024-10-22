#! /usr/bin/Rscript

# Plot the various GxE metrics for publication

library(tidyverse)
library(ggpubr)

# Load data
alpha = read.csv("3_GxE/3b_alpha_diversity.variance_explained.jgw.csv")
beta = read.csv("3_GxE/3c_beta_diversity_GxE.jgw.csv")
# pathways = TODO
# taxonomy = TODO

# Global variables
sig_threshold=0.05 # Significance threshold
colorkey = c("Maize Genotype" = "#74afda", "GXE" = "#00a693", "Environment" = "#da9100")

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

# Alpha plot
alphaplot = div_barplot(alpha, set="Alpha")
betaplot = div_barplot(beta, set="Beta")


# Combine plots and save
top = ggarrange(alphaplot, betaplot, NULL, nrow=1)
bottom = ggarrange(NULL)
combined = ggarrange(top, bottom, nrow=2, common.legend=TRUE)
ggsave(combined, file="3_GxE/3f_GxE_combined_plot.jgw.png", width = 8, height = 6)
