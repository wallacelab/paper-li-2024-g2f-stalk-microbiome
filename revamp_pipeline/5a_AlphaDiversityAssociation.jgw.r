#! /usr/bin/Rscript

# Associate environmental parameters with alpha diversity
# NOTE: WDRF is apparently a shortening of "Woodruff", who developed the specific buffer pH 
#       method used. It measures "total acidity" instead of "active acidity"
#       (see https://www.no-tillfarmer.com/ext/resources/files/nntc/2020/NTF_eBook_Healthier-Soils_2020-002.pdf)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(gridExtra)
library(phyloseq)
library(tidyverse)

# Load data
alpha_diversity = read.csv("2_Diversity/2a_alpha_diversity.jgw.csv")

###########################################################

Pvaluecal <- function(modelsummary) {
  f_value <- modelsummary$fstatistic["value"]
  numdf <- modelsummary$fstatistic["numdf"]
  dendf <- modelsummary$fstatistic["dendf"]
  p <- pf(f_value, numdf, dendf, lower.tail = FALSE)
  return(p)
}

###############################################################



# Investigate the median alpha diversity in a location vs environmental factors
G2F_2019_median_alpha_diversity_by_location_data <- alpha_diversity %>%
  group_by(location) %>%
  summarize(
    median_shannon = median(Shannon, na.rm = TRUE),
    median_observed_features = median(Observed, na.rm = TRUE),
    median_simpson = median(Simpson, na.rm = TRUE)
  )

G2f_2019_soil_data <- read.csv("0_data_files/g2f_2019_soil_data.csv", sep = ",")
G2f_2019_weather_data <- read.csv("0_data_files//G2F_2019_weather_average_and_sum.tsv", sep = "\t")

# Joining soil and weather data with median alpha diversity data
G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data <- G2F_2019_median_alpha_diversity_by_location_data %>%
  left_join(G2f_2019_soil_data, by = c("location" = "Location")) %>%
  left_join(G2f_2019_weather_data, by = "location") %>%
  filter(!is.na(X1.1.Soil.pH)) %>%
  as.data.frame()

#######
# Test associations
#######

# DEPRECATED - This tests everything by everything; resulting models are overparameterized (and useless)
# # Create linear models for each environmental factor against each diversity index
# create_lm <- function(df, formula) {
#   lm_results <- lapply(df, function(x) lm(formula, data = df))
#   lm_summaries <- lapply(lm_results, summary)
#   return(lm_summaries)
# }
# 
# shannon_result_summary <- create_lm(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data, median_shannon ~ .)
# observed_features_result_summary <- create_lm(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data, median_observed_features ~ .)
# simpson_result_summary <- create_lm(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data, median_simpson ~ .)

# Get metrics to test
metrics = names(G2F_2019_median_alpha_diversity_by_location_data) %>% setdiff("location")
# Get environmental covariates to test
covariates = names(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data) %>%
  setdiff(c("location", metrics))  # Remove 'location' and metrics to test

# Test
models=list()
summaries=list()
anovas=list()
pvals=list()
for(metric in metrics) {
  for(covariate in covariates){
    key = paste(metric, covariate, sep="$")
    myformula = paste(metric, "~", covariate)
    mymodel = lm(myformula, data=G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data)
    models[[key]] = mymodel
    summaries[[key]] = summary(mymodel)
    anovas[[key]] = anova(mymodel)
    pvals[[key]] = anovas[[key]]$`Pr(>F)`[1]  # Top p-value from ANOVA table, which is the model term
  }
}

# Extract p-value from ANOVA and collate
pvals.long = lapply(names(pvals), function(myname){
  key = strsplit(myname, split="$", fixed=TRUE)[[1]] # Only 1 item to split, so extract from list
  p = pvals[[myname]]
  result = data.frame(metric=key[1], covariate=key[2], pvalue=p)
}) %>% bind_rows() %>%
  arrange(pvalue)

# Save
write.csv(pvals.long, file="5_Associations/5a_median_alpha_diversity_pvals.jgw.csv", row.names=FALSE)


#########
# Plot median alpha diversity versus pH
##########

# Function to plot median diversity index against soil pH with regression line
plot_diversity_vs_soil_pH <- function(df, x_var, y_var, xlab, ylab) {
  plotdata = data.frame(x=df[,x_var], y=df[,y_var])
  pvalue = pvals.long$pvalue[pvals.long$covariate==x_var & pvals.long$metric==y_var] %>%
    signif(digits=3)
  ggplot(plotdata) +
    aes(x = x, y = y) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
    #stat_cor(size=14, cor.coef.name="R", label.y = max(plotdata$y) * 1.025) + 
    annotate("text" ,  x=-Inf, y = -Inf, label = paste("p =",pvalue), vjust=-1, hjust=-0.1,
             size=12) +
    labs(x = xlab, y = ylab) +
    theme(
      legend.text = element_text(color = "black", size = 30),
      legend.title = element_text(color = "black", size = 30),
      axis.title = element_text(size=30, face="bold"),
      plot.title = element_text(size=30, face="bold"),
      axis.text = element_text(size=30, face="bold"),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    )
}


# Create plots for each diversity index
labellookup = c("median_shannon" = "Median Shannon Index",
                "median_observed_features" = "Median Observed Features",
                "median_simpson" = "Median Simpson Index")
# Soil pH - original
soilplots = lapply(metrics, function(mymetric){
  mylabel = labellookup[mymetric]
  plot_diversity_vs_soil_pH(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data,
                            "X1.1.Soil.pH", mymetric, xlab = "Soil pH", ylab = mylabel)
})
soilplots = ggarrange(plotlist = soilplots, nrow=1, labels=LETTERS[1:3], 
                      font.label=list(size=36, face="bold"), vjust=1)
ggsave(soilplots, file="5_Associations/5a_median_alpha_diversity_plot.soil_ph.jgw.png", height = 8, width = 24, device = "png")

# Soil buffering capacity
bufferplots = lapply(metrics, function(mymetric){
  mylabel = labellookup[mymetric]
  plot_diversity_vs_soil_pH(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data,
                            "WDRF.Buffer.pH", mymetric, xlab = "WDRF buffer pH", ylab=mylabel)
})
bufferplots = ggarrange(plotlist = bufferplots, nrow=1, labels=LETTERS[4:6], 
                        font.label=list(size=36, face="bold"), vjust=1)
ggsave(bufferplots, file="5_Associations/5a_median_alpha_diversity_plot.buffer_ph.jgw.png", height = 8, width = 24, device = "png")


