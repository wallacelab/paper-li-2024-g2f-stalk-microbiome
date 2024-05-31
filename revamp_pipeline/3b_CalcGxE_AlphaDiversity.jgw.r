#! /usr/bin/Rscript

# Calculate GxE for alpha diversity metrics
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(lmerTest) # TODO - This ever used? I don't think so
library(EnvStats)
library(gvlma)

# Parameters
infile="3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds" # Input (filted) data file for GxE analysis
alphafile="2_Diversity/2a_alpha_diversity.jgw.csv"
alpha_id_cols=c("sampleID", "location", "rep_number", "plot_number", "Corrected_pedigree")  # Columns in above file that are just for IDs

#############
# Load data
#############

mydata = readRDS(infile)
alphadiv = read.csv(alphafile)

############
# Merge sample data with alpha diversity
############

metadata = sample_data(mydata) %>%
  data.frame() %>%
  rownames_to_column("sampleID") 
joindata = right_join(alphadiv, metadata)

# Determine the names of the provided alpha diversity metrics
alpha_names = setdiff(names(alphadiv), alpha_id_cols)


##########
# Plot alpha diversity of just these samples
##########

# Visualize the alpha diversity (& check for normality)
plot1 = ggplot(data = joindata, aes(x=Shannon)) + geom_histogram()
plot2 = ggplot(data = joindata, aes(x=Observed)) + geom_histogram()
plot3 = ggplot(data = joindata, aes(x=Simpson)) + geom_histogram()
alphaplots = ggarrange(plotlist=list(plot1, plot2, plot3), nrow=1)
ggsave(alphaplots, file="3_GxE/3b_alpha_diversity_normality_check.jgw.png", width=10, height=4)

#########
# Transform non-normal metrics
#########

# Reflect and log-transform Simpson index
joindata$Simpson.reflect_log = log(1-joindata$Simpson)
alpha_names = c(alpha_names, "Simpson.reflect_log")

# Log-transform Observed ASVs
joindata$Observed.log = log(joindata$Observed+1)
alpha_names = c(alpha_names, "Observed.log")

#########
# Check for outliers with Rosner's test
#########

# Check for max of 10 outliers
rosners = lapply(alpha_names, function(myalpha){
  rosnerTest(joindata[,myalpha], k = 10)
})

outliers = lapply(rosners, function(r){
  r$all.stats$Obs.Num[r$all.stats$Outlier]
})
names(outliers) = alpha_names

# Remove outliers
for(myalpha in alpha_names){
  joindata[outliers[[myalpha]],myalpha] = NA
}

# Visualize the alpha diversity (& check for normality)
postplots = lapply(alpha_names, function(myalpha){
  ggplot(joindata) + 
    aes(x=.data[[myalpha]]) +
    geom_histogram(fill="darkblue", bins=30)
})
alphaplots = ggarrange(plotlist=postplots, nrow=1)
ggsave(alphaplots, file="3_GxE/3b_alpha_diversity_normality_check.postfilter.jgw.png", width=3*length(postplots), height=4)


###########
# Calculate GxE for alpha diversity
###########

# Basic linear model
models = lapply(alpha_names, function(myalpha){
  myformula = paste(myalpha, "~", "location + Corrected_pedigree + location:Corrected_pedigree")
  lm(formula = myformula, data = joindata)
})
names(models) = alpha_names

# Basic type I anova
anovas = lapply(models, anova)
names(anovas) = alpha_names # So can pull them out specifically in the next step

# Calculate % variation explained
variance_explained = lapply(alpha_names, function(alpha_name){
  # Pull out specific ANOVA result
  anova_result = anovas[[alpha_name]]
  
  # Reformat anova results
  result_df = anova_result %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(alpha_diversity = alpha_name, 
           Value = `Sum Sq`/ sum(`Sum Sq`)) %>%
    rename(p=`Pr(>F)`)
  
  # Rename terms so are in G, E, and GxE
  result_df$test_var = sapply(result_df$term, switch, 
                              location="Environment",
                              Corrected_pedigree="Maize Genotype",
                              "location:Corrected_pedigree"="GXE",
                              "Residuals"="Residuals",
                              "unknown")
  
  # Change data frame order for convenience
  result_df = result_df %>%
    relocate(alpha_diversity, test_var, Value, p)
  
  return(result_df)
})

# Bind into a single data frame
alpha_result = bind_rows(variance_explained) %>%
  mutate(category="Alpha Diversity") %>%
  relocate(category)

##############
# Check Linear Regression Assumptions
##############

# Test with GVLMA; have to do some tricks to get rid of interactions that don't appear in the data
assumptions = lapply(alpha_names, function(myalpha){
  mymodel = models[[myalpha]]
  alphavals = mymodel$model[[myalpha]] # Extract Y values from model
  mymatrix = model.matrix(mymodel) # Extract model matrix
  mymatrix.trimmed = mymatrix[,!is.na(mymodel$coefficients)] # Only take non-missing coefficents
  trimmodel = lm(alphavals ~ mymatrix.trimmed + 0) # Fit new full-column-rank model
  
  # Test assumptions and format for printing out
  gvlma(trimmodel) %>%
    display.gvlmatests() %>% 
    as.data.frame() %>%
    rownames_to_column("gvlma_test") %>%
    mutate(metric=myalpha) %>%
    relocate(metric)
}) %>% bind_rows() %>%
  mutate(pass= Decision == "Assumptions acceptable.")

# Plot linear regression assumptions
testplot = ggplot(assumptions) +
  aes(x=metric, y=gvlma_test, label=round(`p-value`, digits=4), fill=pass) +
  geom_tile() +
  geom_text() +
  labs(x="Diversity Metric", y="GVLMA test", title="Test regression assumptions (p-values)") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values=c("FALSE"="lightpink", "TRUE"="dodgerblue")) +
ggsave(testplot, file="3_GxE/3b_alpha_diversity.gvlma_tests.jgw.png", height=6, width=6)

# Write out GVLMA checks
write.csv(assumptions, file="3_GxE/3b_alpha_diversity.gvlma_tests.jgw.csv", row.names=FALSE)

##############
# Output data and plots of variance explained by each term
##############

plotdata = alpha_result %>%
  filter(test_var != "Residuals")

#plot the result in a bar graph 
alphaplot = ggplot(data=plotdata) +
  aes(x=alpha_diversity, y=Value, fill=test_var) +
  geom_col(position="stack") + 
  labs(y = "variance explained", fill = "category") + 
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),limits = c(0,0.8)) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey20"), 
        panel.background = element_blank(),
        legend.title = element_text(size = 15, face = "bold", colour = "grey30"),) + 
  theme(
    legend.text = element_text(color = "black", size = 15),
    axis.title.y = element_text(size=15, face="bold"),
    plot.title = element_text(size=15, face="bold"),
    axis.text = element_text(size=15, face="bold"),
    axis.text.x = element_text(size=10, angle=90, vjust=0.5, hjust=1),
    strip.text.x = element_text(size=15,face="bold"),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    legend.position="bottom",
  )

#output the graph
ggsave(alphaplot, file="3_GxE/3b_alpha_diversity.variance_explained.jgw.png", height=6, width=6, device="png")

# Write data
write.csv(alpha_result, file="3_GxE/3b_alpha_diversity.variance_explained.jgw.csv", row.names=FALSE)

