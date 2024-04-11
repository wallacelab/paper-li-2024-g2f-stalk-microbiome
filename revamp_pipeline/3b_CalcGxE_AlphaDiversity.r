#! /usr/bin/Rscript

# Calculate GxE for alpha diversity metrics
library(phyloseq)
library(tidyverse)
library(gridExtra)
library(lmerTest) # TODO - This ever used? I don't think so
library(EnvStats)
library(gvlma)

# TODO: Is it better to merge prior alpha diversity, or calculate it fresh? Does it matter?
# TODO: Observed and Simpson are not normal. Make sure code accounts for that.
#       TODO: Code just removes outliers. Jason thinks may need to be transformed instead
# TODO: Need to do GVLMA assumption checks for each model. Do with lapply?

# Parameters
infile="3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds" # Input (filted) data file for GxE analysis
alphafile="2_Diversity/2a_alpha_diversity.csv"
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

# # TODO: Are the below necessary? Commented out for now
# joindata$location  <- as.factor(joindata$location)
# joindata$Corrected_pedigree <- as.factor(joindata$Corrected_pedigree)

# Determine the names of the provided alpha diversity metrics
alpha_names = setdiff(names(alphadiv), alpha_id_cols)


##########
# Plot alpha diversity of just these samples
##########

# Visualize the alpha diversity (& check for normality)
plot1 = ggplot(data = joindata, aes(x=Shannon)) + geom_histogram()
plot2 = ggplot(data = joindata, aes(x=Observed)) + geom_histogram()
plot3 = ggplot(data = joindata, aes(x=Simpson)) + geom_histogram()
alphaplots = grid.arrange(grobs=list(plot1, plot2, plot3), nrow=1)
ggsave(alphaplots, file="3_GxE/3b_alpha_diversity_normality_check.png", width=10, height=4)

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

# Remove outliers - TODO - Confirm these are removing the correct values
for(myalpha in alpha_names){
  joindata[outliers[[myalpha]],myalpha] = NA
}

# Visualize the alpha diversity (& check for normality)
plot1 = ggplot(data = joindata, aes(x=Shannon)) + geom_histogram(fill="darkblue")
plot2 = ggplot(data = joindata, aes(x=Observed)) + geom_histogram(fill="darkblue")
plot3 = ggplot(data = joindata, aes(x=Simpson)) + geom_histogram(fill="darkblue")
alphaplots = grid.arrange(grobs=list(plot1, plot2, plot3), nrow=1)
ggsave(alphaplots, file="3_GxE/3b_alpha_diversity_normality_check.postfilter.png", width=10, height=4)


###########
# Calculate GxE for alpha diversity
###########

# Basic linear model
models = lapply(alpha_names, function(myalpha){
  myformula = paste(myalpha, "~", "location + Corrected_pedigree + location:Corrected_pedigree")
  lm(formula = myformula, data = joindata)
})

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
    strip.text.x = element_text(size=15,face="bold"),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    legend.position="bottom",
  )

#output the graph
ggsave(alphaplot, file="3_GxE/3b_alpha_diversity.variance_explained.png", height=5, width=6, device="png")

# Write data
write.csv(alpha_result, file="3_GxE/3b_alpha_diversity.variance_explained.csv", row.names=FALSE)

