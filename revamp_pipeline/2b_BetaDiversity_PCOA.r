# Load necessary libraries for beta diversity analysis
library(phyloseq)
library(ggplot2)
library(dplyr)
library(qiime2R) # Assuming this library provides the read_qza function

# Beta Diversity Analysis Using Weighted and Unweighted UniFrac Distances

# Read in the weighted and unweighted UniFrac distance matrices generated by QIIME2
weighted_unifrac <- read_qza("0_data_files/weighted_unifrac_pcoa_results.qza")
unweighted_unifrac <- read_qza("0_data_files/unweighted_unifrac_pcoa_results.qza")
G2F_Metadata_2019 <- read.table("0_data_files/G2f_2019_sub_corrected_metadata.tsv", sep = "\t", header = TRUE)

# Extract coordinate values from the PCoA results
weighted_unifrac_cordination <- weighted_unifrac$data$Vectors

# Focus on the first two principal coordinates (PC1 and PC2) which usually explain the most variation
weighted_unifrac_cordination <- weighted_unifrac_cordination[, c("SampleID", "PC1", "PC2")]

# Convert location metadata from factor to character for easier manipulation
G2F_Metadata_2019$location <- as.character(G2F_Metadata_2019$location)


# Enhance visibility by grouping locations into broader state groups
G2F_Metadata_2019_by_region <- G2F_Metadata_2019 %>%
  mutate(State_group = case_when(
    location %in% c("MOH1", "IAH2", "IAH4", "MNH1", "NEH1", "NEH2") ~ "Mid West",
    location %in% c("WIH1", "INH1", "MIH1", "OHH1") ~ "East Mississippi River",
    location %in% c("GAH1", "GAH2", "SCH1", "NCH1") ~ "South",
    location %in% c("NYH2", "NYH3", "DEH1") ~ "North East"
  ))

# Combine PCoA coordinates with metadata
input_pcoa_cordination_metadata <- inner_join(weighted_unifrac_cordination, G2F_Metadata_2019_by_region , by = "SampleID")

# Extract and format the percentage of variation explained by the first two PCs
pc1 <- as.character(round(weighted_unifrac$data$ProportionExplained[1, 1], 3) * 100)
pc2 <- as.character(round(weighted_unifrac$data$ProportionExplained[1, 2], 3) * 100)
pc1_label <- paste0("PC1 (", pc1, "%)")
pc2_label <- paste0("PC2 (", pc2, "%)")

# Plotting Weighted UniFrac Distances

# Generate a scatter plot of the weighted UniFrac distances, colored by state group
input_pcoa_plot <- ggplot(input_pcoa_cordination_metadata, aes(x = PC1, y = PC2, color = State_group)) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(pc1_label) +
  ylab(pc2_label) +
  scale_colour_manual(values = c("purple", "red", "steelblue", "darkgreen")) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold")
  )

# Display the plot
# input_pcoa_plot

# Save the plot to a file
ggsave("2_Diversity/2b_beta_diversity_weighted_unifrac.png", plot = input_pcoa_plot, height = 10, width = 10, device = "png")

# Additional analyses, such as dbRDA and heritability calculations, follow here...
