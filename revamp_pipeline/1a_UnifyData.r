#! /usr/bin/Rscript

# Load in and unify the data from QIIME, doing some minor filtering/correcting along the way
library(phyloseq)
library(tidyverse)
library(gridExtra)

rarefy_depth=1500  # Depth to rarify samples to
rarefy_seed=1 # Random number seed for rarefaction
yellow_stripes = c("PHW52/PHM49","B73/PHM49","F42/H95","F42/MO17","OH43/B37","PHW52/PHN82","B14A/OH43","B14A/H95","B14A/MO17",
                   "LH74/PHN82","PHG39/PHN82","B73/MO17","B73/PHN82","B37/H95","F42/OH43","B37/MO17","B37/OH43","CG119/CG108",
                   "CG44/CGR01","2369/LH123HT")  # "Yellow-stripe" genotypes that are the focus of this study

#############
# Load data
#############


#Read in the ASV table after mitochondira and chloroplast filteration
otu_table	=	read.table("0_data_files/dada2_table-no-mitochondria-no-chloroplast.txt.gz",
                       row.names =1, check.names = FALSE, sep="\t", skip=1, header = TRUE)

#taxonomy come from qiime2 artifact and needs to seperate in domain, phylum,order,......
taxonomy = read.csv("0_data_files/taxonomy.tsv",sep="\t",row.names =1)

#import phylogenetic tree table
phy_tree	=	read_tree("0_data_files/tree.nwk")

#import metadata table and convert to matrix 
metadata	=	read.table("0_data_files/G2f_2019_sub_corrected_metadata.tsv",sep ="\t",header = 1,row.names=1)


#############
# Remove taxa in blank samples
#############

# Identify the blank taxa 
is_blank = grepl(names(otu_table), pattern="^B[0-9]+$", perl=TRUE) # Blanks are B1, B2, B3, etc.
blank_taxa <- otu_table[rowSums(otu_table[,is_blank]) > 1,] # Taxa that show up in at least 1 blank

#filter out the blank taxa 
otu_table_blank_taxa_filtered <- otu_table[rowSums(otu_table[,is_blank]) < 1,] # Keep only those taxa that have 0 reads in blank samples

#convert to matrix 
otu_table	=	as.matrix(otu_table_blank_taxa_filtered)

#############
# Clean up taxonomy table
#############

#revise the taxonomy so that uncultured family have previous taxonomy information
revised_taxonomy <- taxonomy

#change the column to character 
revised_taxonomy$Family <- as.character(revised_taxonomy$Family)
revised_taxonomy$Order <- as.character(revised_taxonomy$Order)
revised_taxonomy$Class <- as.character(revised_taxonomy$Class)

#select family taxa that have no taxonomy information until phylumn level and paste the phylumn level information with family information
uncultured_class = revised_taxonomy$Family ==" f__uncultured"& revised_taxonomy$Order ==" o__uncultured" & revised_taxonomy$Class ==" c__uncultured"
revised_taxonomy$Family[uncultured_class] <-  paste(revised_taxonomy$Family[uncultured_class],revised_taxonomy$Phylum[uncultured_class])

#select family taxa that have no taxonomy information until class level and paste the class level information with family information
uncultured_order = revised_taxonomy$Family ==" f__uncultured"& revised_taxonomy$Order ==" o__uncultured"
revised_taxonomy$Family[uncultured_order] <-  paste(revised_taxonomy$Family[uncultured_order],revised_taxonomy$Class[uncultured_order])

#select family taxa that have no taxonomy information until order level and paste the order level information with family information
uncultured_family = revised_taxonomy$Family ==" f__uncultured"
revised_taxonomy$Family[uncultured_family] <-  paste(revised_taxonomy$Family[uncultured_family],revised_taxonomy$Order[uncultured_family])

# revise factor to character to avoid problem in column operation  
revised_taxonomy$Family <- as.factor(revised_taxonomy$Family)
revised_taxonomy$Order <- as.factor(revised_taxonomy$Order)
revised_taxonomy$Class <- as.factor(revised_taxonomy$Class)
taxonomy <- as.matrix(revised_taxonomy)

###############
# Convert to Phylsoeq object
###############

#format the input data 
OTU	=	otu_table(otu_table,taxa_are_rows	=	TRUE)
TAX	=	tax_table(taxonomy)
META = sample_data(metadata)

#build phyloseq object 
ps	=	phyloseq(OTU,TAX,phy_tree,META)

# Rarefy to specific depth using the seed value from earlier
ps.rarefied = rarefy_even_depth(ps, sample.size=rarefy_depth, replace=F, rngseed=rarefy_seed)

# Limit to just the "yellow-stripe" maize varieties
ps.rarefied_ys_filtered <- subset_samples(ps.rarefied, Corrected_pedigree %in% yellow_stripes)


###############
# Write out unified phyloseq objects
###############

saveRDS(ps, file="1_parsed_files/1a_asv_table_no_taxa_from_blanks.phyloseq.rds")
saveRDS(ps.rarefied, file=paste("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy",rarefy_depth, ".phyloseq.rds", sep=""))
saveRDS(ps.rarefied_ys_filtered, file=paste("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy",rarefy_depth, ".yellow_stripe.phyloseq.rds", sep=""))

###############
# Graph presence of genotypes in each location
###############

plot_heatmap = function(mydata){
  toplot = sample_data(mydata) %>%
    group_by(Corrected_pedigree, location) %>%
    summarize(count=n(), .groups="drop")
  myplot = ggplot(toplot) + 
    aes(x=location,y=Corrected_pedigree, fill=count) +
    geom_tile() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))
  return(myplot)
}

# Make plots
rawplot = plot_heatmap(ps) + ggtitle("Raw data")
rareplot = plot_heatmap(ps.rarefied) + ggtitle("Rarefied")
ysplot = plot_heatmap(ps.rarefied_ys_filtered) + ggtitle("Rarefied + Yellow Stripe")

output = grid.arrange(grobs=list(rawplot, rareplot, ysplot), nrow=1)
ggsave(output, file="1_parsed_files/1a_genotype_by_location_heatmap.png", width=16, height=5)


