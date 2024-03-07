#! /usr/bin/Rscript

# Load in and unify the data from QIIME, doing some minor filtering/correcting along the way
library(phyloseq)


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

##build phyloseq object 
ps	=	phyloseq(OTU,TAX,phy_tree,META)


#remove NCH1-1-159 samples for m
ps <- subset_samples(ps, rownames(sample_data(ps)) != "NCH1-1-159")

