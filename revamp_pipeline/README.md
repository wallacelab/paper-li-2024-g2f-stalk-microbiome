# Genomes to Fields Maize Stalk Microbiome Analysis

This is the workflow for Li et al. 2024's analysis of the maize stalk microbiome across the Genomes to Fields initiative. This is specifically a revamped one to be more modular and reproducible than the original one. Currently it is being stored in a subfolder to make revision easier; eventually it will be moved to the main folder.

# Step 0 - ASV calling and Picrust analysis

This is the basic analysis step, going from raw sequencing data to files ready to be imported into R for analysis.

Importantly, these steps are not really meant to be run from a cloned repo, as they require the fastq files (which are WAY too big for Github). They are instead included so all the paramters used can be looked up (if necessary), and the key output files are included in 0_data_files

**0a_qiime_workflow.sh**

This script uses QIIME2 to process raw data, build ASV table, and perform alpha and beta diversity analysis.

**0b_picrust2_workflow.sh**

This bash script was used to run picrust2 on UGA's high-performance computing cluster, sapelo2.

**0_data_files**

This directory holds all the output files from the above. Subsequent scripts will load from here to access the original data.

# Step 1 - Data unification and filtering



# Step 2 - TODO


# Step 3 - TODO




**Explanation for files**


**2019_G2F_metadata_enviromental.tsv**

Metadata of 2019 G2F samples


**G2F_metadata_2019_duplicate_pedigree.txt**

SampleID of G2F duplicate samples 



**aligned-rep-seqs.qza**

representative sequence 


**alpha_diversity_metrics.txt**

alpha diversity name for loop operation 


**beta_diversity_metrics.txt**

beta diversity name for loop operation 


**column_name_list.txt**

Metadata column name 


**dada2_table-no-mitochondria-no-chloroplast-blank-filtered-yellow-stripe-duplicate-pedigree.qza**

ASV table after filtering mitochondria, chloroplast, and blank, keeping samples with duplicate 


**dada2_table-no-mitochondria-no-chloroplast-blank-filtered.qza**

ASV table after filtering mitochondria, chloroplast, and blank


**dada2_table-no-mitochondria-no-chloroplast.qza**

ASV table after filtering mitochondria and chloroplast 


**dada2_table.qza**

original ASV table 


**path_abun_unstrat_descrip.tsv**

Pathway estimation based on ASV table with only duplicated samples of YS pedigree 






**yellow_strip_gbs_data_orinigal.tsv**

The file that contains the original yellow stripe species name 
