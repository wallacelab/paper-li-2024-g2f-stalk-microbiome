# Genomes to Fields Maize Stalk Microbiome Analysis

This is the workflow for Li et al. 2024's analysis of the maize stalk microbiome across the Genomes to Fields initiative. This is specifically a revamped one to be more modular and reproducible than the original one. Currently it is being stored in a subfolder to make revision easier; eventually it will be moved to the main folder.

# Step 0 - ASV calling

This is the basic analysis step, going from raw sequencing data to files ready to be imported into R for analysis.

Importantly, these steps are not really meant to be run from a cloned repo, as they require the fastq files (which are WAY too big for Github). They are instead included so all the paramters used can be looked up (if necessary), and the key output files are included in 0_data_files

- **0a_qiime_workflow.sh**: This script uses QIIME2 to process raw data, build ASV table, and perform alpha and beta diversity analysis.
- **0_data_files**: This directory holds all the output files from the above. Subsequent scripts will load from here to access the original data.

# Post-ASV workflow

After calling ASVs, all subsequent steps are part of an R workflow. Unless you really want the raw analysis, this is probably the best place to pick up

- **0b_RunRScripts.sh**: Wrapper script to run all the R scripts in subsequent steps

## Step 1 - Data unification and filtering

This is a preparatory step to get the data combined and ready to use for subsequent analysis. A big part of this is creating RDS (R data structure) files with the ASV table pre-loaded as a Phyloseq object.

- **1a_UnifyData.r**: Unify the ASV table, metadata, taxonomy, & phylogenetic tree into a phyloseq object and write out. (Raw, rarefied, and rarefied + yellow stripe filtered)
    - Also makes a heatmap of genotype presence across locations
- **1b_PlotFieldLocations.r**: Make a plot of where the field locations are
- **1c_CalculateOtuTableStats.r**: Calculate basic stats on OTU tables so can easily get them for the manuscript
 
## Step 2 - Diversity Metrics

This step takes the ASV table from Step 1 and performs basic alpha and beta diversity analysis

- **2a_AlphaDiversity.r**: Calculate and plot basic alpha diversity measures
- **2b_BetaDiversity_calc.r**: Calculate beta diversity measures. Can take a while, so put in own script to be able to redo others (e.g., plotting) quickly
- **2c_BetaDiversity_Significance.r**: Calculate significance of beta diversity with PERMANOVA against location and genotype
- **2d_BetaDiversity_PCOA.r**: Plot PCoA plots of beta diversity

## Step 3 - GxE analysis

This step does variance analysis on various metrics to determine the relative importance of genetics, environment, and GxE to each of them.

- **3a_FilterSamplesForGxE.r**: Filter the ASV table down to just the samples to be used for GxE analysis
- **3b_CalcGxE_AlphaDiversity.r**: Perform GxE variance breakdown for alpha diversity
- **3c_CalcGxE_BetaDiversity.r**: Perform GxE variance breakdown for alpha diversity
- **3d_RunPicrust2OnGxEData.sh**: Run PICRUST2 (via a conda environment) to get pathway predictions
- **3e_CalcGxe_Pathways.r**: Perform GxE variance breakdown for PICRUST2 pathway predictions
- **3f_CaclGxE_Taxonomy.r**: Perform GxE variance breakdown for bacterial taxonomic levels
- **3g_PlotGxE.r**: Combine the individual outputs of the above scripts into a single figure for publication

## Step 4 - Core Microbiome

This step determines the core microbiome across samples

- **4a_CoreMicrobiome.r**: Calculate the "core" microbiome
- **4d_SummarizeCoreTaxa.r**: Make summary tables of the core microbiome

## Step 5 - Association

This step associates various metrics with environmental variables to try to determine what could be driving some of the community makeup.

- **5a_AlphaDiversityAssociation.r**: Test environmental variables against alpha diversity
- **5b_BetaDiversityAssociation.r**: Test environmental variables against beta diversity
- **5c_BetaDiversityAssociation_GenotypeLocation.r**: More testing of beta diversity, this time of specific genotypes across locations
- **5d_PathwayAssociation.r**: Test environmental variables against PICRUST2 predicted pathways
