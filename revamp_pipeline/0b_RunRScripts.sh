#! /bin/bash

# This is just a wrapper script to run all the R scripts without having to do them manually.

# Make needed directories
if [[ ! -e 1_parsed_files ]] ; then mkdir 1_parsed_files; fi
if [[ ! -e 2_Diversity ]] ; then mkdir 2_Diversity; fi
if [[ ! -e 3_GxE ]] ; then mkdir 3_GxE; fi
if [[ ! -e 4_CoreMicrobiome ]] ; then mkdir 4_CoreMicrobiome; fi
if [[ ! -e 5_Associations ]] ; then mkdir 5_Associations; fi
if [[ ! -e 5_Associations/5d_pathways ]] ; then mkdir 5_Associations/5d_pathways; fi

# Have script exit if encounters an error
set -e

# # Basic data unification and cleaning
# Rscript 1a_UnifyData.r
# Rscript 1b_PlotFieldLocations.r
# 
# # Alpha and beta diversity
# Rscript 2a_AlphaDiversity.r
# Rscript 2b_BetaDiversity_calc.r
# Rscript 2c_BetaDiversity_Significance.r
# Rscript 2d_BetaDiversity_PCOA.r
# 
# # GxE Breakdown
# Rscript 3a_FilterSamplesForGxE.r
# Rscript 3b_CalcGxE_AlphaDiversity.r
# Rscript 3c_CalcGxE_BetaDiversity.r
bash 3d_RunPicrust2OnGxEData.sh
Rscript 3e_CalcGxe_Pathways.r
Rscript 3f_CaclGxE_Taxonomy.r
Rscript 3g_PlotGxE.r

# Core microbiome
Rscript 4a_CoreMicrobiome.r
Rscript 4d_SummarizeCoreTaxa.r

# Environmental associations
Rscript 5a_AlphaDiversityAssociation.r
Rscript 5b_BetaDiversityAssociation.r
Rscript 5c_BetaDiversityAssociation_GenotypeLocation.r
Rscript 5d_PathwayAssociation.r

