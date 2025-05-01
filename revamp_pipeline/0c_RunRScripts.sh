#! /bin/bash

# This is just a wrapper script to run all the R scripts without having to do them manually.


# Basic data unification and cleaning
Rscript 1a_UnifyData.r
Rscript 1b_PlotFieldLocations.r

# Alpha and beta diversity
Rscript 2a_AlphaDiversity.r
Rscript 2b_BetaDiversity_PCOA.r

# GxE Breakdown
Rscript 3a_FilterSamplesForGxE.r
Rscript 3b_CalcGxE_AlphaDiversity.r
Rscript 3c_CalcGxE_BetaDiversity.r
Rscript 3d_CalcGxe_Pathways.r
Rscript 3e1_CompareGxE_Taxonomy.r
Rscript 3e_CaclGxE_Taxonomy.r
Rscript 3f_PlotGxE.r

# Core microbiome
Rscript 4a_CoreMicrobiome.r
Rscript 4d_SummarizeCoreTaxa.r

# Environmental associations
Rscript 5a_AlphaDiversityAssociation.r
Rscript 5b_BetaDiversityAssociation.r
Rscript 5c_BetaDiversityAssociation_GenotypeLocation.r
Rscript 5d_PathwayAssociation.r

