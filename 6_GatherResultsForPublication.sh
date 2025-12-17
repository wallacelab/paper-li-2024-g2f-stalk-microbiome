#! /bin/bash

# Gather everything together for publication. (May be slight tweaks after this for some figures, but generally all here)

# Set up directory
pubdir=6_PublicationFiles
if [ ! -e $pubdir ]; then mkdir $pubdir; fi

##############
# Main Text
##############

# Figure 1 - Location map
cp 1_parsed_files/1b_g2f_locations.publication.png "$pubdir/figure.location_map.png"

# Figure 2 - Core heatmap
cp 4_CoreMicrobiome/4b_genotype_Class_core_microbiome_heatmap.all.clustered.png "$pubdir/figure.core_classes_genotype.png"
cp 4_CoreMicrobiome/4b_location_Class_core_microbiome_heatmap.all.clustered.png "$pubdir/figure.core_classes_location.png"

# Figure 3 - Alpha diversity
cp 2_Diversity/2a_alpha_diversity.plots.jgw.png "$pubdir/figure.alpha_diversity.png"

# Figure 4 - Beta diversity
cp 2_Diversity/2d_beta_diversity.weighted.publication.png "$pubdir/figure.beta_diversity.png"

# Figure 5 - GxE Breakdown
cp 3_GxE/3g_GxE_combined_plot.png "$pubdir/figure.gxe_breakdown.png"

# Figure 6 - Soil pH
cp 5_Associations/5a_median_alpha_diversity_plot.soil_ph.jgw.png "$pubdir/figure.soil_ph.png"
cp 5_Associations/5a_median_alpha_diversity_plot.buffer_ph.jgw.png "$pubdir/figure.buffer_ph.png"

# Figure 7 - Potassium
cp 5_Associations/5b_mantel_tests.pub.png "$pubdir/figure.potassium.png"

# Table 1 - Core taxa (this table gets heavily reformatted in the acutal paper)
cp 4_CoreMicrobiome/4d_short_summary.csv "$pubdir/table.core_taxa.csv"

#############
# Supplemental
#############

# Figure S1 - Core taxa % 
cp 4_CoreMicrobiome/4d_core_percent_by_sample.png "$pubdir/figure.core_taxa_percents.supplemental.png"

# Figure S2 - Weighted/unweighted unifrac
cp 2_Diversity/2d_beta_diversity.genotypes.supplemental.png "$pubdir/figure.beta_diversity.supplemental.png"

# Figure S3 - Metabolism only significant GxE
cp 3_GxE/3e_MetaCyc_pathway_GXE.significant_only.png "$pubdir/figure.metabolism_gxe_sig_only.supplemental.png"

# Supplemental tables 
Rscript 6a_MakeOverviewTables.r -o $pubdir
cp 4_CoreMicrobiome/4c_prevalence.genotypes.csv "$pubdir/table.core_by_genotype.supplemental.csv"
cp 4_CoreMicrobiome/4c_prevalence.locations.csv "$pubdir/table.core_by_location.supplemental.csv"
cp 2_Diversity/2a_alpha_diversity.kruskal_test_results.jgw.csv "$pubdir/table.alpha_diversity_significance.supplemental.csv"
cp 2_Diversity/2c_beta_diversity.permanova_results.csv "$pubdir/table.beta_diversity_significance.supplemental.csv"
cp 3_GxE/3b_alpha_diversity.variance_explained.jgw.csv "$pubdir/table.gxe_alpha.supplemental.csv"
cp 3_GxE/3c_beta_diversity_GxE.jgw.csv "$pubdir/table.gxe_beta.supplemental.csv"
cp 3_GxE/3e_MetaCyc_pathway_GXE.csv "$pubdir/table.gxe_pathways.supplemental.csv"
cp 3_GxE/3f_taxon_GXE.csv "$pubdir/table.gxe_taxonomy.supplemental.csv"
