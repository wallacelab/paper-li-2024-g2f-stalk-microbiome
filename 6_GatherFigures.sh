#! /bin/bash

# Gather the figures for publication; not adding numbers since can change during revisions, etc.

# Set up output directory
figdir=6_Figures
if [[ ! -e "6_Figures" ]]; then mkdir $figdir; fi

# Overall 
cp 1_parsed_files/1b_g2f_locations.publication.png "$figdir/Figure - Location map.png"

# Core taxa
cp 4_CoreMicrobiome/4c_prevalence_summary.overall_pub.csv "$figdir/Table - Core Microbiome Summary.csv"
cp 4_CoreMicrobiome/4c_prevalence.genotypes.csv "$figdir/Table - Core Microbiome (by genotype).csv"
cp 4_CoreMicrobiome/4c_prevalence.locations.csv "$figdir/Table - Core Microbiome (by location).csv"

# Alpha diversity

# Beta diversity

# GxE
