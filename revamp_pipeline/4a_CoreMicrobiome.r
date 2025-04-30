#! /usr/bin/Rscript

# Identify the core microbiome of the samples

library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(phyloseq)
library(qiime2R)
library(parallel)
library(ggnewscale)


# Global variables
min_fraction = 0.001 # Min fraction of a sample an OTU has to consist of to "count" for core calculations
min_prevalence = 0.6 # Fraction of samples a taxon has to be in to be "core"
ranks=c("Phylum", "Class", "Order", "Family", "Genus", "Species") # Taxonomic ranks to test at
ncores=7 # Number of CPU cores to use for mclapply()

# Load data 

#asvs = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")
asvs = readRDS("3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds") # To be consistent with GxE analysis
metadata = asvs %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")

# Precompute collapsed taxon profiles (takes a while)
collapsed_taxa = mclapply(ranks, function(myrank){
  tax_glom(asvs, taxrank = myrank)
}, mc.cores=ncores)
names(collapsed_taxa) = ranks

###############
# Helper functions
###############

# Function to construct a name key
get_namekey = function(mydata){
  namekey = tax_table(mydata) %>%
    as.data.frame() %>%
    apply(MARGIN=1, FUN=paste, collapse=";")
  namekey = data.frame(Taxaname = namekey) %>%
    rownames_to_column("taxon")
  return(namekey)
}

# Preparse the ASV data down to just locations and genotypes
core_microbiome_preparse <- function(tax_level, mydata, collapsed_taxa){

  # Grab collapsed data & name key at given taxonomic level
  collapsed = collapsed_taxa[[tax_level]]
  namekey = get_namekey(collapsed)
  mytaxa = rownames(tax_table(collapsed))
  
  # Calculate prevalence at each location
  # Start by reformatting data into tidy data with location, taxon, and count
  collapsed_counts = otu_table(collapsed) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("SampleID") %>%
    left_join(metadata, by="SampleID") %>%
    #relocate(SampleID, Corrected_pedigree, location) %>%
    select(SampleID, Corrected_pedigree, location, all_of(mytaxa)) %>%
    pivot_longer(cols=all_of(mytaxa), names_to="taxon", values_to="count") %>%
    group_by(SampleID) %>%
    mutate(fraction=count / sum(count), present = fraction >= min_fraction)%>% # True/false for if taxon present in a sample at at least min_fraction abundance
    #mutate(present = count > 0) %>% # Add a true/false column for if taxon present
    left_join(namekey, by="taxon")
  
  return(collapsed_counts)
}

# Plot heatmap
plot_heatmap_custom = function(mycore, myrank, xval, outfile, plot_all=FALSE){
  # Format taxonomy
  taxonomy = str_split_fixed(mycore$Taxaname, ';', n=Inf)
  colnames(taxonomy) = rank_names(asvs)
  mycore$plot_taxon = taxonomy[,myrank]
  
  # Set which variable is on X axis
  mycore$xval = mycore[,xval] %>% unlist()

  # Set significance marker
  mycore$significant = ifelse(mycore$prevalence >= min_prevalence,
                              yes="*", no="")
    
  # Set mapping aesthetics
  my_aes = aes(x=xval, y=plot_taxon, label=significant)
  if(plot_all){ # Add prevalence as a fill if specified
    my_aes = aes(x=xval, y=plot_taxon, label=significant, fill=prevalence)
  }
  
  # Make plot
  #colorscale = c("#000000","#2b2b2b", "#000080","#0000FF")
  #colorbreaks = c(0, min_prevalence*0.999, min_prevalence, 1)
  colorscale = c("#FFFFFF","#0000dd")
  colorbreaks = c(0, 1)
  myplot = ggplot(mycore) +
    my_aes +
    geom_tile() + 
    geom_text(color="white") +
    theme(
      axis.text = element_text(size=15, face="bold"),
      axis.text.x = element_text(angle=90),
      axis.title = element_text(size=15,face="bold")
    ) +
    labs(y = myrank, x=xval, 
         title=paste("Core microbiome of samples (only counts when >=",min_fraction,"of a sample)") )+
  scale_fill_gradientn(colors=colorscale, values=colorbreaks) +
    scale_y_discrete(limits=rev)   # Reverse Y axis
  
  # Save plot
  ggsave(myplot, file=outfile, height=8, width=12)
}

# Make a quick lookup key from 2 vectors (e.g., location and taxon) (putting here so is consistent)
make_lookup = function(x, y){
  paste(x,y, sep="_")
}

##########
# Preparse data
##########

# Preparse the ASV data to make it quicker to plot
preparse = lapply(ranks, FUN=core_microbiome_preparse, mydata = asvs, collapsed_taxa=collapsed_taxa)
names(preparse) = ranks

###############
# Prevalence by location
###############

# General prevalence
location_prevalence = lapply(ranks, function(myrank){
  prevalence = preparse[[myrank]] %>%
    #filter(fraction >= min_fraction) %>%
    group_by(location, taxon, Taxaname) %>%
    summarize(prevalence = sum(present) / n(), .groups="drop")
  return(prevalence)
})
names(location_prevalence) = ranks

# Filter for core microbiome
location_cores = lapply(location_prevalence, function(myprev){
  myprev %>%
    filter(prevalence >= min_prevalence)
})

# Core but with actual prevalence values retained
location_cores_all = lapply(ranks, function(myrank){
  location_prevalence[[myrank]] %>%
    filter(taxon %in% location_cores[[myrank]]$taxon)
})
names(location_cores_all) = ranks

# Plot and save
locplots = lapply(ranks, function(myrank){
  # Basic binary plot
  mycore = location_cores[[myrank]]
  outfile=paste("4_CoreMicrobiome/4a_location_",myrank,"_core_microbiome_heatmap.jgw.png", sep="")
  plot_heatmap_custom(mycore, myrank, xval="location", outfile=outfile)

  # More thorough heatmap
  mycore = location_cores_all[[myrank]]
  outfile=paste("4_CoreMicrobiome/4b_location_",myrank,"_core_microbiome_heatmap.all.jgw.png", sep="")
  plot_heatmap_custom(mycore, myrank, xval="location", outfile=outfile, plot_all=TRUE)
  
})

# Make a quick lookup key to check for what is core (for later)
location_core_key = lapply(ranks, function(r){
  mycore = location_cores[[r]]
  make_lookup(mycore$location, mycore$taxon)
  mykey = data.frame(location=mycore$location, taxon=mycore$taxon,
                     key = make_lookup(mycore$location, mycore$taxon))
})
names(location_core_key) = ranks


##################
# Core by Genotype
##################

genotype_prevalence = lapply(ranks, function(myrank){
  prevalence = preparse[[myrank]] %>%
    group_by(Corrected_pedigree, taxon, Taxaname) %>%
    summarize(prevalence = sum(present) / n(), .groups="drop")
  return(prevalence)
})
names(genotype_prevalence) = ranks

# Filter for core
genotype_cores = lapply(genotype_prevalence, function(myprev){
  myprev %>%
    filter(prevalence >= min_prevalence)
})

# Core but with actual prevalence values retained
genotype_cores_all = lapply(ranks, function(myrank){
  genotype_prevalence[[myrank]] %>%
    filter(taxon %in% genotype_cores[[myrank]]$taxon)
})
names(genotype_cores_all) = ranks



# Plot and save
genoplots = lapply(ranks, function(myrank){
  mycore = genotype_cores[[myrank]]
  outfile=paste("4_CoreMicrobiome/4a_genotype_",myrank,"_core_microbiome_heatmap.jgw.png", sep="")
  plot_heatmap_custom(mycore, myrank, xval="Corrected_pedigree", outfile=outfile)
  
  mycore = genotype_cores_all[[myrank]]
  outfile=paste("4_CoreMicrobiome/4b_genotype_",myrank,"_core_microbiome_heatmap.all.jgw.png", sep="")
  plot_heatmap_custom(mycore, myrank, xval="Corrected_pedigree", outfile=outfile, plot_all=TRUE)
  
})

# Make a quick lookup key to check for what is core (for later)
geno_core_key = lapply(ranks, function(r){
  mycore = genotype_cores[[r]]
  mykey = data.frame(geno=mycore$Corrected_pedigree, taxon=mycore$taxon,
                     key = make_lookup(mycore$Corrected_pedigree, mycore$taxon))
})
names(geno_core_key) = ranks

####################
# Overall heatmap of all samples
####################

make_overall_heatmap = function(collapsed, rank="Phylum", facet_by="none", 
                                min_fraction, min_prevalence, mybreaks){
  # Get data and make relative abundance
  mydata=collapsed[[rank]] 
  myfraction = mydata %>%
    transform_sample_counts(fun=function(x) x / sum(x) )
  
  # Determine which OTUs count as core
  mytable=otu_table(myfraction)
  fraction_present = rowSums(mytable >= min_fraction, na.rm=TRUE) / ncol(mytable)
  is_core = fraction_present >= min_prevalence
  core_taxa = names(fraction_present)[is_core]

  # Subset to just the core taxa
  if(length(core_taxa)==0){
    cat("No remaining taxa at rank",rank,"\n")
    return(NULL)
  }
  mycore = prune_taxa(taxa_names(myfraction) %in% core_taxa, myfraction)
  # taxa_names(mycore) = tax_table(mycore)[,rank]
  
  # Make taxa names
  taxlevels = mycore %>% tax_table() %>% colnames()
  thisrank = which(taxlevels==rank)
  phyrank = which(taxlevels=="Phylum")
  mylevels = tax_table(mycore)[,phyrank:thisrank]
  ## Paste together
  mynames = apply(mylevels, MARGIN=1, FUN=paste, collapse=" - ") %>%
    gsub(pattern=".__", repl="") %>%
    gsub(pattern=" ", repl="")
  taxa_names(mycore) = mynames

  # Reformat for plotting
  samplekey = sample_data(mycore) %>%
    data.frame() %>%
    rownames_to_column("SampleID")
  plotdata = otu_table(mycore) %>%
    as.data.frame() %>%
    rownames_to_column("taxon") %>%
    pivot_longer(cols=-taxon, names_to="SampleID", values_to="abundance") %>%
    left_join(samplekey, by="SampleID") %>%
    mutate(abundance_adj = abundance + min_fraction,
           abundance_log = log(abundance_adj))
  

  # Make heatmap
  overall_map = ggplot(plotdata) +
    #aes(x=SampleID, y=taxon, fill=abundance_log) +
    aes(x=SampleID, y=taxon, fill=abundance_adj) +
    geom_tile() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.spacing=unit(0.2, "mm"),
          axis.title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold", size=12)) +
    scale_y_discrete(limits=rev)  + # Reverse Y axis
    scale_fill_gradient(trans = "log", breaks=mybreaks, limits=range(mybreaks)) +
    labs(fill=str_wrap("Relative Abundance", width=10),
         y=paste("Taxonomic",rank))

  # Facet plot by location or genotype if requested
  if(facet_by=="none"){ 
    overall_map = overall_map +
      labs(x="Samples")
    ggsave(overall_map, file=paste("4_CoreMicrobiome/4b_overall_",rank,"_core_microbiome_heatmap.jgw.png", sep=""), width=8, height=5)
  } else if (facet_by=="location"){
    overall_map = overall_map +
      facet_wrap(~location, scales="free_x", nrow=1, strip.position = "bottom") +
      labs(x="Samples by location")
    ggsave(overall_map, file=paste("4_CoreMicrobiome/4b_overall_",rank,"_core_microbiome_heatmap.by_location.jgw.png", sep=""), width=13, height=5)
  } else if(facet_by=="genotype"){
    overall_map = overall_map +
      facet_wrap(~Corrected_pedigree, scales="free_x", nrow=1, strip.position = "bottom") +
      theme(strip.text = element_text(size=7)) +
      labs(x="Samples by Genotype")
    ggsave(overall_map, file=paste("4_CoreMicrobiome/4b_overall_",rank,"_core_microbiome_heatmap.by_genotype.jgw.png", sep=""), width=18, height=5)
  }
  #return(overall_map)

}

mybreaks = c(0.001, 0.01, 0.1, 1)
lapply(ranks, make_overall_heatmap, facet_by="none", collapsed=collapsed_taxa, min_fraction=min_fraction, min_prevalence=min_prevalence, mybreaks=mybreaks)
lapply(ranks, make_overall_heatmap, facet_by="location", collapsed=collapsed_taxa, min_fraction=min_fraction, min_prevalence=min_prevalence, mybreaks=mybreaks)
lapply(ranks, make_overall_heatmap, facet_by="genotype", collapsed=collapsed_taxa, min_fraction=min_fraction, min_prevalence=min_prevalence, mybreaks=mybreaks)




####################
# Overall heatmap with 2-tone coloring
####################

make_overall_heatmap = function(collapsed, rank="Phylum", facet_by="none", lookup_key, mybreaks){
  # Get data and make relative abundance
  mydata=collapsed[[rank]] 
  myfraction = mydata %>%
    transform_sample_counts(fun=function(x) x / sum(x) )
  
  # Determine which OTUs count as core
  mylookup = lookup_key[[rank]]
  core_taxa = unique(mylookup$taxon)
  
  # Subset to just the core taxa
  if(length(core_taxa)==0){
    cat("No remaining taxa at rank",rank,"\n")
    return(NULL)
  }
  mycore = prune_taxa(taxa_names(myfraction) %in% core_taxa, myfraction)
  # taxa_names(mycore) = tax_table(mycore)[,rank]
  
  # Make taxa names
  taxlevels = mycore %>% tax_table() %>% colnames()
  thisrank = which(taxlevels==rank)
  phyrank = which(taxlevels=="Phylum")
  mylevels = tax_table(mycore)[,phyrank:thisrank]
  ## Paste together
  mynames = apply(mylevels, MARGIN=1, FUN=paste, collapse=" - ") %>%
    gsub(pattern=".__", repl="") %>%
    gsub(pattern=" ", repl="") %>%
    sub(pattern="-.+-(.+)", repl=" (\\1)")
  taxa_names(mycore) = mynames
  ## Save original taxon names to pair with loojup for 2-tone plotting
  orig_names = names(mynames)
  names(orig_names) = mynames
  
  # Reformat for plotting
  samplekey = sample_data(mycore) %>%
    data.frame() %>%
    rownames_to_column("SampleID")
  plotdata = otu_table(mycore) %>%
    as.data.frame() %>%
    rownames_to_column("taxon") %>%
    pivot_longer(cols=-taxon, names_to="SampleID", values_to="abundance") %>%
    left_join(samplekey, by="SampleID") %>%
    mutate(abundance_adj = abundance + min_fraction,
           abundance_log = log(abundance_adj),
           orig_name = orig_names[taxon])
  
  # Take lookup key into account (for two-tone coloring)
  if(facet_by == "genotype" && !is.null(lookup_key)){
    plotdata$key = make_lookup(plotdata$Corrected_pedigree, plotdata$orig_name)
    plotdata$is_core = plotdata$key %in% mylookup$key
  } else if(facet_by == "location" && !is.null(lookup_key)){
    plotdata$key = make_lookup(plotdata$location, plotdata$orig_name)
    plotdata$is_core = plotdata$key %in% mylookup$key
  } else{
    plotdata$is_core = TRUE  # Set to all true by default
  }
  
  # Make heatmap
  overall_map = ggplot(plotdata) +
    #aes(x=SampleID, y=taxon, fill=abundance_log) +
    aes(x=SampleID, y=taxon, fill=abundance_adj) +
    geom_tile() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.spacing=unit(0.2, "mm"),
          axis.title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold", size=12)) +
    scale_y_discrete(limits=rev)  + # Reverse Y axis
    scale_fill_gradient(trans = "log", breaks=mybreaks, limits=range(mybreaks),
                        low="#000000", high="#56B1F7",
                        guide=guide_legend(order=1, reverse=TRUE,
                                           title=str_wrap("Relative Abundance (core)", width=10))) +
    labs(#fill=str_wrap("Relative Abundance (core)", width=10),
         y=paste("Taxonomic",rank))
  
  # Deal with coloring for core/noncore (but only if there are noncore things to plot)
  if(any(!plotdata$is_core)){
    overall_map = overall_map +
      ggnewscale::new_scale_fill() +
      aes(x=SampleID, y=taxon, fill=abundance_adj) +
      geom_tile(mapping=aes(alpha=as.integer(!is_core))) +
      scale_fill_gradient(trans = "log", breaks=mybreaks, limits=range(mybreaks),
                          low="#000000", high="#f75656",
                          guide=guide_legend(order=2, reverse=TRUE,
                                             title=str_wrap("Relative Abundance (noncore)", width=10))) +
      scale_alpha_continuous(guide="none")
    
  }
  
  # Facet plot by location or genotype if requested
  if(facet_by=="none"){ 
    overall_map = overall_map +
      labs(x="Samples by location")
    ggsave(overall_map, file=paste("4_CoreMicrobiome/4b_twotone_",rank,"_core_microbiome_heatmap.jgw.png", sep=""), width=8, height=5)
  } else if (facet_by=="location"){
    overall_map = overall_map +
      facet_wrap(~location, scales="free_x", nrow=1, strip.position = "bottom") +
      labs(x="Samples by location") +
      theme(strip.text = element_text(size=7))
    ggsave(overall_map, file=paste("4_CoreMicrobiome/4b_twotone_",rank,"_core_microbiome_heatmap.by_location.jgw.png", sep=""), width=8, height=5)
  } else if(facet_by=="genotype"){
    overall_map = overall_map +
      facet_wrap(~Corrected_pedigree, scales="free_x", nrow=1, strip.position = "bottom") +
      theme(strip.text = element_text(size=7, angle=90)) +
      labs(x="Samples by Genotype")
    ggsave(overall_map, file=paste("4_CoreMicrobiome/4b_twotone_",rank,"_core_microbiome_heatmap.by_genotype.jgw.png", sep=""), width=8, height=5)
  }
  #return(overall_map)
  
}

mybreaks = c(0.001, 0.01, 0.1, 1)
lapply(ranks, make_overall_heatmap, facet_by="location", lookup_key = location_core_key, collapsed=collapsed_taxa, mybreaks=mybreaks)
lapply(ranks, make_overall_heatmap, facet_by="genotype", lookup_key = geno_core_key, collapsed=collapsed_taxa, mybreaks=mybreaks)


##############
# Output tables of results
##############

#The location_cores_all and genotype_cores_all lists have prevalence for everything

# Helper to turn the *_core_all into an output table
make_output_table = function(mydata){
  lapply(ranks, function(myrank){
    mydata[[myrank]] %>%
      mutate(taxon_level=myrank, is_core = prevalence >= min_prevalence) %>%
      relocate(taxon_level) 
  }) %>%
    bind_rows() %>%
    rename(taxon_name = Taxaname) %>%
    mutate(taxon_name = sub(taxon_name, pattern="(;NA)+", repl=""),
           prevalence = round(prevalence, digits=3),
           taxon_level = factor(taxon_level, levels=ranks))
}

summarize_cores = function(mydata){
  mydata %>%
    group_by(taxon_level, taxon_name) %>%
    summarize(fraction_core = sum(prevalence >= min_prevalence) / n()) %>%
    relocate(taxon_level, fraction_core) %>%
    mutate(taxon_level = factor(taxon_level, levels=ranks),
           fraction_core = round(fraction_core, digits=3)) %>%
    arrange(as.numeric(taxon_level), desc(fraction_core), taxon_name)
}

# Make location table
location_cores_output = make_output_table(location_cores_all) %>%
  select(taxon_level, location, prevalence, is_core, taxon_name) %>%
  arrange(as.numeric(taxon_level), taxon_name, location)
write.csv(location_cores_output, file="4_CoreMicrobiome/4c_prevalence.locations.csv", row.names=FALSE)

# Calculate how often is core across locations
location_summary = summarize_cores(location_cores_output)
write.csv(location_summary, file="4_CoreMicrobiome/4c_prevalence_summary.locations.csv", row.names=FALSE)

# Make genotype table
genotype_cores_output = make_output_table(genotype_cores_all) %>%
  rename(genotype="Corrected_pedigree") %>%
  select(taxon_level, genotype, prevalence, is_core, taxon_name) %>%
  arrange(as.numeric(taxon_level), taxon_name, genotype)
write.csv(location_cores_output, file="4_CoreMicrobiome/4c_prevalence.genotypes.csv", row.names=FALSE)

# Calculate how often is core across locations
location_summary = summarize_cores(location_cores_output)
write.csv(location_summary, file="4_CoreMicrobiome/4c_prevalence_summary.genotypes.csv", row.names=FALSE)


# Now need to make the overall table
overall_prevalence = lapply(ranks, function(myrank){
  prevalence = preparse[[myrank]] %>%
    #filter(fraction >= min_fraction) %>%
    group_by(taxon, Taxaname) %>%
    summarize(prevalence = sum(present) / n(), .groups="drop")
  return(prevalence)
})
names(overall_prevalence) = ranks

# And write out (this is already basically at the summary level, since is by individual samples)
overall_cores_output = make_output_table(overall_prevalence) %>%
  select(taxon_level, prevalence, is_core, taxon_name) %>%
  arrange(as.numeric(taxon_level), desc(prevalence), taxon_name)
write.csv(overall_cores_output, file="4_CoreMicrobiome/4c_prevalence_summary.overall.csv", row.names=FALSE)


