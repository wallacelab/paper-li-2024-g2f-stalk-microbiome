#! /usr/bin/Rscript

# Graph the core taxa for locations and genotypes using Metacoder

library(tidyverse)
library(metacoder)

# Load prevalence tables
genotable = read.csv("4_CoreMicrobiome/4c_prevalence_summary.genotypes.csv")
loctable= read.csv("4_CoreMicrobiome/4c_prevalence_summary.locations.csv")

# Function to make actual tree
make_heat_tree = function(mydata, hide_threshold=0, outfile="99_tmp.png"){
  # Make taxa table
  mytaxa = parse_tax_data(tax_data=mydata,
                            class_cols="taxon_name",
                            class_sep=";",
                            class_regex = "^(.+)__(.+)$",
                            class_key = c(tax_rank = "info", # A key describing each regex capture group
                                        tax_name = "taxon_name"),
                            named_by_rank=FALSE)
  
  # Find missing nodes & add so can plot properly
  allnodes = names(mytaxa$taxa)
  missing = setdiff(allnodes, mytaxa$data$tax_data$taxon_id)
  if(length(missing)>0){
    fillin = lapply(missing, function(m){
      myname = mytaxa$taxa[[m]]$get_name()
      data.frame(taxon_id=m, taxon_level="ROOT?", fraction_core=1,
                 taxon_name=myname)
    }) %>%
      bind_rows()
    mytaxa$data$tax_data = 
      bind_rows(fillin, mytaxa$data$tax_data)
  }
  
  # Hide names of taxa at too low of prevalence
  treedata = mytaxa$data$tax_data
  defaults = sapply(treedata$taxon_id, function(i){
    mytaxa$taxa[[i]]$name$name
  })
  defaultkey = data.frame(taxon_id = names(defaults), 
                          plotname = as.character(defaults))
  treedata = left_join(treedata, defaultkey, by="taxon_id")
  labels = ifelse(treedata$fraction_core >= hide_threshold,
                  yes = treedata$plotname,
                  no = "")
  # if(nrow(tohide) > 0){
  #   for(node in tohide$taxon_id){
  #     mytaxa$taxa[[node]]$name = "." # Doesn't work. Need to be unique?
  #   }
  # }
  
  # Plot
  set.seed(1)
  treeplot = heat_tree(
    mytaxa,
    node_label = labels,
    #node_label   = taxon_names,           # show taxon names
    node_color   = fraction_core,   # map prevalence to color
    node_size    = n_subtaxa,             # size by number of descendants (change if desired)
    #  node_color_axis_label = "Fraction where is core",
    #node_size_axis_label  = "Number of subtaxa",
    layout = "davidson-harel", initial_layout = "reingold-tilford",
    title = "Taxonomic prevalence heat tree"
  )
  
  ggsave(treeplot, file=outfile)
}

make_heat_tree(genotable, hide_threshold=0.6,
                          outfile="4_CoreMicrobiome/4e_prevalence_heat_tree.genotypes.png")
make_heat_tree(loctable, hide_threshold=0.6,
                        outfile="4_CoreMicrobiome/4e_prevalence_heat_tree.locations.png")

