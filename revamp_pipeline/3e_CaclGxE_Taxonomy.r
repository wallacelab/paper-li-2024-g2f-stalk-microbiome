#! /usr/bin/Rscript

# Original code copied from Roy's file /media/jgwall/Expansion/hl46161/R/G2F_analysis/second_scratch.R

library(qiime2R)
library(phyloseq)
library(tidyverse)

# Load data
metadata = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.phyloseq.rds") %>% 
  sample_data() %>%
  data.frame() %>%
  rownames_to_column("SampleID")
taxa_levels=c("phylum", "class","order","family","genus","species")
core_taxa = lapply(taxa_levels, function(t){
  infile=paste("0_data_files/core_taxa/dada2_table-no-mitochondria-no-chloroplast-blank-filtered-YS-only-",t,"_25_filtered.qza", sep="")
  read_qza(infile)$data
})
names(core_taxa) = taxa_levels

############
# ANOVA calculations
############

# Helper function to calculate ANOVA results for taxa
anova_calculation_taxa <- function(input_data,input_metadata){
  
  # Reformat core taxa table
  core_taxa_table <- as.data.frame(t(input_data))
  
  # Add pseudocoount of 1
  pseudo_core_taxa_table <- core_taxa_table + 1
  
  # Create sample ID column 
  pseudo_core_taxa_table$SampleID <- rownames(pseudo_core_taxa_table)
  
  # Join with metadata
  input_table_core_taxa_meatadata_results <- dplyr::inner_join(input_metadata,pseudo_core_taxa_table,by="SampleID")

  
  # Run linear regression 
  taxa = rownames(input_data)
  models = lapply(taxa, function(taxon){
    #myformula = paste(taxon, "~ location + Corrected_pedigree + location:Corrected_pedigree")
    lm(input_table_core_taxa_meatadata_results[,taxon] ~ location + Corrected_pedigree + location:Corrected_pedigree,
       data = input_table_core_taxa_meatadata_results)
  })
  anovas <- lapply(models, anova)
  names(anovas) = taxa
  return(anovas)
  
  # input_table_core_taxa_meatadata_results <-lapply(taxa_group_start:ncol(input_table_core_taxa_meatadata_results), 
  #                                                  function(x) lm(input_table_core_taxa_meatadata_results[,x] ~ location + Corrected_pedigree + location:Corrected_pedigree,
  #                                                                 data = input_table_core_taxa_meatadata_results))
  # input_table_core_taxa_annova_results <- lapply(input_table_core_taxa_meatadata_results,anova)
  # return(input_table_core_taxa_annova_results)
}

anovas = lapply(core_taxa, anova_calculation_taxa, input_metadata=metadata)

# annova_calculation_taxa(species_table_core_taxa,G2F_metadata_2019_NO_sub_yellow_stripe_sub,52)
# 
# 
# #run linear 
# species_table_core_taxa_annova_results <- 
# species_table_core_taxa_annova_results
# 
# genus_table_core_taxa_annova_results <- annova_calculation_taxa(genus_table_core_taxa,G2F_metadata_2019_NO_sub_yellow_stripe_sub,52)
# genus_table_core_taxa_annova_results
# 
# family_table_core_taxa_annova_results <-annova_calculation_taxa(family_table_core_taxa,G2F_metadata_2019_NO_sub_yellow_stripe_sub,52)
# family_table_core_taxa_annova_results
# 
# class_table_core_taxa_annova_results <-annova_calculation_taxa(class_table_core_taxa,G2F_metadata_2019_NO_sub_yellow_stripe_sub,52)
# class_table_core_taxa_annova_results
# 
# order_table_core_taxa_annova_results <-annova_calculation_taxa(order_table_core_taxa,G2F_metadata_2019_NO_sub_yellow_stripe_sub,52)
# order_table_core_taxa_annova_results
# 
# phylumn_table_core_taxa_annova_results <-annova_calculation_taxa(phylumn_table_core_taxa,G2F_metadata_2019_NO_sub_yellow_stripe_sub,52)
# phylumn_table_core_taxa_annova_results
# 

#######################
# Plot results
#######################

# 
# heribitility_calculation_taxa_location <- function(anova_result_list){
#   
#   location_heritibility_list <- list()
#   
#   for (x in 1:length(anova_result_list)) {
#     temp_anova_result <- as.data.frame(anova_result_list[x])
#     ssq_residue <- (temp_anova_result$Sum.Sq[4])
#     ssq_location <- (temp_anova_result$Sum.Sq[1])
#     ssq_pedigree <- (temp_anova_result$Sum.Sq[2])
#     ssq_pedigree_location <- (temp_anova_result$Sum.Sq[3])
#     heritbality_location <- (ssq_location/(ssq_residue + ssq_location + ssq_pedigree +ssq_pedigree_location))
#     location_heritibility_list[x] <- heritbality_location
#   }
#   return(location_heritibility_list)
# }
# ###############################################
# heribitility_calculation_taxa_pedigree <- function(anova_result_list){
#   
#   pedigree_heritibility_list <- list()
#   
#   for (x in 1:length(anova_result_list)) {
#     temp_anova_result <- as.data.frame(anova_result_list[x])
#     ssq_residue <- (temp_anova_result$Sum.Sq[4])
#     ssq_location <- (temp_anova_result$Sum.Sq[1])
#     ssq_pedigree <- (temp_anova_result$Sum.Sq[2])
#     ssq_pedigree_location <- (temp_anova_result$Sum.Sq[3])
#     heritbality_pedigree <- (ssq_pedigree/(ssq_residue + ssq_location + ssq_pedigree + ssq_pedigree_location))
#     pedigree_heritibility_list[x] <- heritbality_pedigree
#     
#   }
#   return(pedigree_heritibility_list)
#   
# }
# #####################################################################
# 
# heribitility_calculation_taxa_pedigree_location <- function(anova_result_list){
#   
#   GE_heritibility_list <- list()
#   
#   for (x in 1:length(anova_result_list)) {
#     temp_anova_result <- as.data.frame(anova_result_list[x])
#     ssq_residue <- (temp_anova_result$Sum.Sq[4])
#     ssq_location <- (temp_anova_result$Sum.Sq[1])
#     ssq_pedigree <- (temp_anova_result$Sum.Sq[2])
#     ssq_pedigree_location <- (temp_anova_result$Sum.Sq[3])
#     heritbality_pedigree_location <- (ssq_pedigree_location/(ssq_residue + ssq_location + ssq_pedigree +ssq_pedigree_location))
#     GE_heritibility_list[x] <- heritbality_pedigree_location
#   }
#   return(GE_heritibility_list)
# }


# Helper function to plot results
calc_heritability <- function(anova_results){
  
  taxa = names(anova_results)
  herits = lapply(taxa, function(taxon){
    herit = anova_results[[taxon]] %>% 
      data.frame() %>%
      rename("SS" = "Sum.Sq", "pval" = "Pr..F.") %>%  # Rename
      dplyr::select(SS, pval) %>%
      rownames_to_column("term") %>%
      mutate(herit = SS / sum(SS)) %>% # Calculate heritability
      filter(term != "Residuals") %>% # Remove residuals
      mutate(term = sapply(term, switch, "location"="Environment",
                           "Corrected_pedigree" = "Maize Genotype",
                           "location:Corrected_pedigree" = "GXE"),
             taxon=taxon) %>%
      relocate(taxon)
  }) %>% bind_rows() 
  return(herits)
  
  # # Get SS for each component
  # location_heritibility_list  <- heribitility_calculation_taxa_location(annova_results)
  # pedigree_heritibility_list  <- heribitility_calculation_taxa_pedigree(annova_results)
  # GE_heritibility_list  <- heribitility_calculation_taxa_pedigree_location(annova_results)
  # 
  # location_heritibility_list
  # pedigree_heritibility_list
  # GE_heritibility_list 
  # 
  # #organize the results from list to data table and add taxa name and location 
  # location_heritibility_table <- data.table(heritibility=location_heritibility_list)
  # location_heritibility_table$type <- "Environment"
  # location_heritibility_table$taxa <- taxa_name
  # 
  # #get the top five taxa with highest location heritibility
  # location_heritibility_table <- as.data.frame(lapply(location_heritibility_table,unlist))
  # location_heritibility_table <- location_heritibility_table[order(location_heritibility_table$heritibility,decreasing = TRUE),]
  # 
  # #organize the results from list to data table and add taxa name and location
  # pedigree_heritibility_table <- data.table(heritibility=pedigree_heritibility_list)
  # pedigree_heritibility_table$type <- "Maize genotype" 
  # pedigree_heritibility_table$taxa <- taxa_name
  # 
  # #get the top five taxa with highest location heritibility
  # pedigree_heritibility_table <- as.data.frame(lapply(pedigree_heritibility_table,unlist))
  # pedigree_heritibility_table <- pedigree_heritibility_table[order(pedigree_heritibility_table$heritibility,decreasing = TRUE),]
  # 
  # #organize the results from list to data table and add taxa name and location
  # GE_heritibility_table <- data.table(heritibility=GE_heritibility_list)
  # GE_heritibility_table$type <- "GXE" 
  # GE_heritibility_table$taxa <- taxa_name
  # 
  # GE_heritibility_table <- as.data.frame(lapply(GE_heritibility_table,unlist))
  # GE_heritibility_table <- GE_heritibility_table[order(GE_heritibility_table$heritibility,decreasing = TRUE),]
  # 
  # #get the top five taxa with highest location heritibility
  # location_heritibility_table <- as.data.frame(lapply(location_heritibility_table,unlist))
  # location_heritibility_table <- location_heritibility_table[order(location_heritibility_table$heritibility,decreasing = TRUE),]
  # 
  # #build the heritibility table 
  # heritibility_table <- dplyr::bind_rows(location_heritibility_table,pedigree_heritibility_table)
  # heritibility_table <- dplyr::bind_rows(heritibility_table,GE_heritibility_table)
  # 
  # heritibility_table <- as.data.frame(lapply(heritibility_table, unlist))
  
  ## plot the heritibility of each term on family level core taxa 
  
  # return(heritibility_table)
}

# taxa_level_heritibility <- ggplot(heritibility_table, aes(x=type, y=heritibility,fill = type)) + 
#   geom_violin(trim = FALSE) +  geom_boxplot(width=0.1) + 
#   theme(
#     legend.position='none',
#     axis.title.x = element_blank(),
#     panel.border = element_rect(fill = NA, colour = "grey28"), 
#     panel.background = element_blank(),
#     axis.title.y = element_text(size=30, face="bold"),
#     plot.title = element_text(size=30, face="bold"),
#     axis.text.y = element_text(size=30, face="bold"),
#     axis.text.x = element_text(size=30, face="bold"),
#   )
##############

# Calculate heritability
herits = lapply(anovas, calc_heritability)

# Add taxon levels & combine
for(level in names(herits)){
  herits[[level]] = herits[[level]] %>%
    mutate(level=level) %>%
    relocate(level)
}
herits = bind_rows(herits)

# species_taxa_name <- rownames(species_table_core_taxa$data)
# genus_taxa_name <- rownames(genus_table_core_taxa$data)
# family_taxa_name <- rownames(family_table_core_taxa$data)  
# class_taxa_name <- rownames(class_table_core_taxa$data)  
# order_taxa_name <- rownames(order_table_core_taxa$data)  
# phylumn_taxa_name <- rownames(phylumn_table_core_taxa$data)  

# species_table_core_taxa_heritibility_table <- heritibility_plotting(species_table_core_taxa_annova_results,species_taxa_name)
# genus_table_core_taxa_heritibility_table <- heritibility_plotting(genus_table_core_taxa_annova_results,genus_taxa_name)
# family_table_core_taxa_heritibility_table <- heritibility_plotting(family_table_core_taxa_annova_results,family_taxa_name)
# class_table_core_taxa_heritibility_table <- heritibility_plotting(class_table_core_taxa_annova_results,class_taxa_name)
# order_table_core_taxa_heritibility_table <- heritibility_plotting(order_table_core_taxa_annova_results,order_taxa_name)
# phylumn_table_core_taxa_heritibility_table <- heritibility_plotting(phylumn_table_core_taxa_annova_results,phylumn_taxa_name)

# species_table_core_taxa_heritibility_table$level <- "species"
# genus_table_core_taxa_heritibility_table$level <- "genus"
# family_table_core_taxa_heritibility_table$level <- "family"
# class_table_core_taxa_heritibility_table$level <- "class"
# order_table_core_taxa_heritibility_table$level <- "order"
# phylumn_table_core_taxa_heritibility_table$level <- "phylumn"

# core_taxa_heritibility_table <- dplyr::bind_rows(species_table_core_taxa_heritibility_table,genus_table_core_taxa_heritibility_table)
# core_taxa_heritibility_table <- dplyr::bind_rows(core_taxa_heritibility_table,family_table_core_taxa_heritibility_table)
# core_taxa_heritibility_table <- dplyr::bind_rows(core_taxa_heritibility_table,class_table_core_taxa_heritibility_table)
# core_taxa_heritibility_table <- dplyr::bind_rows(core_taxa_heritibility_table,order_table_core_taxa_heritibility_table)
# core_taxa_heritibility_table <- dplyr::bind_rows(core_taxa_heritibility_table,phylumn_table_core_taxa_heritibility_table)

#######################################################################

# core_taxa_plot_first <- ggplot(core_taxa_heritibility_table[core_taxa_heritibility_table$level %in% c('species','genus','family'),], aes(x=type, y=heritibility,fill = type)) + 
#   geom_violin(trim = FALSE) + geom_boxplot(width=0.1) + ylab("variance explained") + scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) +
#   theme(
#     legend.position='top',
#     axis.title.x = element_blank(),
#     panel.border = element_rect(fill = NA, colour = "grey28"), 
#     panel.background = element_blank(),
#     axis.title.y = element_text(size=30, face="bold"),
#     plot.title = element_text(size=25, face="bold"),
#     axis.text.y = element_text(size=25, face="bold"),
#     axis.text.x = element_blank(),
#   ) + facet_grid(~factor(level, levels=c("species",'genus','family','order',"class",'phylumn'))) +
#   theme(
#     strip.text.x = element_text(
#       size = 35
#     ),
#     strip.text.y = element_text(
#       size = 35
#     )
#   ) + scale_fill_discrete(name = "",)
# 
# core_taxa_plot_first
# 
# #ggsave("/home/hl46161/new_G2F_dada2/exported_table_phyloseq/core_taxa_heritibility_plot_0.25.png", height=15, width=45, device="png")
# 
# 
# core_taxa_plot_second <- ggplot(core_taxa_heritibility_table[core_taxa_heritibility_table$level %in% c('order','class','phylumn'),], aes(x=type, y=heritibility,fill = type)) + 
#   geom_violin(trim = FALSE) + geom_boxplot(width=0.1) + ylab("variance explained") + scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) +
#   theme(
#     legend.position='none',
#     axis.title.x = element_blank(),
#     panel.border = element_rect(fill = NA, colour = "grey28"), 
#     panel.background = element_blank(),
#     axis.title.y = element_text(size=30, face="bold"),
#     plot.title = element_text(size=25, face="bold"),
#     axis.text.y = element_text(size=25, face="bold"),
#     axis.text.x = element_blank(),
#   ) + facet_grid(~factor(level, levels=c("species",'genus','family','order',"class",'phylumn'))) +
#   theme(
#     strip.text.x = element_text(
#       size = 35
#     ),
#     strip.text.y = element_text(
#       size = 35
#     )
#   ) 
# 
# core_taxa_plot_second
# 
# core_taxa_plot <- grid.arrange(core_taxa_plot_first,core_taxa_plot_second)
# 
# core_taxa_plot <- ggarrange(core_taxa_plot_first,core_taxa_plot_second,nrow = 2)
# 
# ggsave("/home/hl46161/new_G2F_dada2/exported_table_phyloseq/core_taxa_heritibility_plot_0.25.png",core_taxa_plot,height=15, width=15, device="png")
# 


# core_taxa_plot <- ggplot(core_taxa_heritibility_table, aes(x=factor(level, levels=c("species",'genus','family','order',"class",'phylumn')), y=heritibility,fill = factor(level, levels=c("species",'genus','family','order',"class",'phylumn')))) + 
#   geom_violin(trim = FALSE) + geom_boxplot(width=0.1) + ylab("variance explained") + scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) +
#   theme(
#     legend.position='top',
#     axis.title.x = element_blank(),
#     panel.border = element_rect(fill = NA, colour = "grey28"), 
#     panel.background = element_blank(),
#     axis.title.y = element_text(size=30, face="bold"),
#     plot.title = element_text(size=25, face="bold"),
#     axis.text.y = element_text(size=25, face="bold"),
#     axis.text.x = element_blank(),
#     legend.text = element_text(size = 35)
#   )  + scale_fill_discrete(name = "",) + facet_grid(~factor(type, levels=c("Environment","Maize genotype","GXE"))) +
#   theme(
#     strip.text.x = element_text(
#       size = 35
#     ),
#     strip.text.y = element_text(
#       size = 35
#     )
#   )

# Plot
herits$level = factor(herits$level, levels=taxa_levels)
core_taxa_plot <- ggplot(herits) +
  aes(x=level, y=herit, fill=term) +
  geom_violin(trim = TRUE) +
  #geom_boxplot(width=0.1) + 
  ylab("variance explained") + 
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey28"), 
    panel.background = element_blank(),
    axis.title.y = element_text(size=30, face="bold"),
    plot.title = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=15,angle=90, face="bold"),
    legend.position="none"
  )  + 
  scale_fill_discrete(name = "",) +
  facet_grid(~term) +
  theme(strip.text.x = element_text(size = 20,face="bold"),
    strip.text.y = element_text(size = 20,face="bold")
  ) + 
  stat_summary(fun = "mean", geom = "crossbar",  width = 0.5,colour = "black")
ggsave(core_taxa_plot, file="3_GxE/3e_taxon_GXE.png", height=8, width=12)

#ggsave("/home/hl46161/new_G2F_dada2/exported_table_phyloseq/core_taxa_heritibility_plot_0.25.png",core_taxa_plot,height=15, width=15, device="png") 
