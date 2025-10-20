#! /usr/bin/Rscript

# Update the table stats with the new GxE table

source("1c_CalculateOtuTableStats.r")


gxe = readRDS("3_GxE/3a_asv_table.rarefied.filt_for_gxe.rds")
stats(gxe, "GxE_subset")
