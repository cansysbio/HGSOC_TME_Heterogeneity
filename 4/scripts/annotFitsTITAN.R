rm(list=ls())

library(data.table)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

titan = fread("data/titanCNA_optimalClusterSolutions.csv")

qc_filter = readRDS("data/QC/ovct_filtered_samples.rds")

titan$QC = titan$barcode %in% qc_filter

write.csv(titan, "data/QC/titanCNA_optimalClusterSolutionsQC.csv", row.names=FALSE)
	