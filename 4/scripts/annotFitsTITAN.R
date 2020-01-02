rm(list=ls())

library(data.table)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

titan = fread("data/titanCNA_optimalClusterSolutions.csv")

qc_filter = readRDS("data/QC/ovct_filtered_samples.rds")

titan$QC = titan$barcode %in% qc_filter

write.csv(titan, "data/QC/titanCNA_optimalClusterSolutionsQC.csv", row.names=FALSE)
	


# Mannually annotated table format
# -------------------------
titan_clean = fread("data/QC/titanCNA_optimalClusterSolutionsQCmanual_clean.csv")

filter_matrix_geoff = readRDS("data/QC/filter_matrix.rds")
filter_matrix_geoff$barcode = rownames(filter_matrix_geoff)

tab = merge(titan_clean, filter_matrix_geoff, by="barcode")
write.csv(tab, "data/QC/titanCNA_QC.csv", row.names=FALSE)



tab$QC_manual[is.na(tab$QC_manual)] = TRUE

!tab$facets_titan_disagreement & !tab$under_powered & !tab$titan_bug 
!tab$facets_titan_disagreement & !tab$under_powered & !tab$titan_bug & tab$QC_manual == TRUE

sum(!tab$facets_titan_disagreement & !tab$under_powered & !tab$titan_bug & tab$QC_manual == TRUE)