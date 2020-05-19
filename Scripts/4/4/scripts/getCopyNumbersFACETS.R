rm(list=ls())

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")

# Data dir containing FASCETS output files
data_dir = "~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/data/pier/New Segment Files_FACETS"

# Load all copy number segments
file_names = list.files(data_dir, "*cncf.txt")
segs = lapply(file_names, function(file_name) {
	fread(file.path(data_dir, file_name))
})
names(segs) = sapply(strsplit(file_names, "_"), function(x) x[1])

# Format segment tables
segs_format = lapply(segs, function(s) {
	s$chrom = paste0("chr", s$chrom)
	s$start = s$loc.start
	s$end = s$loc.end
	return(s)
})

cna_mat = segmentGeneFeatures(segs_list=segs_format, feature="tcn.em")  # total copy number

write.csv(cna_mat, file="data/cna_FACETS.csv", row.names=FALSE)
