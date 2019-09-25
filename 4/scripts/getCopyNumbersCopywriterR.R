rm(list=ls())

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")

dirs = list.dirs("data/copywriteR", recursive=FALSE)

# Load all genomic segments
segs = lapply(dirs, function(dir) {
	env = new.env()
	load(paste0(dir, "/CNAprofiles/segment.RData"), env)
	return(env$segment.CNA.object)
})
names(segs) = basename(dirs)


segs_format = lapply(segs, function(x) {
	s = x$output

	# Exclude 'none' comparisons 
	idx = grep(".vs.none", s$ID, invert=TRUE)

	s = s[idx, ]

	s$chrom = paste0("chr", s$chrom)

	s$start = floor(s$loc.start)
	s$end = ceiling(s$loc.end)

	# Exclude large average segments for each chromosome
	exclude = s$seg.mean == 0.0

	message("Excluding zero mean segments: ", sum(exclude))

	s = s[!exclude,]

	return(s)
})

cna_mat = segmentGeneFeatures(segs_list=segs_format, feature="seg.mean")

write.csv(cna_mat, file="data/cna_copywriteR.csv", row.names=FALSE)

# cna_mat[which(cna_mat$hgnc_symbol == "MYC"), ]
# cna_mat[which(cna_mat$hgnc_symbol == "KRAS"), ]
# cna_mat[which(cna_mat$hgnc_symbol == "CDH1"), ]
