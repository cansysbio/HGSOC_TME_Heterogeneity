# Contains code for both 10kb and 30kb. Manual selections
# 
rm(list=ls())

library(data.table)
library(TitanCNA)
library(stringr)
library(IRanges)
library(GenomicRanges)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/plots.R")
source("lib/parse.R")

results_dir = "/Users/koplev01/mnt/mmlab_data/projects/OVCT/titanCNA/TitanCNA/scripts/snakemake"

# opt_clust = fread("data/titanCNA_optimalClusterSolution_TCGA.txt")
# opt_clust = fread("data/titanCNA_optimalClusterSolution_TCGA30kb.txt")
opt_clust = fread("data/titanCNA_optimalClusterSolution_TCGA30kb_extend.txt")

scatterPlot(opt_clust$ploidy, opt_clust$purity,
	xlab="Ave. tumor ploidy", ylab="Tumor purity")


# # Load all optimal fits
# d = loadTitanResultsEnv(
# 	paste0("TCGA_", opt_clust$path),  # output folder was manually renamed
# 	results_dir)

# # Load all optimal fits
# d = loadTitanResultsEnv(
# 	paste0("TCGA30kb_", opt_clust$path),  # output folder was manually renamed
# 	results_dir)


d = loadTitanResultsEnv(
	paste0("TCGA30kb_extend_", opt_clust$path),  # output folder was manually renamed
	results_dir)



# evaluate segments of copy-number calls from TitanCNA
segs = lapply(1:nrow(opt_clust), function(i) {
	message(i)
	s = outputTitanSegments(d[[i]]$results, id=opt_clust$id[i], d[[i]]$convergeParams)
	return(s)
})


# Calculate non-integer copy numbers
segs = lapply(segs, function(s) {
	# Assumes that remaining cancer cells for subclonal calls are diploid
	s$Copy_Number_Noninteger = s$Cellular_Frequency * s$Copy_Number + 2 * (1 - s$Cellular_Frequency)

	# (Rare) segments without clonal calls, use copy number calls
	idx = is.na(s$Clonal_Cluster)
	s$Copy_Number_Noninteger[idx] = s$Copy_Number[idx]
	message("Calls wihtout subclones: ", sum(idx))
	return(s)
})

# save(segs, file="data/segs_TCGA.RData")
# save(segs, file="data/segs_TCGA30kb.RData")
save(segs, file="data/segs_TCGA30kb_extend.RData")



data.frame(segs[[7]])

# Check if ranges are overlapping
k = 5
segs[[k]]

# Convert to GRanges object
gr = makeGRangesFromDataFrame(segs[[k]],
	keep.extra.columns=TRUE,
	start.field="Start_Position.bp.",
	end.field="End_Position.bp.",
)


reduce(gr)

genome_coverage = coverage(gr)

for (i in 1:length(genome_coverage)) {
	print(table(coverage(gr)[[i]]))
}
