rm(list=ls())

library(RColorBrewer)
setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

# Load TitanCNA results table
titan = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")

# Load ASCAT
ascat = fread("data/ASCAT/4_Summary.txt")
ascat$barcode = sapply(strsplit(ascat$name, "_"), function(x) x[1])

# ascat$CCF[ascat$CCF == 1] = NA   # remove possible outlier
stopifnot(all(ascat$barcode == titan$barcode))

# Compare ASCAT with TITAN, purity and ploidy
pdf("plots/ASCAT_TITAN_ploidy_purity.pdf", width=8, height=4.5)
par(mfrow=c(1, 2))
scatterPlot(ascat$Ploidy, titan$ploidy, xlab="ASCAT ploidy", ylab="TITAN ploidy",
	alternative="greater")
scatterPlot(ascat$CCF, titan$purity, xlab="ASCAT purity", ylab="TITAN purity",
	method="kendall",
	alternative="greater")
dev.off()

# scatterPlot(ascat$Ploidy, ascat$CCF)
# cor.test(ascat$Ploidy, ascat$CCF)
