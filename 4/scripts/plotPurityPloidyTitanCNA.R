# Plots for the default copy-number alteration models from TitanCNA
# Only considers summary statistics: average ploidy and purity
rm(list=ls())

library(RColorBrewer)
library(stringr)
library(data.table)


setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
# source("lib/plots.R")

# Load TitanCNA results table
opt_clust = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")

colors = loadPatientColors(opt_clust)


pdf("plots/titanCNA_ploidy_purity.pdf", width=5.5, height=5.5)

pts_col = colors[as.integer(opt_clust$patient_id)]
pts_cex = opt_clust$numClust

cor_test = cor.test(opt_clust$ploidy, opt_clust$purity)

plot(opt_clust$ploidy, opt_clust$purity,
	main="Copy-number alterations",
	xlab="Ave. tumor ploidy",
	ylab="Tumor purity",
	bty="n",
	pch=21,
	bg=pts_col,
	cex=pts_cex
)

legend("bottomright", legend=levels(opt_clust$patient_id), pch=21, pt.bg=colors)

legend("topleft", legend=c(
	paste0("r=", format(cor_test$estimate, digits=3)),
		paste0("P=", format(cor_test$p.value, digits=3))
	),
	bty="n"
)
dev.off()
