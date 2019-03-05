# log ratio test taking into account the compositional nature of the copy-number signature fractions.
rm(list=ls())

library(data.table)
# library(RColorBrewer)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

titan = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")

patient_colors = loadPatientColors(titan)

cn_sig = readRDS("data/miller_cnsignatures.rds")

# Test if tumor labels are aligned
stopifnot(all(titan$barcode == colnames(cn_sig)))

# One vs all log ratios
log_ratio = sapply(1:nrow(cn_sig), function(k) {
	log10(cn_sig[k, ] / apply(cn_sig[-k, ], 2, sum))
})

log_ratio[is.infinite(log_ratio)] = NA
log_ratio = t(log_ratio)


pdf("plots/copyNumberSignaturesCompositional.pdf", height=2.8)
par(mfrow=c(1, 3))
for (k in c(1, 6, 4)) {
	scatterPlot(log_ratio[k, ], titan$purity,
		xlab=paste0("log s", k, " / s"),
		ylab="Purity (TitanCNA)"
	)
}
dev.off()
