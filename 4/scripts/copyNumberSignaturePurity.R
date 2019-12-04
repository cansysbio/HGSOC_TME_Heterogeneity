# Analysis of association between TitanCNA purity and copy-number signatures.

rm(list=ls())

library(data.table)
library(RColorBrewer)
library(gmodels)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

titan = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")

patient_colors = loadPatientColors(titan)

# cn_sig = readRDS("data/signatureExposure/miller_cnsignatures.rds")  # TITAN signatures
cn_sig = readRDS("data/signatureExposure/copywriteR_ovct_sigs_20191126.rds")  # copywriteR signatures

# Fraction of samples with more than 20% exposure by signature
apply(cn_sig > 0.2, 1, mean)

apply(cn_sig > 0.05, 1, mean)

# Mean and confidence intervals 
apply(cn_sig, 1, ci)

# Fraction of samples with particular highest (maximum) signature
prop.table(table(apply(cn_sig, 2, which.max)))

cn_sig[, grep("RG4", colnames(cn_sig))]

cn_sig[, grep("RG6", colnames(cn_sig))]

# Test if tumor labels are aligned
stopifnot(all(titan$barcode == colnames(cn_sig)))


pdf("plots/cn_signature_scatter.pdf", height=2.7, width=5.5)
pt_col = patient_colors[as.integer(titan$patient_id)]

par(mfrow=c(1, 3))
for (i in c(1, 6, 4)) {
	scatterPlot(cn_sig[i, ], titan$purity,
		xlab=rownames(cn_sig)[i],
		ylab="Tumour purity (TitanCNA)",
		pch=21,
		bg=pt_col
	)
}
dev.off()

# Stacked barplot ordered by purity
# pdf("plots/cn_signature_purity_order.pdf", height=3.5, width=9)
pdf("plots/cn_signature_purity_order_copywriteR.pdf", height=3.5, width=9)

idx = order(titan$purity, decreasing=TRUE)

sig_colors = brewer.pal(8, "Set2")

bp = barplot(cn_sig[, idx],
	col=sig_colors,
	las=2,
	legend.text=row.names(cn_sig),
	ylab="Exposure, purity",
	xlab="Tumours",
	cex.names=0.8
)

lines(bp, titan$purity[idx], col=brewer.pal(9, "Set1")[1], lwd=2)

points(bp, rep(-0.05, length(bp)),
	pch=21,
	bg=pt_col[idx],
	xpd=TRUE)
dev.off()
