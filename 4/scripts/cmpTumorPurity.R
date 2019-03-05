rm(list=ls())

library(RColorBrewer)
library(data.table)

library(gplots)
library(viridis)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

# Load TitanCNA results table
titan = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")


# Load mRNA purity estimates from Alejandro
# sample annotation
sample_annot = loadSampleAnnot("data/OVCT_Tnaive_WES_labels.txt")

conv_scores = loadConvScores("data/OVCTp_log2exp_loess_norm_estimate_score.txt")

# combine purity scores with sample annotation data
conv = merge(conv_scores, sample_annot, by="exp_label")

# Match conv to titan results
conv_match = conv[match(titan$barcode, conv$tumor), ]


pdf("plots/purity_comparison_mRNA_CNA.pdf", height=2.5, width=6)

colors = brewer.pal(9, "Set1")
par(mfrow=c(1, 3))

scatterPlot(conv_match$TumorPurity, titan$purity,
	xlab="Tumor purity (mRNA deconvolution)",
	ylab="Tumor purity (TitanCNA)",
	alternative="greater"
)


scatterPlot(conv_match$ImmuneScore, titan$purity,
	xlab="Immune score",
	ylab="Tumor purity (TitanCNA)",
	alternative="less",
	corner="topright",
	col=colors[2]
)

scatterPlot(conv_match$StromalScore, titan$purity,
	xlab="Stromal score",
	ylab="Tumor purity (TitanCNA)",
	alternative="less",
	corner="topright",
	col=colors[3]
)
dev.off()


# Confusion matrix when classifying tumors into cold and hot using median value
# -------------------------------------------------------
median_conf_mat = table(
	paste("RNA", median(conv_match$TumorPurity, na.rm=TRUE) > conv_match$TumorPurity),
	paste("CNA", median(titan$purity, na.rm=TRUE) > titan$purity)
)

pdf("plots/median_tumor_purity_confusion_matrix_mRNA_CNA.pdf", width=6, height=6)
heatmap.2(median_conf_mat, trace="none",
	main="Confusion matrix",
	xlab="CNA tumor purity",
	ylab="mRNA tumor purity",
	Rowv=FALSE, Colv=FALSE,
	col=colorRampPalette(brewer.pal(9, "Greens"))(100),
	mar=c(20, 22),
	cellnote=median_conf_mat,
	notecol="black",
	cexCol=1.0, cexRow=1.0
)
dev.off()
