rm(list=ls())

library(data.table)
library(RColorBrewer)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")


# Load ssGSEA hallmark enrichment results
hallmarks = fread("data/HallmarksStromaImmune_NES.txt")
terms = hallmarks$Term
hallmarks = data.matrix(hallmarks[, -1])
rownames(hallmarks) = terms

# Load Ki67 IF data
if_frac = fread("data/Ki67/TableS1_TxNaiveData.csv")

if_frac_matched = if_frac[match(colnames(hallmarks), if_frac$RNA_WELL), ]


pathway = "HALLMARK_G2M_CHECKPOINT"

nes = hallmarks[which(rownames(hallmarks) == pathway), ]

pdf(paste0("plots/ssGSEA_", pathway, "_Ki67.pdf"), width=4.0, height=4.5)
scatterPlot(if_frac_matched$Ki67, nes,
	xlab="Ki67 (% positive cells)",
	ylab=paste0(pathway, " (NES, mRNA ssGSEA)"),
	method="kendall",
	bg=brewer.pal(9, "Pastel1")[1],
	pch=21)
dev.off()


# idx = nes > 0.45
# sum(idx)

# scatterPlot(if_frac_matched$Ki67[!idx], nes[!idx], method="spearman")
# if_frac_matched[idx, ]
# scatterPlot(if_frac_matched$Ki67, nes, method="kendall")
