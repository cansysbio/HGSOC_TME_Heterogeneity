rm(list=ls())

library(data.table)
library(GSA)
library(gplots)
library(RColorBrewer)
library(WGCNA)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

# Load CNA
cna = loadCNA("data/cna.csv")

cna_mat = data.matrix(cna[, c(-1, -2)])
rownames(cna_mat) = cna$hgnc_symbol

colnames(cna_mat)[colnames(cna_mat) == "RG13T12"] = "RG13T122"  # Reintroduce typo for consistency with mRNA annotation


# load hallmark gene sets from msigDB
hallmark = GSA.read.gmt("data/msigDB/h.all.v6.2.symbols.gmt")

names(hallmark$genesets) = hallmark$geneset.names


# Load Titan table
titan = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")



# Analysis of MYC and MYC target CNA mutual exclusivity
# ------------------------------------------
#genes = hallmark$genesets$HALLMARK_MYC_TARGETS_V2
genes = hallmark$genesets$HALLMARK_MYC_TARGETS_V1

myc_genes = c("MYC", "MYCN", "MYCL")

# genes = hallmark$genesets$HALLMARK_MYC_TARGETS_V1

genes = c(genes, myc_genes)

idx = rownames(cna_mat) %in% genes


cp_colors = c(
	rev(brewer.pal(5, "Blues")[c(4, 5)]), # 0, 1
	rgb(220, 220, 220, maxColorValue=255),
	brewer.pal(8, "Reds")[c(-1, -2)]
)

colors = brewer.pal(9, "Set1")
colors = c("white", "black")
row_cols = colors[as.integer(rownames(cna_mat[idx, ]) %in% myc_genes) + 1]

heatmap.2(cna_mat[idx, ],
	trace="none",
	col=cp_colors,
	denscol="black",
	RowSideColors=row_cols
)


idx = rownames(cna_mat) %in% genes

cor_tests = corAndPvalue(t(cna_mat[idx, ]), t(cna_mat[rownames(cna_mat) %in% myc_genes, ]))
head(sort(cor_tests$cor))
cor_tests$cor




stopifnot(colnames(cna_mat) == titan$barcode)

# pdf("plots/mutualExclusivityCna.pdf", width=3.3)
# pdf("plots/mutualExclusivityCna.pdf", width=7.5, height=3.0)
pdf("plots/mutualExclusivityCna.pdf", width=6.3, height=2.5)
par(mfrow=c(1, 4))

colors = loadPatientColors(titan)
pts_col = colors[as.integer(titan$patient_id)]

gene1 = "MYC"
gene2 = "MRTO4"

j = which(rownames(cna_mat) == gene1)
i = which(rownames(cna_mat) == gene2)
scatterPlotJit(cna_mat[i, ], cna_mat[j, ],
	corner="topright",
	col=pts_col,
	jit_amount=0.2,
	xlab=paste(rownames(cna_mat)[i], "CN"),
	ylab=paste(rownames(cna_mat)[j], "CN")
)

abline(v=4.5, lty=3, col="grey")
abline(h=4.5, lty=3, col="grey")

gene1 = "MYC"
gene2 = "MAP3K6"

j = which(rownames(cna_mat) == gene1)
i = which(rownames(cna_mat) == gene2)
scatterPlotJit(cna_mat[i, ], cna_mat[j, ],
	corner="topright",
	col=pts_col,
	jit_amount=0.2,
	xlab=paste(rownames(cna_mat)[i], "CN"),
	ylab=paste(rownames(cna_mat)[j], "CN")
)

abline(v=4.5, lty=3, col="grey")
abline(h=4.5, lty=3, col="grey")


gene1 = "MYC"
gene2 = "HNRNPR"

j = which(rownames(cna_mat) == gene1)
i = which(rownames(cna_mat) == gene2)
scatterPlotJit(cna_mat[i, ], cna_mat[j, ],
	corner="topright",
	col=pts_col,
	jit_amount=0.2,
	xlab=paste(rownames(cna_mat)[i], "CN"),
	ylab=paste(rownames(cna_mat)[j], "CN")
)

abline(v=4.5, lty=3, col="grey")
abline(h=4.5, lty=3, col="grey")


dev.off()


# Mean MYC target copy numbers
idx = rownames(cna_mat) %in% hallmark$genesets$HALLMARK_MYC_TARGETS_V1
mean_myc_target_v1 = apply(cna_mat[idx, ], 2, mean, na.rm=TRUE)

idx = rownames(cna_mat) %in% hallmark$genesets$HALLMARK_MYC_TARGETS_V2
mean_myc_target_v2 = apply(cna_mat[idx, ], 2, mean, na.rm=TRUE)


pdf("plots/mutualExclusivityMeanMyc.pdf", height=5, width=6.5)
colors = loadPatientColors(titan)
pts_col = colors[as.integer(titan$patient_id)]

par(mfcol=c(2, 3))
for (gene in myc_genes) {
	i = which(rownames(cna_mat) == gene)
	scatterPlot(cna_mat[i, ], mean_myc_target_v1,
		pch=21,
		bg=pts_col,
		xlab=paste(rownames(cna_mat)[i], "(CN)"),
		ylab="MYC target v1 (mean CN)"
	)
	scatterPlot(cna_mat[i, ], mean_myc_target_v2,
		pch=21,
		bg=pts_col,
		xlab=paste(rownames(cna_mat)[i], "(CN)"),
		ylab="MYC target v2 (mean CN)"
	)
}
dev.off()
