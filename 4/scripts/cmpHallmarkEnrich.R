# Compares GSEA enrichment results for CNA-mRNA, mRNA-purity, and CNA-purity correlation signatures.

rm(list=ls())

library(gplots)
library(RColorBrewer)
library(rjson)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

# Load hallmark color definitions
hallmark_class = fromJSON(file="data/msigDb/hallmark_class.json")

# Load GSEA enrichment result tables
enrich_cna_purity = rbind(
	read.csv("data/GSEA/CNA_purity/GseaPreranked_all_m10000/gsea_report_for_na_pos_1548262603112.csv"),
	read.csv("data/GSEA/CNA_purity/GseaPreranked_all_m10000/gsea_report_for_na_neg_1548262603112.csv")
)

enrich_cna_mrna = rbind(
	read.csv("data/GSEA/CNA_mRNA/GseaPreranked_hallmarks_all_m10000/gsea_report_for_na_pos_1555074277249.csv"),
	read.csv("data/GSEA/CNA_mRNA/GseaPreranked_hallmarks_all_m10000/gsea_report_for_na_neg_1555074277249.csv")
)

enrich_mrna_purity = rbind(
	read.csv("data/GSEA/mRNA_purity/GseaPreranked_hallmarks_all_m10000/gsea_report_for_na_pos_1555074918665.csv"),
	read.csv("data/GSEA/mRNA_purity/GseaPreranked_hallmarks_all_m10000/gsea_report_for_na_neg_1555074918665.csv")
)


# Order and match hallmark terms
enrich_cna_mrna = enrich_cna_mrna[order(enrich_cna_mrna$NES, decreasing=TRUE), ]  # pivot
enrich_cna_purity = enrich_cna_purity[match(enrich_cna_mrna$NAME, enrich_cna_purity$NAME), ]
enrich_mrna_purity = enrich_mrna_purity[match(enrich_cna_mrna$NAME, enrich_mrna_purity$NAME), ]

stopifnot(enrich_cna_mrna$NAME == enrich_cna_purity$NAME)
stopifnot(enrich_cna_mrna$NAME == enrich_mrna_purity$NAME)


# Format NES (normalizaed enrichment score) and FDR matrices
nes_mat = rbind(enrich_cna_mrna$NES, enrich_mrna_purity$NES, enrich_cna_purity$NES)
fdr_mat = rbind(enrich_cna_mrna$FDR, enrich_mrna_purity$FDR, enrich_cna_purity$FDR)

rownames(nes_mat) = c("CNA-mRNA", "mRNA-purity", "CNA-purity")  # used for labeling heatmap

# Format hallmark term names
colnames(nes_mat) = as.character(enrich_cna_mrna$NAME)
colnames(nes_mat) = sapply(strsplit(colnames(nes_mat), "HALLMARK_"), function(x) x[2])
colnames(nes_mat) = gsub("_", " ", colnames(nes_mat))



# Heatmap comparison of enrichment results
pdf("plots/cor_ranking_hallmark_GSEA_enrich_v4.pdf", width=7.35)

# Filter terms to include
fdr = 0.05

include_terms = unique(c(
	as.character(enrich_cna_mrna$NAME[enrich_cna_mrna$FDR < fdr]),
	as.character(enrich_cna_purity$NAME[enrich_cna_purity$FDR < fdr]),
	as.character(enrich_mrna_purity$NAME[enrich_mrna_purity$FDR < fdr])
))


idx = enrich_cna_mrna$NAME %in% include_terms
nes_mat_filter = nes_mat[, idx]

# Evaluate significance
sig = matrix("", nrow=nrow(fdr_mat), ncol=sum(idx))

sig[fdr_mat[, idx] < 0.5] = "."
sig[fdr_mat[, idx] < 0.25] = "*"
sig[fdr_mat[, idx] < 0.1] = "**"
sig[fdr_mat[, idx] < 0.05] = "***"
sig[fdr_mat[, idx] < 0.01] = "****"

heatmap.2(t(nes_mat_filter), trace="none",
	mar=c(12, 28),
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
	Rowv=FALSE, Colv=FALSE,
	cexRow=0.8, cexCol=1,
	cellnote=t(sig),
	notecol="black",
	ylab="GSEA hallmarks",
	xlab="Spearman cor.",
	key.title="",
	key.xlab="Normalized enrichment score (NES)",
	density.info="none"
)
dev.off()


# CNA-mRNA vs mRNA-purity enrichment scatterplot
pdf("plots/scatter_GSEA_enrich_v2.pdf", width=5, height=6)
i = 1
j = 2

pt_col = rep("grey", ncol(fdr_mat))
colors = brewer.pal(9, "Set1")

fdr = 0.05

pt_col[fdr_mat[i, ] < fdr] = colors[1]
pt_col[fdr_mat[j, ] < fdr] = colors[2]
pt_col[fdr_mat[i, ] < fdr & fdr_mat[j, ] < fdr] = colors[3]

legend_idx = fdr_mat[i, ] < fdr | fdr_mat[j, ] < fdr

plot(nes_mat[i, ], nes_mat[j, ],
	xlab=paste0(rownames(nes_mat)[i], " Spearman cor. enrichment (GSEA NES)"),
	ylab=paste0(rownames(nes_mat)[j], " Spearman cor. enrichment (GSEA NES)"),
	pch=21,
	bg=pt_col,
	bty="n"
)

abline(v=0, col="grey", lty=2)
abline(h=0, col="grey", lty=2)

text(x=nes_mat[i, legend_idx],
	y=nes_mat[j, legend_idx],
	labels=colnames(nes_mat)[legend_idx],
	pos=3,
	cex=0.8,
	offset=0.3,
	col=pt_col[legend_idx],
	xpd=TRUE
)

legend("topleft",
	legend=c(
		paste0(rownames(nes_mat)[i], " (FDR < ", fdr, ")"),
		paste0(rownames(nes_mat)[j], " (FDR < ", fdr, ")"),
		paste0(rownames(nes_mat)[i], " & ", rownames(nes_mat)[j], " (FDR < ", fdr, ")")
	),
	pt.bg=colors,
	cex=0.8,
	pch=21
)
dev.off()


# Enrichment scatterplot, colored by hallmark classification
pdf("plots/scatter_GSEA_enrich_v3.pdf", width=5, height=6)
i = 1
j = 2

# Hallmark gene sets color, from Alejandro
hallmark_colors = brewer.pal(9, "Set1")
names(hallmark_colors) = c("RED", "BLUE", "GREEN", "purple", "ORANGE", "YELLOW", "BROWN", "PINK", "GREY")
hallmark_colors = as.list(hallmark_colors)

pt_col = sapply(as.character(enrich_cna_mrna$NAME), function(term) {
	group = hallmark_class[[term]]
	col = hallmark_colors[[group]]
	return(col)
})

fdr = 0.05

pt_pch = rep(21, ncol(fdr_mat))

pt_pch[fdr_mat[i, ] < fdr] = 24  # triangle up
pt_pch[fdr_mat[j, ] < fdr] = 25  # triangle down
pt_pch[fdr_mat[i, ] < fdr & fdr_mat[j, ] < fdr] = 23  # diamond

legend_idx = fdr_mat[i, ] < fdr | fdr_mat[j, ] < fdr

plot(nes_mat[i, ], nes_mat[j, ],
	xlab=paste0(rownames(nes_mat)[i], " Spearman cor. enrichment (GSEA NES)"),
	ylab=paste0(rownames(nes_mat)[j], " Spearman cor. enrichment (GSEA NES)"),
	pch=pt_pch,
	bg=pt_col,
	bty="n"
)

abline(v=0, col="grey", lty=2)
abline(h=0, col="grey", lty=2)

text(x=nes_mat[i, legend_idx],
	y=nes_mat[j, legend_idx],
	labels=colnames(nes_mat)[legend_idx],
	pos=3,
	cex=0.6,
	offset=0.3,
	col=pt_col[legend_idx],
	xpd=TRUE
)

legend("topleft",
	legend=c(
		paste0(rownames(nes_mat)[i], " (FDR < ", fdr, ")"),
		paste0(rownames(nes_mat)[j], " (FDR < ", fdr, ")"),
		paste0(rownames(nes_mat)[i], " & ", rownames(nes_mat)[j], " (FDR < ", fdr, ")"),
		"Oncogenic", "Stromal", "Immune", "Cellular stress", "Other"
	),
	pt.bg=c(
		rep("grey", 3),
		unlist(hallmark_colors[c("RED", "BLUE", "GREEN", "purple", "GREY")])),
	cex=0.8,
	pch=c(
		c(24, 25, 23),
		rep(21, 5))
)
dev.off()
