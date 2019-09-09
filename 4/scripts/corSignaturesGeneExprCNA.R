# Gene-specific correlations with copy-numbers and purity estimates.
# 
# Writes gene signatures files for use with GSEA.

rm(list=ls())

library(RColorBrewer)
library(data.table)

library(gplots)
library(viridis)
library(devtools)  # for installing heatmap.3
library(metap)  # for combining p-values
library(GSA)
library(rjson)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

# load hallmark gene sets from msigDB
hallmark = GSA.read.gmt("data/msigDB/h.all.v6.2.symbols.gmt")

names(hallmark$genesets) = hallmark$geneset.names

# Load hallmark color definitions
hallmark_class = fromJSON(file="data/msigDB/hallmark_class.json")

hallmark_colors = brewer.pal(9, "Set1")
names(hallmark_colors) = c("RED", "BLUE", "GREEN", "purple", "ORANGE", "YELLOW", "BROWN", "PINK", "GREY")
hallmark_colors = as.list(hallmark_colors)


# Load CNA count matrix
# ------------------------------------
cna = loadCNA("data/cna.csv")  # TITAN
# cna = loadCNA("data/cna_FACETS.csv")

cna_mat = data.matrix(cna[, c(-1, -2)])
rownames(cna_mat) = cna$hgnc_symbol

colnames(cna_mat)[colnames(cna_mat) == "RG13T12"] = "RG13T122"  # Reintroduce typo for consistency with mRNA annotation

# Load gene expression matrix
expr = fread("data/OVCTp_log2exp_loess_norm.txt")

# expr_samples = loadSampleAnnot("data/OVCT_Tnaive_WES_labels.txt")
expr_samples = loadSampleAnnot(file_path="data/TreatmentNaive_SampleLabels_WESTumourCellularity_mRNAtumourCellularity_MAPPINGS.csv")

# Load Titan table
titan = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")

# expression matrix
expr_mat = data.matrix(expr[, -1])

colnames(expr_mat) = expr_samples$wes_label[match(colnames(expr_mat), expr_samples$Well)]
rownames(expr_mat) = expr$Hugo_Symbol

message("Matching samples: ", sum(colnames(expr_mat) %in% colnames(cna_mat)))
message("Mismatches: ",
	paste(colnames(expr_mat)[!colnames(expr_mat) %in% colnames(cna_mat)], collapse=", "))

# Match expression matrix to cna_mat
expr_mat_match = expr_mat[
	match(rownames(cna_mat), rownames(expr_mat)),
	match(colnames(cna_mat), colnames(expr_mat))
]

# Include genes where we have at least one measurement from microarray data
include_genes = apply(expr_mat_match, 1, function(row) {
	any(!is.na(row))
})

cna_mat_filter = cna_mat[include_genes, ]
expr_mat_match_filter = expr_mat_match[include_genes, ]

# Correlate CNA with gene expression
cor_tests = lapply(1:nrow(cna_mat_filter), function(i) {
	out = tryCatch({
		cor.test(cna_mat_filter[i, ], expr_mat_match_filter[i, ], method="spearman")
	}, error=function(err) {
		warning(err)
		return(NA)
	})

	return(out)
})
names(cor_tests) = rownames(cna_mat_filter)


# Collect correlation statistics
cna_mrna_stats = data.frame(
	gene_symbol=rownames(cna_mat_filter),
	getCorStats(cor_tests)
)

cna_mrna_stats$mean_cna = apply(cna_mat_filter, 1, mean, na.rm=TRUE)

# Remove missing correlation coefficient entries
cna_mrna_stats = cna_mrna_stats[!is.na(cna_mrna_stats$r), ]

# Correct for multiple hypotheses
cna_mrna_stats$padj = p.adjust(cna_mrna_stats$p)

cna_mrna_stats = cna_mrna_stats[order(cna_mrna_stats$r, decreasing=TRUE), ]

# Write signature for running GSEA
write.table(cna_mrna_stats[, c("gene_symbol", "r")],
	"data/GSEA/CNA_mRNA_spearman_rank.rnk",
	# "data/GSEA/CNA_mRNA_spearman_rank_FACETS.rnk",
	quote=FALSE,
	sep="\t",
	row.names=FALSE)	

write.csv(cna_mrna_stats,
	file="data/corSig/CNA_mRNA_cor_spearman.csv",
	# file="data/corSig/CNA_mRNA_cor_spearman_FACETS.csv",
	row.names=FALSE)


# CNA-mRNA volcano plot
pdf("plots/cna_mRNA_volcano_v2.pdf", width=3.8, height=4)
cna_mrna_stats = cna_mrna_stats[order(cna_mrna_stats$p), ]

fdr = 0.05
sum(cna_mrna_stats$padj < fdr)

n_label = 10  # top points to label
highlight_col = brewer.pal(9, "Set1")[1]
pts_col = rep("grey", nrow(cna_mrna_stats))
pts_col[cna_mrna_stats$padj < fdr] = highlight_col
plot(cna_mrna_stats$r, -log10(cna_mrna_stats$p),
	col=pts_col,
	pch=16,
	cex=0.7,
	xlab="Spearman cor.",
	ylab=expression("-log"[10] * " P")
)

text(cna_mrna_stats$r[1:n_label], -log10(cna_mrna_stats$p[1:n_label]),
	labels=cna_mrna_stats$gene_symbol[1:n_label],
	pos=2,  # to the left
	cex=0.8
)
dev.off()


# CNA-purity signature
# ------------------------------------------------
cna_purity_tests = lapply(1:nrow(cna_mat), function(i) {
	out = tryCatch({
		cor.test(cna_mat[i, ], titan$purity, method="spearman")
	}, error=function(err) {
		warning(err)
		return(NA)
	})

	return(out)
})
names(cna_purity_tests) = rownames(cna_mat)

cna_purity_stats = data.frame(
	gene_symbol=rownames(cna_mat),
	getCorStats(cna_purity_tests)
)

cna_purity_stats$padj = p.adjust(cna_purity_stats$p)

# Remove missing
cna_purity_stats = cna_purity_stats[!is.na(cna_purity_stats$r), ]

cna_purity_stats = cna_purity_stats[order(cna_purity_stats$r, decreasing=TRUE), ]

write.table(cna_purity_stats[, c("gene_symbol", "r")],
	"data/GSEA/CNA_purity_spearman_rank.rnk",
	# "data/GSEA/CNA_purity_spearman_rank_FACETS.rnk",
	quote=FALSE,
	sep="\t",
	row.names=FALSE)	


# mRNA-purity signature
# -------------------------------

# Check if expression matrix matches
stopifnot(all(colnames(expr_mat_match) == titan$barcode, na.rm=TRUE))

mrna_purity_tests = lapply(1:nrow(expr_mat_match), function(i) {
	tryCatch({
		cor.test(expr_mat_match[i, ], titan$purity, method="spearman")
	}, error=function(err) {
		warning(err)
		return(NA)
	})
})
names(mrna_purity_tests) = rownames(expr_mat_match)

mrna_purity_stats = getCorStats(mrna_purity_tests)

# Remove missing
mrna_purity_stats = mrna_purity_stats[!is.na(mrna_purity_stats$r), ]

# Order by correlation coefficients
mrna_purity_stats = mrna_purity_stats[order(mrna_purity_stats$r, decreasing=TRUE), ]

mrna_purity_stats$padj = p.adjust(mrna_purity_stats$p)

write.table(mrna_purity_stats[, c("gene_symbol", "r")],
	"data/GSEA/mRNA_purity_spearman_rank.rnk",
	quote=FALSE,
	sep="\t",
	row.names=FALSE)	

write.csv(mrna_purity_stats,
	file="data/corSig/mRNA_purity_cor_spearman.csv",
	# file="data/corSig/CNA_mRNA_cor_spearman_FACETS.csv",
	row.names=FALSE)


# mRNA-CNA, mRNA-purity scatter
# -------------------------------------------------------------
cna_mrna_stats_matched = cna_mrna_stats[
	match(mrna_purity_stats$gene_symbol, cna_mrna_stats$gene_symbol), ]

stopifnot(all(cna_mrna_stats_matched$gene_symbol == mrna_purity_stats$gene_symbol, na.rm=TRUE))


# Gene-level mRNA and CNA vs purity scatter plot
pdf("plots/mRNA_CNA_purity_cor_coeff_scatter_v2.pdf", width=4.7, height=5)

pts_col = rep("grey", nrow(cna_mrna_stats_matched))

for (pathway in names(hallmark$genesets)) {
	genes = hallmark$genesets[[pathway]]

	# Get color of hallmark gene set
	color_class = hallmark_class[[pathway]]
	col = hallmark_colors[[color_class]]

	pts_col[cna_mrna_stats_matched$gene_symbol %in% genes] = col
}

r_cutoff = 0.5

idx = cna_mrna_stats_matched$gene_symbol %in% unlist(hallmark$genesets)

cor_test = cor.test(cna_mrna_stats_matched$r[idx], mrna_purity_stats$r[idx])

plot(cna_mrna_stats_matched$r[idx], mrna_purity_stats$r[idx],
	xlab="CNA-mRNA Spearman cor.",
	ylab="mRNA-purity Spearman cor.",
	main=paste0(
		"Hallmark genes",
		", r=", format(cor_test$estimate, digits=3),
		", P=", format(cor_test$p.value, digits=3),
		", n=", sum(idx)),
	cex.main=0.8,
	col=pts_col[idx],
	lwd=0.5,
	pch=16,
	cex=0.5,
	bty="n"
)

abline(v=r_cutoff, col="grey", lty=2)
abline(v=-r_cutoff, col="grey", lty=2)
abline(h=r_cutoff, col="grey", lty=2)
abline(h=-r_cutoff, col="grey", lty=2)

# idx_label = cna_mrna_stats_matched$gene_symbol %in% unlist(hallmark$genesets) & (cna_mrna_stats_matched$r > r_cutoff) & (abs(mrna_purity_stats$r) > r_cutoff)

idx_label = cna_mrna_stats_matched$gene_symbol %in% unlist(hallmark$genesets) & (abs(cna_mrna_stats_matched$r) > r_cutoff) & (abs(mrna_purity_stats$r) > r_cutoff)

text(cna_mrna_stats_matched$r[idx_label], mrna_purity_stats$r[idx_label],
	labels=cna_mrna_stats_matched$gene_symbol[idx_label],
	cex=0.8,
	offset=0.0,
	pos=3)

legend("topleft",
	legend=c("Oncogenic", "Stromal", "Immune", "Cellular stress", "Other"),
	col=unlist(hallmark_colors[c("RED", "BLUE", "GREEN", "purple", "GREY")]),
	cex=0.8,
	pch=16
)

dev.off()
