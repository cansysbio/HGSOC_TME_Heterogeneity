rm(list=ls())

library(stringr)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/plots.R")

titan_cna = fread("data/cna.csv")
titan_cna = data.frame(titan_cna)

facets_cna = fread("data/cna_FACETS.csv")
facets_cna = data.frame(facets_cna)

# Match column names
facets_cna = facets_cna[, match(colnames(titan_cna), colnames(facets_cna))]

# Match by entrez IDs
facets_cna = facets_cna[match(titan_cna$entrez_id, facets_cna$entrez_id), ]

# Test if gene symbols are aligned
stopifnot(facets_cna$gene_symbol == titan_cna$gene_symbol)
stopifnot(colnames(facets_cna) == colnames(titan_cna))
stopifnot(dim(facets_cna) == dim(titan_cna))


# Variance explained per sample
# -----------------------------
R2 = sapply(3:ncol(titan_cna), function(i) {
	cor(titan_cna[, i], facets_cna[, i], use="pairwise.complete.obs")^2
})
names(R2) = colnames(titan_cna)[3:ncol(titan_cna)]

pdf("plots/cna_TITAN_FACETS_var_explained.pdf", width=2.5)
boxplot(R2,
	# frame=FALSE,
	outline=FALSE,
	ylim=range(R2),
	ylab="Per-sample variance explained R2 (TITAN vs FACETS gene copy numbers)"
)
points(
	jitter(rep(1, length(R2)), amount=0.1),
	R2,
	pch=21,
	bg="grey"
)
dev.off()

# Join principal component analysis
# -----------------------------------

titan_cna_mat = data.matrix(titan_cna[, c(-1, -2)])
colnames(titan_cna_mat) = paste0("TITAN_", colnames(titan_cna_mat))

facets_cna_mat = facets_cna[, c(-1, -2)]
colnames(facets_cna_mat) = paste0("FACETS_", colnames(facets_cna_mat))

cna_comb = cbind(titan_cna_mat, facets_cna_mat)

cna_comb[is.na(cna_comb)] = 2  # crude assumption, imputing to diploidy 



# Get sample annotation data frame
samples = str_split_fixed(colnames(cna_comb), "_", 2)
samples = data.frame(samples)
colnames(samples) = c("method", "id")


pca = prcomp(t(cna_comb), scale=TRUE, center=TRUE)

pdf("plots/cna_TITAN_FACETS_PCA.pdf", width=4, height=4.5)
pts_pch = rep(NA, nrow(samples))
pts_pch[samples$method == "TITAN"] = 21
pts_pch[samples$method == "FACETS"] = 24

summary(pca)
plot(pca$x[, 1], pca$x[, 2],
	main="Gene copy numbers",
	type="n",
	xlab=paste0(
		"PC1 (",
		format(summary(pca)$importance[2, 1] * 100, digits=3), "%)"
	),
	ylab=paste0(
		"PC2 (",
		format(summary(pca)$importance[2, 2] * 100, digits=3), "%)"
	)
)
# text(pca$x[, 1], pca$x[, 2], labels=colnames(cna_comb), cex=0.5)

titan_idx = samples$method == "TITAN"
facets_idx = samples$method == "FACETS"

segments(
	pca$x[titan_idx, 1],
	pca$x[titan_idx, 2],
	pca$x[facets_idx, 1],
	pca$x[facets_idx, 2]
)

points(pca$x[, 1], pca$x[, 2],
	pch=pts_pch,
	bg="grey"
)

legend("bottomright",
	legend=c("TITAN", "FACETS"),
	pch=c(21, 24),
	pt.bg="grey"
)

dev.off()


# Formal test of 'parallel lines' in PCA plot
pdf("plots/cna_TITAN_FACETS_PC2_scatter.pdf", width=4, height=4.5)
scatterPlot(pca$x[titan_idx, 2], pca$x[facets_idx, 2],
	xlab="PC2 (TITAN)",
	ylab="PC2 (FACETS)"
)
dev.off()


# CNA-mRNA enrichment comparison
# ---------------------------

enrich_cna_mrna = rbind(
	read.csv("data/GSEA/CNA_mRNA/GseaPreranked_hallmarks_all_m10000/gsea_report_for_na_pos_1555074277249.csv"),
	read.csv("data/GSEA/CNA_mRNA/GseaPreranked_hallmarks_all_m10000/gsea_report_for_na_neg_1555074277249.csv")
)


enrich_cna_mrna_facets = rbind(
	read.csv("data/GSEA/CNA_mRNA_FACETS/my_analysis.GseaPreranked.1562337229186/CNA-mRNA_FACETS_gsea_report_for_na_pos_1562337229186.csv"),
	read.csv("data/GSEA/CNA_mRNA_FACETS/my_analysis.GseaPreranked.1562337229186/CNA-mRNA_FACETS_gsea_report_for_na_neg_1562337229186.csv")
)

# Match by term
enrich_cna_mrna_facets = enrich_cna_mrna_facets[match(enrich_cna_mrna$NAME, enrich_cna_mrna_facets$NAME), ]

stopifnot(enrich_cna_mrna$NAME == enrich_cna_mrna_facets$NAME)


pdf("plots/scatter_GSEA_enrich_TITAN_FACETS.pdf", width=5.0, height=5)
scatterPlot(enrich_cna_mrna$NES, enrich_cna_mrna_facets$NES,
	xlab="GSEA NES (CNA-mRNA correlation, TITAN)",
	ylab="GSEA NES (CNA-mRNA correlation, FACETS)",
	method="kendall"
)
dev.off()

