# Copy-number mRNA expression boxplots
#
rm(list=ls())

library(RColorBrewer)
library(data.table)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")

# Load CNA count matrix
cna = fread("data/cna.csv")

cna_mat = data.matrix(cna[, c(-1, -2)])
rownames(cna_mat) = cna$hgnc_symbol

colnames(cna_mat)[colnames(cna_mat) == "RG13T12"] = "RG13T122"  # Reintroduce typo for consistency with mRNA annotation

# Load gene expression matrix
expr = fread("data/OVCTp_log2exp_loess_norm.txt")

# expr_samples = loadSampleAnnot("data/OVCT_Tnaive_WES_labels.txt")
expr_samples = loadSampleAnnot("data/TreatmentNaive_SampleLabels_WESTumourCellularity_mRNAtumourCellularity_MAPPINGS.csv")


# expression matrix
expr_mat = data.matrix(expr[, -1])

colnames(expr_mat) = expr_samples$wes_label[match(colnames(expr_mat), expr_samples$Well)]
rownames(expr_mat) = expr$Hugo_Symbol

# Match expression matrix to cna_mat
expr_mat_match = expr_mat[
	match(rownames(cna_mat), rownames(expr_mat)),
	match(colnames(cna_mat), colnames(expr_mat))
]

message("Matching samples: ", sum(!is.na(match(colnames(cna_mat), colnames(expr_mat)))))

# Include genes where we have at least one measurement from microarray data
include_genes = apply(expr_mat_match, 1, function(row) {
	any(!is.na(row))
})

cna_mat_filter = cna_mat[include_genes, ]
expr_mat_match_filter = expr_mat_match[include_genes, ]


pdf("plots/cna_mRNA_boxplots.pdf", width=5, height=6)
# gene_list = c("RAD51AP1", "PDCD10", "ETV7", "CEBPG")
# par(mfrow=c(2, 2))


pdf("plots/cna_mRNA_boxplots.pdf", width=5, height=6)
gene_list = c("MYC", "MYCN")
par(mfrow=c(1, 2))
for (gene in gene_list) {
	i = which(rownames(cna_mat_filter) == gene)

	stopifnot(rownames(cna_mat_filter)[i] == rownames(expr_mat_match_filter)[i])

	# Group expression by copy number
	cna_range = 0:8
	expr_by_cna = lapply(cna_range, function(n) {
		idx = cna_mat_filter[i, ] == n
		values = expr_mat_match_filter[i, idx]
		values = na.omit(values)
		return(values)
	})
	names(expr_by_cna) = cna_range

	boxplot(expr_by_cna,
		col=brewer.pal(9, "Greens"),
		ylab=paste(gene, "mRNA"), xlab="Copy number",
		frame=FALSE
	)

	for (i in 1:length(expr_by_cna)) {
		values = expr_by_cna[[i]]
		points(rep(i, length(values)), values,
			bg=brewer.pal(9, "Set1")[3],
			pch=21)
	}
}
dev.off()
