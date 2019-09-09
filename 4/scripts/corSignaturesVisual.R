rm(list=ls())

library(data.table)
library(RColorBrewer)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

# Load correlation tables
cna_mrna_titan = fread("data/corSig/CNA_mRNA_cor_spearman.csv")
mrna_purity_titan = fread("data/corSig/mRNA_purity_cor_spearman.csv")
cna_mrna_facets = fread("data/corSig/CNA_mRNA_cor_spearman_FACETS.csv")


# load hallmark gene sets from msigDB
hallmark = GSA.read.gmt("data/msigDB/h.all.v6.2.symbols.gmt")

names(hallmark$genesets) = hallmark$geneset.names


# Plot correlation coefficient histogram
pdf("plots/cna_mRNA_histogram_TITAN_ASCAT.pdf")
par(mfrow=c(2, 1))
hist(cna_mrna_titan$r, xlim=c(-1, 1),
	main="TITAN",
	xlab="CNA-mRNA correlation",
	col=brewer.pal(9, "Pastel1")[1],
	breaks=seq(-1, 1, length.out=80)
)
hist(cna_mrna_facets$r, xlim=c(-1, 1),
	main="FACETS",
	xlab="CNA-mRNA correlation",
	col=brewer.pal(9, "Pastel1")[2],
	breaks=seq(-1, 1, length.out=80)
)
dev.off()


# Myc targets tables
# -----------------------
myc_genes = unique(c(
	hallmark$genesets$HALLMARK_MYC_TARGETS_V2,
	hallmark$genesets$HALLMARK_MYC_TARGETS_V1
))

write.csv(
	cna_mrna_titan[cna_mrna_titan$gene_symbol %in% myc_genes, c(-2, -5)],
	file="data/corSig/MycTargets/cna_mRNA_cor_spearman_myc.csv",
	row.names=FALSE
)

write.csv(
	mrna_purity_titan[mrna_purity_titan$gene_symbol %in% myc_genes, ],
	file="data/corSig/MycTargets/mRNA_purity_cor_spearman_myc.csv",
	row.names=FALSE
)

# Wnt beta catening signaling tables
# ---------------------------
wnt_genes = hallmark$genesets$HALLMARK_WNT_BETA_CATENIN_SIGNALING

write.csv(
	cna_mrna_titan[cna_mrna_titan$gene_symbol %in% wnt_genes, c(-2, -5)],
	file="data/corSig/WntBetaCat/CNA_mRNA_cor_spearman_wnt.csv",
	row.names=FALSE
)

write.csv(
	mrna_purity_titan[mrna_purity_titan$gene_symbol %in% wnt_genes, ],
	file="data/corSig/WntBetaCat/mRNA_purity_cor_spearman_wnt.csv",
	row.names=FALSE
)


# Other genes
# ---------------------------------
genes = c("IL23A", "CCL23", "CXCL10", "CXCL9", "CCL4")

cna_mrna_titan[cna_mrna_titan$gene_symbol %in% genes, -2]
mrna_purity_titan[mrna_purity_titan$gene_symbol %in% genes, ]
