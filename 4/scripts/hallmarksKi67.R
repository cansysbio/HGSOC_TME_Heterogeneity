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


# Load CNA count matrix
cna = loadCNA("data/cna.csv")

cna_mat = data.matrix(cna[, c(-1, -2)])
rownames(cna_mat) = cna$hgnc_symbol

colnames(cna_mat)[colnames(cna_mat) == "RG13T12"] = "RG13T122"  # Reintroduce typo for consistency with mRNA annotation


# Load sample annotation, for matching
expr_samples = loadSampleAnnot(file_path="data/TreatmentNaive_SampleLabels_WESTumourCellularity_mRNAtumourCellularity_MAPPINGS.csv")

cna_mat_matched = cna_mat[,
	match(
		expr_samples$wes_label[match(colnames(hallmarks), expr_samples$Well)],
		colnames(cna_mat)
	)
]


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



# Scatterplot with indicator for amplification of gene
pathway = "HALLMARK_G2M_CHECKPOINT"
gene = "KRAS"
amp_idx = cna_mat_matched[which(rownames(cna_mat_matched) == gene), ] > 5

nes = hallmarks[which(rownames(hallmarks) == pathway), ]

pt_col = rep("white", ncol(hallmarks))
pt_col[amp_idx] = brewer.pal(9, "Set1")[1]

pdf(paste0("plots/ssGSEA_", pathway, "_Ki67_amp", gene, ".pdf"), width=4.0, height=4.5)
scatterPlot(if_frac_matched$Ki67, nes,
	xlab="Ki67 (% positive cells)",
	ylab=paste0(pathway, " (NES, mRNA ssGSEA)"),
	method="kendall",
	# bg=brewer.pal(9, "Pastel1")[1],
	bg=pt_col,
	pch=21)
dev.off()


# idx = nes > 0.45
# sum(idx)

# scatterPlot(if_frac_matched$Ki67[!idx], nes[!idx], method="spearman")
# if_frac_matched[idx, ]
# scatterPlot(if_frac_matched$Ki67, nes, method="kendall")
