
rm(list=ls())

library(data.table)


setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")


# Load ssGSEA hallmark enrichment results
hallmarks = fread("data/HallmarksStromaImmune_NES.txt")
terms = hallmarks$Term
hallmarks = data.matrix(hallmarks[, -1])
rownames(hallmarks) = terms

# Load CNA count matrix
cna = loadCNA("data/cna.csv")

cna_mat = data.matrix(cna[, c(-1, -2)])
rownames(cna_mat) = cna$hgnc_symbol

colnames(cna_mat)[colnames(cna_mat) == "RG13T12"] = "RG13T122"  # Reintroduce typo for consistency with mRNA annotation


# Sample annotation
expr_samples = loadSampleAnnot(file_path="data/TreatmentNaive_SampleLabels_WESTumourCellularity_mRNAtumourCellularity_MAPPINGS.csv")


# Match hallmark labels to WES IDs
colnames(hallmarks) = expr_samples$wes_label[match(colnames(hallmarks), expr_samples$Well)]

hallmarks_matched = hallmarks[,
	match(colnames(cna_mat), colnames(hallmarks))
]


# Select association to analyze
# --------------------------------------
rownames(hallmarks_matched)
gene = "KRAS"
pathway = "HALLMARK_G2M_CHECKPOINT"

# scatterPlot(
# 	cna_mat[which(rownames(cna_mat) == gene), ],
# 	hallmarks_matched[which(rownames(hallmarks_matched) == pathway), ],
# 	xlab=gene,
# 	ylab=paste0(pathway, " (NES)"),
# 	method="spearman"
# 	# method="pearson"
# )


# Collect NES statistics for amplified vs non-amplified copy number samples
nes = hallmarks_matched[which(rownames(hallmarks_matched) == pathway), ]
amp_idx = cna_mat[which(rownames(cna_mat) == gene), ] > 5

amp_nes = list(
	"-"=nes[!amp_idx],
	"+"=nes[amp_idx]
)

t_test = t.test(amp_nes[[1]], amp_nes[[2]])


pdf(paste0("plots/ssGSEA_copyNumber_", gene, "_", pathway, ".pdf"), width=2.3, height=5)
# Boxplot
bp = boxplot(amp_nes,
	frame=FALSE,
	ylab=paste0(pathway, " (NES, mRNA ssGSEA)"),
	xlab=paste0(gene, " amp. (>5 CN)"),
	col=rgb(240, 240, 240, maxColorValue=255),
	main=paste0("P=", format(t_test$p.value, digits=3))
)

# Jitter points
colors = c("white", brewer.pal(9, "Set1")[1])
for (i in 1:length(amp_nes)) {
	points(
		jitter(rep(i, length(amp_nes[[i]])), amount=0.2),
		amp_nes[[i]],
		pch=21,
		bg=colors[i]
	)
}
dev.off()
