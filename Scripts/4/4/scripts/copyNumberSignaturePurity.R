# Analysis of association between TitanCNA purity and copy-number signatures.

rm(list=ls())

library(data.table)
library(RColorBrewer)
library(gmodels)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

titan = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")

patient_colors = loadPatientColors(titan)

# cn_sig = readRDS("data/signatureExposure/miller_cnsignatures.rds")  # TITAN signatures
cn_sig = readRDS("data/signatureExposure/copywriteR_ovct_sigs_20191126.rds")  # copywriteR signatures


# Fitler all
# ---------------------
# exclude_samples = as.character(read.table("data/QC/QCexclude.txt")[, 1])

# Filter signature data
# cn_sig = cn_sig[, !colnames(cn_sig) %in% exclude_samples]
# titan = titan[!titan$barcode %in% exclude_samples, ]
# --------------------------------------------

# Fraction of samples with more than 20% exposure by signature
apply(cn_sig > 0.2, 1, mean)

apply(cn_sig > 0.05, 1, mean)

# Mean and confidence intervals 
apply(cn_sig, 1, ci)

# Fraction of samples with particular highest (maximum) signature
prop.table(table(apply(cn_sig, 2, which.max)))

cn_sig[, grep("RG4", colnames(cn_sig))]

cn_sig[, grep("RG6", colnames(cn_sig))]

# Test if tumor labels are aligned
stopifnot(all(titan$barcode == colnames(cn_sig)))


pdf("plots/cn_signature_scatter.pdf", height=2.7, width=5.5)
pt_col = patient_colors[as.integer(titan$patient_id)]

par(mfrow=c(1, 3))
for (i in c(1, 6, 4)) {
	scatterPlot(cn_sig[i, ], titan$purity,
		xlab=rownames(cn_sig)[i],
		ylab="Tumour purity (TitanCNA)",
		pch=21,
		bg=pt_col
	)
}
dev.off()

# Stacked barplot ordered by purity
# pdf("plots/cn_signature_purity_order.pdf", height=3.5, width=9)
pdf("plots/cn_signature_purity_order_copywriteR.pdf", height=3.5, width=9)

# pdf("plots/cn_signature_purity_order_copywriteR_QCfilter.pdf", height=3.5, width=9)
idx = order(titan$purity, decreasing=TRUE)

sig_colors = brewer.pal(8, "Set2")

bp = barplot(cn_sig[, idx],
	col=sig_colors,
	las=2,
	legend.text=row.names(cn_sig),
	ylab="Exposure, purity",
	xlab="Tumours",
	cex.names=0.8
)

lines(bp, titan$purity[idx], col=brewer.pal(9, "Set1")[1], lwd=2)

points(bp, rep(-0.05, length(bp)),
	pch=21,
	bg=pt_col[idx],
	xpd=TRUE)
dev.off()



# Combination of signatures
# -----------------------------

boxplotTest = function(values, ...) {
	t_test1 = t.test(values[[1]], values[[2]]) 
	t_test2 = t.test(values[[3]], values[[4]]) 
	t_test3 = t.test(values[[5]], values[[6]]) 

	boxplot(values,
		las=2,
		main=paste0(
			"P=", format(t_test1$p.value, digits=3),
			", P=", format(t_test2$p.value, digits=3),
			", P=", format(t_test3$p.value, digits=3)
			),
		frame=FALSE,
		...
	)

	colors = c("white",
		brewer.pal(9, "Set1")[1],
		"white",
		brewer.pal(9, "Set1")[2],
		"white",
		brewer.pal(9, "Set1")[3]
	)

	for (k in 1:length(values)) {
		points(
			jitter(rep(k, length(values[[k]])), amount=0.2),
			values[[k]],
			bg=colors[k],
			pch=21
		)
	}
}


cn_above_median = apply(cn_sig, 1, function(exposure) {
	return(exposure > median(exposure))
})


i = 4
j = 7
pdf("plots/cn_signature_purity_comb_boxplot.pdf", height=4.0, width=4)
idx_comb = cn_above_median[, i] & cn_above_median[, j]  # both high
idx_ctrl = !cn_above_median[, i] & !cn_above_median[, j]  # both low

idx_i = cn_above_median[, i]  # i high
idx_j = cn_above_median[, j]  # j high

values = list(
	"Low s4"=titan$purity[!idx_i],
	"High s4"=titan$purity[idx_i],
	"Low s7"=titan$purity[!idx_j],
	"High s7"=titan$purity[idx_j],
	"Low s4 and s7"=titan$purity[idx_ctrl],
	"High s4 and s7"=titan$purity[idx_comb]
)

boxplotTest(values, ylab="Tumor cellularity")

dev.off()



# Load IF cell fractions for tumor slices
# -------------------------------------------------
cell_frac = read.table("data/IF_quantification_ajs.txt", sep="\t", header=TRUE)

# Calculate fractions
cell_frac$cd4_frac = cell_frac$CD4 / cell_frac$cell_counts
cell_frac$cd8_frac = cell_frac$CD8 / cell_frac$cell_counts
cell_frac$foxp3_frac = cell_frac$FoxP3 / cell_frac$cell_counts
cell_frac$cd4foxp3_frac = cell_frac$CD4FoxP3 / cell_frac$cell_counts

cell_frac$sample_id = paste(cell_frac$case, cell_frac$site, cell_frac$color, sep="_")

# Combined fractions
cell_frac$cd4cd8_comb_frac = cell_frac$cd4_frac + cell_frac$cd8_frac  # assumes that there are no CD4+CD8+ double positive cells
cell_frac$cd4cd8foxp3_comb_frac = cell_frac$cd4_frac + cell_frac$cd8_frac + cell_frac$foxp3_frac - cell_frac$cd4foxp3_frac  # corrects for double counting

cell_frac$T_help = cell_frac$cd4_frac - cell_frac$cd4foxp3_frac


# Load WES site labels, specifying which site and color (area) that the WES samples were extracted from
# --------------------------------
wes_label = read.csv("data/Radiogenomics IDs.csv")

# Remove underscores from tumor IDs, for compatability with labels from WES files.
# barcode is ID column from TitanCNA results table
wes_label$barcode = gsub("_", "", wes_label$Sample.ID)

wes_label$case = str_match(wes_label$barcode, "^RG([0-9]+)")[, 2]

# Clean up labels
wes_label$site = "other"
wes_label$site[grep("ovary", tolower(wes_label$Paper.ID))] = "ovary"
wes_label$site[grep("omentum", tolower(wes_label$Paper.ID))] = "omentum"

wes_label$color = NA
wes_label$color[grep("blue", tolower(wes_label$Paper.ID))] = "blue"
wes_label$color[grep("yellow", tolower(wes_label$Paper.ID))] = "yellow"
wes_label$color[grep("green", tolower(wes_label$Paper.ID))] = "green"

# Define how samples are matched between IF and WES data, taking into account the site and color
wes_label$sample_id = paste(wes_label$case, wes_label$site, wes_label$color, sep="_")


# Add sample IDs to TitanCNA output table
# -----------------------------------------------
titan = merge(titan, wes_label[, c("barcode", "case", "site", "color", "sample_id")],
	by="barcode", all.x=TRUE)


selection = c("cd4_frac", "cd4foxp3_frac", "foxp3_frac", "cd8_frac", "cd4cd8_comb_frac", "cd4cd8foxp3_comb_frac", "T_help", "sample_id")

mean_cell_frac = aggregate(
	.~sample_id,
	data=cell_frac[, selection],
	mean)


# Match the purity tumors to IF data
idx = match(as.character(titan$sample_id), as.character(mean_cell_frac$sample_id))

sum(!is.na(idx), na.rm=TRUE)
cat(unique(as.character(titan$tumor[is.na(idx)])), sep="\n")

mean_cell_frac_matched = mean_cell_frac[idx, ]



# Boxplot of low high signatures vs IF CD4 and CD8 cell fractions

i = 4
j = 7

idx_comb = cn_above_median[, i] & cn_above_median[, j]  # both high
idx_ctrl = !cn_above_median[, i] & !cn_above_median[, j]  # both low

idx_i = cn_above_median[, i]  # i high
idx_j = cn_above_median[, j]  # j high

pdf("plots/cn_signature_CD8_comb_boxplot.pdf", height=4.0, width=4)

values = list(
	"Low s4"=mean_cell_frac_matched$cd8_frac[!idx_i],
	"High s4"=mean_cell_frac_matched$cd8_frac[idx_i],
	"Low s7"=mean_cell_frac_matched$cd8_frac[!idx_j],
	"High s7"=mean_cell_frac_matched$cd8_frac[idx_j],
	"Low s4 and s7"=mean_cell_frac_matched$cd8_frac[idx_ctrl],
	"High s4 and s7"=mean_cell_frac_matched$cd8_frac[idx_comb]
)
boxplotTest(values, ylab="CD8 fraction")

dev.off()


pdf("plots/cn_signature_CD4_comb_boxplot.pdf", height=4.0, width=4)

values = list(
	"Low s4"=mean_cell_frac_matched$cd4_frac[!idx_i],
	"High s4"=mean_cell_frac_matched$cd4_frac[idx_i],
	"Low s7"=mean_cell_frac_matched$cd4_frac[!idx_j],
	"High s7"=mean_cell_frac_matched$cd4_frac[idx_j],
	"Low s4 and s7"=mean_cell_frac_matched$cd4_frac[idx_ctrl],
	"High s4 and s7"=mean_cell_frac_matched$cd4_frac[idx_comb]
)
boxplotTest(values, ylab="CD4 fraction")

dev.off()


pdf("plots/cn_signature_cd4cd8foxp3_comb_boxplot.pdf", height=4.0, width=4)

values = list(
	"Low s4"=mean_cell_frac_matched$cd4cd8foxp3_comb_frac[!idx_i],
	"High s4"=mean_cell_frac_matched$cd4cd8foxp3_comb_frac[idx_i],
	"Low s7"=mean_cell_frac_matched$cd4cd8foxp3_comb_frac[!idx_j],
	"High s7"=mean_cell_frac_matched$cd4cd8foxp3_comb_frac[idx_j],
	"Low s4 and s7"=mean_cell_frac_matched$cd4cd8foxp3_comb_frac[idx_ctrl],
	"High s4 and s7"=mean_cell_frac_matched$cd4cd8foxp3_comb_frac[idx_comb]
)
boxplotTest(values, ylab="CD4 or CD8 or Foxp3 fraction")

dev.off()



# Load mRNA purity estimates from Alejandro
# sample annotation
# sample_annot = loadSampleAnnot("data/OVCT_Tnaive_WES_labels.txt")
sample_annot = loadSampleAnnot("data/TreatmentNaive_SampleLabels_WESTumourCellularity_mRNAtumourCellularity_MAPPINGS.csv")
sample_annot$exp_label = sample_annot$Well

conv_scores = loadConvScores("data/OVCTp_log2exp_loess_norm_estimate_score.txt")

# combine purity scores with sample annotation data
conv = merge(conv_scores, sample_annot, by="exp_label")

# Match conv to titan results
# conv_match = conv[match(titan$barcode, conv$tumor), ]
message("Samples matched: ", sum(!is.na(match(as.character(titan$barcode), conv$wes_label))))

conv_match = conv[match(as.character(titan$barcode), conv$wes_label), ]


pdf("plots/cn_signature_immuneScore_comb_boxplot.pdf", height=4.0, width=4)

values = list(
	"Low s4"=conv_match$ImmuneScore[!idx_i],
	"High s4"=conv_match$ImmuneScore[idx_i],
	"Low s7"=conv_match$ImmuneScore[!idx_j],
	"High s7"=conv_match$ImmuneScore[idx_j],
	"Low s4 and s7"=conv_match$ImmuneScore[idx_ctrl],
	"High s4 and s7"=conv_match$ImmuneScore[idx_comb]
)
boxplotTest(values, ylab="ImmuneScore")

dev.off()


pdf("plots/cn_signature_ESTIMATEScore_comb_boxplot.pdf", height=4.0, width=4)

values = list(
	"Low s4"=conv_match$ESTIMATEScore[!idx_i],
	"High s4"=conv_match$ESTIMATEScore[idx_i],
	"Low s7"=conv_match$ESTIMATEScore[!idx_j],
	"High s7"=conv_match$ESTIMATEScore[idx_j],
	"Low s4 and s7"=conv_match$ESTIMATEScore[idx_ctrl],
	"High s4 and s7"=conv_match$ESTIMATEScore[idx_comb]
)
boxplotTest(values, ylab="ESTIMATEScore")

dev.off()



# Sample analysis as above filtered for QC
# ----------------------------------------
stopifnot(all(titan$barcode == colnames(cn_sig)))

# included_samples = readRDS("data/QC/ovct_filtered_samples.rds")  # From Geoff MacIntyre
# idx_filter = colnames(cn_sig) %in% included_samples


exclude_samples = as.character(read.table("data/QC/QCexclude.txt")[, 1])
idx_filter = !colnames(cn_sig) %in% exclude_samples

cn_above_median_filter = apply(cn_sig[, idx_filter], 1, function(exposure) {
	return(exposure > median(exposure))
})


i = 4
j = 7

idx_comb = cn_above_median_filter[, i] & cn_above_median_filter[, j]  # both high
idx_ctrl = !cn_above_median_filter[, i] & !cn_above_median_filter[, j]  # both low

idx_i = cn_above_median_filter[, i]  # i high
idx_j = cn_above_median_filter[, j]  # j high


pdf("plots/cn_signature_immuneScore_comb_boxplot_QCfilter.pdf", height=4.0, width=4)
values = list(
	"Low s4"=conv_match$ImmuneScore[idx_filter][!idx_i],
	"High s4"=conv_match$ImmuneScore[idx_filter][idx_i],
	"Low s7"=conv_match$ImmuneScore[idx_filter][!idx_j],
	"High s7"=conv_match$ImmuneScore[idx_filter][idx_j],
	"Low s4 and s7"=conv_match$ImmuneScore[idx_filter][idx_ctrl],
	"High s4 and s7"=conv_match$ImmuneScore[idx_filter][idx_comb]
)
boxplotTest(values, ylab="ImmuneScore")
dev.off()


pdf("plots/cn_signature_ESTIMATEScore_comb_boxplot_QCfilter.pdf", height=4.0, width=4)

values = list(
	"Low s4"=conv_match$ESTIMATEScore[idx_filter][!idx_i],
	"High s4"=conv_match$ESTIMATEScore[idx_filter][idx_i],
	"Low s7"=conv_match$ESTIMATEScore[idx_filter][!idx_j],
	"High s7"=conv_match$ESTIMATEScore[idx_filter][idx_j],
	"Low s4 and s7"=conv_match$ESTIMATEScore[idx_filter][idx_ctrl],
	"High s4 and s7"=conv_match$ESTIMATEScore[idx_filter][idx_comb]
)
boxplotTest(values, ylab="ESTIMATEScore")

dev.off()



# G2M and KRAS amp associations, filtered signatures
# -------------------------------------------

boxplotTest1 = function(values, colors=NA, ...) {
	t_test1 = t.test(values[[1]], values[[2]]) 

	boxplot(values,
		las=2,
		main=paste0(
			"P=", format(t_test1$p.value, digits=3)
			),
		frame=FALSE,
		...
	)

	if (is.na(colors)) {
		colors = list(
			"white",
			brewer.pal(9, "Set1")[1]
		)
	}

	for (k in 1:length(values)) {
		points(
			jitter(rep(k, length(values[[k]])), amount=0.2),
			values[[k]],
			bg=colors[[k]],
			pch=21
		)
	}
}


# Load ssGSEA hallmark enrichment results
# -----------------------------------
hallmarks = fread("data/HallmarksStromaImmune_NES.txt")
terms = hallmarks$Term
hallmarks = data.matrix(hallmarks[, -1])
rownames(hallmarks) = terms


# Sample annotation
expr_samples = loadSampleAnnot(file_path="data/TreatmentNaive_SampleLabels_WESTumourCellularity_mRNAtumourCellularity_MAPPINGS.csv")

# Match hallmark labels to WES IDs
colnames(hallmarks) = expr_samples$wes_label[match(colnames(hallmarks), expr_samples$Well)]

hallmarks_match = hallmarks[, match(rownames(cn_above_median_filter), colnames(hallmarks))]



# Load CNA count matrix
# --------------------
cna = loadCNA("data/cna.csv")

cna_mat = data.matrix(cna[, c(-1, -2)])
rownames(cna_mat) = cna$hgnc_symbol

colnames(cna_mat)[colnames(cna_mat) == "RG13T12"] = "RG13T122"  # Reintroduce typo for consistency with mRNA annotation

cna_mat_match = cna_mat[, match(rownames(cn_above_median_filter), colnames(cna_mat))]
cna_mat_match = data.matrix(cna_mat_match)


pdf("plots/s1_associations_v2.pdf", width=2.4)
par(mfrow=c(2, 1))
# s1 vs G2M enrichment
# ---------------------
pathway = "HALLMARK_G2M_CHECKPOINT"
idx = cn_above_median_filter[, 1]  # s1


# Colors based on KRAS amplification
pts_col = list(
	rep("white", sum(!idx)),
	rep("white", sum(idx))
)

i = which(rownames(cna_mat_match) == gene)
pts_col[[1]][cna_mat_match[i, !idx] >= 6] = brewer.pal(9, "Set1")[1]
pts_col[[2]][cna_mat_match[i, idx] >= 6] = brewer.pal(9, "Set1")[1]


# Test if labels are matched
stopifnot(all(rownames(cn_above_median_filter) == colnames(hallmarks_match), na.rm=TRUE))

values = list(
	"Low s1"=hallmarks_match[rownames(hallmarks_match) == pathway, !idx],
	"High s1"=hallmarks_match[rownames(hallmarks_match) == pathway, idx]
)

boxplotTest1(values, colors=pts_col, ylab="G2M erichment")


# S1 vs KRAS amp
# ----------------------
gene = "KRAS"
idx = cn_above_median_filter[, 1]  # s1

stopifnot(all(rownames(cn_above_median_filter) == colnames(cna_mat_match), na.rm=TRUE))

i = which(rownames(cna_mat_match) == gene)
values = list(
	"Low s1"=cna_mat_match[i, !idx],
	"High s1"=cna_mat_match[i, idx]
)

boxplotTest1(values, ylab=paste(gene, " copy number"))
dev.off()