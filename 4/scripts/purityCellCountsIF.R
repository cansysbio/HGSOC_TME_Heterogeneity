# Compares tumor purity estimated by TitanCNA with quantification of infiltrating leukocytes 
# from immuno flourescence microscopy data
#

rm(list=ls())

library(RColorBrewer)
library(stringr)
library(data.table)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

# Load TitanCNA results table
# ---------------------------------------
opt_clust = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")

colors = loadPatientColors(opt_clust)


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
opt_clust = merge(opt_clust, wes_label[, c("barcode", "case", "site", "color", "sample_id")],
	by="barcode", all.x=TRUE)


# TitanCNA purity or ploidy associtions with TIL (tumor infiltrating leukocytes)
# Based on imaging data.
# -------------------------------------------------------------

# selection = c("cd4_frac", "cd4foxp3_frac", "foxp3_frac", "cd8_frac", "cd4cd8_comb_frac", "sample_id")
# selection = c("cd4_frac", "cd4foxp3_frac", "foxp3_frac", "cd8_frac", "cd4cd8_comb_frac", "T_help", "sample_id")
selection = c("cd4_frac", "cd4foxp3_frac", "foxp3_frac", "cd8_frac", "cd4cd8_comb_frac", "cd4cd8foxp3_comb_frac", "T_help", "sample_id")

mean_cell_frac = aggregate(
	.~sample_id,
	data=cell_frac[, selection],
	mean)

sd_cell_frac = aggregate(
	.~sample_id,
	data=cell_frac[, selection],
	sd)

n_cell_frac = aggregate(
	.~sample_id,
	data=cell_frac[, selection],
	length)

# 95% confidence interval
ci_cell_frac = 1.96 * sd_cell_frac[, -1] / sqrt(n_cell_frac[, -1])

# Match the purity tumors to IF data
idx = match(as.character(opt_clust$sample_id), as.character(mean_cell_frac$sample_id))

sum(!is.na(idx), na.rm=TRUE)
cat(unique(as.character(opt_clust$tumor[is.na(idx)])), sep="\n")

mean_cell_frac_matched = mean_cell_frac[idx, ]
ci_cell_frac_matched = ci_cell_frac[idx, ]


pdf("plots/TIL_purity_ploidy_cor_v6.pdf", width=10.5, height=2.7)
par(mfcol=c(1, 6))

feat = "cd8_frac"
scatterPlotCI(mean_cell_frac_matched[[feat]] * 100, opt_clust$purity,
	ci_cell_frac_matched[[feat]] * 100,
	ylab="Tumour purity (TitanCNA)", xlab="CD8+ (%)",
	bg=brewer.pal(9, "Set1")[1]
)

feat = "cd4_frac"
scatterPlotCI(mean_cell_frac_matched[[feat]] * 100, opt_clust$purity,
	ci_cell_frac_matched[[feat]] * 100,
	ylab="Tumour Purity (TitanCNA)", xlab="CD4+ (%)",
	bg=brewer.pal(9, "Set1")[2]
)

feat = "cd4cd8_comb_frac"
scatterPlotCI(mean_cell_frac_matched[[feat]] * 100, opt_clust$purity,
	ci_cell_frac_matched[[feat]] * 100,
	ylab="Tumour Purity (TitanCNA)", xlab="CD4+ or CD8+ (%)",
	bg=brewer.pal(9, "Set1")[3]
)

feat = "foxp3_frac"
scatterPlotCI(mean_cell_frac_matched[[feat]] * 100, opt_clust$purity,
	ci_cell_frac_matched[[feat]] * 100,
	ylab="Tumour Purity (TitanCNA)", xlab="Foxp3+ (%)",
	bg=brewer.pal(9, "Set1")[4]
)

feat = "cd4foxp3_frac"
scatterPlotCI(mean_cell_frac_matched[[feat]] * 100, opt_clust$purity,
	ci_cell_frac_matched[[feat]] * 100,
	ylab="Tumour Purity (TitanCNA)", xlab="CD4+Foxp3+ (%)",
	bg=brewer.pal(9, "Set1")[5]
)

feat = "cd4cd8foxp3_comb_frac"
scatterPlotCI(mean_cell_frac_matched[[feat]] * 100, opt_clust$purity,
	ci_cell_frac_matched[[feat]] * 100,
	ylab="Tumour Purity (TitanCNA)", xlab="CD4+ or CD8+ or Foxp3+ (%)",
	bg=brewer.pal(9, "Set1")[7]
)

dev.off()
