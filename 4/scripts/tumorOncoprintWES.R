# Plots oncoprint representation of whole-exome sequencing data. Both mutation and copy-number calls.
rm(list=ls())

library(ComplexHeatmap)
library(data.table)
library(RColorBrewer)
library(stringr)
library(circlize)

library(compiler)
enableJIT(3)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")


# Load TitanCNA results table
# ----------------------------------------
opt_clust = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")

patient_colors = loadPatientColors(opt_clust)

names(patient_colors) = levels(opt_clust$patient_id)  # named colors used as input to oncoPrint

opt_clust$patient = str_extract(opt_clust$patient_id, "[0-9]+")

# Load CNA count matrix from whole-exome sequencing data
# -------------------------------------
cna = fread("data/cna.csv")

cna_mat = data.matrix(cna[, c(-1, -2)])
rownames(cna_mat) = cna$hgnc_symbol

tumor_ids = colnames(cna_mat)  # global specification of tumor ids
patient_ids = str_extract(tumor_ids, "^RG([0-9]+)")  # patient ids

# Load SNP calls from whole-exome sequencing data
# Columns of interest: Variant_Type, Variant_Classification, Tumor_Sample_Barcode, Hugo_Symbol
# ----------------------------------------
mut = fread("data/oncotator_maf_ovct_filtered_updated.maf")  # using the padded exome bed file (@Danish)

# Mutation classes prepended by SNP/DEL/INS...
mut$Variant_Code = paste(mut$Variant_Type, mut$Variant_Classification, sep="_")


# Load clinical data
# --------------------------------
clinical = read.csv("data/clinical/TableS1.csv")  # contains BRCA2 germline status


# oncoprint from ovarian tumor samples
# ------------------------------------------

# Ext list 2
genes = c("BRCA1", "BRCA2", "NF1", "PTEN", "RB1", "TP53",
    "KRAS", "MYC", "SRC",
    "ERBB2",  # HER2
    "BRAF"
)

# Ext list 2, including all Myc family genes
genes = c("BRCA1", "BRCA2", "NF1", "PTEN", "RB1", "TP53",
    "KRAS",
    "MYC",
    "MYCL",
    "MYCN",
    "SRC",
    "ERBB2",  # HER2
    "BRAF"
)



mut_type_all = getVariantClassIndicatorMatrices(genes, tumor_ids, mut)

# Manually added mutation annotations. Based on manual inspection of DNA aligments.
mut_type_all[["Frame_Shift_Del"]][genes == "TP53", grep("RG4T", tumor_ids)] = 1

# Select and group mutation indicator matrices
mut_type = c(
	mut_type_all[c("Missense_Mutation", "Nonsense_Mutation", "Intron", "5Flank", "3UTR")],
	# composite indicators, assumes that <1 is ok
	list(
		Indel_Frameshift=mut_type_all$Frame_Shift_Del + mut_type_all$Frame_Shift_Ins,
		Indel_Inframe=mut_type_all$In_Frame_Del + mut_type_all$In_Frame_Ins
	)
)


# Known germline mutation status
mut_type$germline = matrix(0, nrow=length(genes), ncol=nrow(opt_clust))
rownames(mut_type$germline) = genes
colnames(mut_type$germline) = opt_clust$barcode

brca2_patients = clinical$CASE_ID[clinical$BRCA_STATUS == "BRCA2"]
mut_type$germline[genes == "BRCA2", opt_clust$patient %in% brca2_patients] = 1

# Parse copy-numbers into list of indicator matrices
# ---------------------------------
cna_mat_sub = cna_mat[match(genes, rownames(cna_mat)), ]


# Generate list of annotation matrices
cna_mat_list = list(
	deletion0=(cna_mat_sub == 0) * 1,
	deletion1=(cna_mat_sub == 1) * 1,
	amplification34=(cna_mat_sub > 2 & cna_mat_sub <=4) * 1,
	amplification5=(cna_mat_sub > 4) * 1
)

# set missing values to zero
cna_mat_list = lapply(cna_mat_list, function(mat) {
	mat[is.na(mat)] = 0
	return(mat)
})


# Oncoprint based on lists of indicator matrices for the gene set
# Dependencies: cna_mat_list, mut_type
# ------------------------------------

# pdf("plots/oncoprint_extended_list2_v4_manualTP53.pdf", height=3.5, width=10)  # with deletion1
pdf("plots/oncoprint_extended_list2_v4_manualTP53_MYC.pdf", height=3.5, width=10)  # with deletion1

colors=c("black", brewer.pal(8, "Accent"))

box_col = c(
    deletion0=brewer.pal(9, "Set1")[2],
    deletion1=brewer.pal(9, "Blues")[4],
    # amplification34=brewer.pal(9, "Reds")[4],
    amplification5=brewer.pal(9, "Set1")[1],
    Missense_Mutation=colors[1],
    Nonsense_Mutation=colors[2],
    Indel_Frameshift=colors[3],
    Indel_Inframe=colors[4],
    Intron=colors[5],
    '5Flank'=colors[6],
    '3UTR'=colors[7],
    germline="grey"
)

outline_col = "black"

# Tests if rownames and column names are consistent across indicator matrices
stopifnot(rownames(cna_mat_list[[1]]) == rownames(cna_mat_list[[2]]))
stopifnot(colnames(cna_mat_list[[1]]) == colnames(cna_mat_list[[2]]))

stopifnot(rownames(cna_mat_list[[1]]) == rownames(mut_type[[1]]))
stopifnot(colnames(cna_mat_list[[1]]) == colnames(mut_type[[1]]))

# CNA + Oncotator mutation calls
oncoPrint(c(
    # cna_mat_list[c("deletion0", "amplification5")],
    cna_mat_list[c("deletion0", "deletion1", "amplification5")],
    mut_type[c("Missense_Mutation", "Nonsense_Mutation", "Indel_Frameshift", "Indel_Inframe", "germline")]
    ),
    alter_fun = list(
        background = function(...) NULL,  # background box for each tile
        deletion0 = function(x, y, w, h) grid.rect(
            x, y, w*0.9, h*0.9, 
            gp=gpar(fill=box_col["deletion0"], col=outline_col)),
        deletion1 = function(x, y, w, h) grid.rect(
            x, y, w*0.9, h*0.9, 
            gp=gpar(fill=box_col["deletion1"], col=outline_col)),
        # amplification34 = function(x, y, w, h) grid.rect(
        #     x, y, w*0.9, h*0.9, 
        #     gp=gpar(fill=box_col["amplification34"], col=outline_col)),
        amplification5 = function(x, y, w, h) grid.rect(
            x, y, w*0.9, h*0.9, 
            gp=gpar(fill=box_col["amplification5"], col=outline_col)),
        Missense_Mutation = function(x, y, w, h) grid.rect(
            x, y, w*0.9, h*0.25,
            gp=gpar(fill=box_col["Missense_Mutation"], col=outline_col)),
        germline = function(x, y, w, h) grid.rect(
            x, y, w*0.9, h*0.1,
            gp=gpar(fill=box_col["germline"], col=outline_col)),
        Nonsense_Mutation = function(x, y, w, h) grid.rect(
            x, y + h*0.25, w*0.9, h*0.25,
            gp=gpar(fill=box_col["Nonsense_Mutation"], col=outline_col)),
        Indel_Frameshift = function(x, y, w, h) grid.rect(
            x, y - h*0.25, w*0.9, h*0.25,
            gp=gpar(fill=box_col["Indel_Frameshift"], col=outline_col)),
        Indel_Inframe = function(x, y, w, h) grid.rect(
            x, y, w*0.4, h*0.9,
            gp=gpar(fill=box_col["Indel_Inframe"], col=outline_col)),
        Intron = function(x, y, w, h) grid.rect(
            x, y, w*0.25, h*0.25,
            gp=gpar(fill=box_col["Intron"], col=outline_col)),
        '5Flank' = function(x, y, w, h) grid.rect(
            x, y + h*0.25, w*0.25, h*0.25,
            gp=gpar(fill=box_col["5Flank"], col=outline_col)),
        '3UTR' = function(x, y, w, h) grid.rect(
            x, y - h*0.25, w*0.25, h*0.25,
            gp=gpar(fill=box_col["3UTR"], col=outline_col))
    ),
    col=box_col,  # legend?
    top_annotation = HeatmapAnnotation(
        cbar=anno_oncoprint_barplot(),
        patient=patient_ids,
        purity=opt_clust$purity,
        col=list(
            patient=patient_colors,
            purity=colorRamp2(
                range(opt_clust$purity),
                c("white", brewer.pal(9, "Greens")[8])
            )
        )
    )
)
dev.off()
