rm(list=ls())

library(stringr)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

# TCGA whole exome data downloaded by Danish Memon, stored at
# /data/memon01/tcga/TCGA-OV
# including the 'manifest.txt' file describing available .bam files

samples = read.table("data/TCGA/manifest.txt",
	stringsAsFactors=FALSE,
	header=TRUE)
# samples$filename


# Parse file names and add to sample table
file_name_annot = data.frame(
	str_split_fixed(
		samples$filename,
		"[.]",  # separator
		4),  # number of extracted columns
	stringsAsFactors=FALSE)
colnames(file_name_annot) = c("group", "tcga_id", "file_info", "file_ext")  # name columns

samples = cbind(samples, file_name_annot)


# Samples encoding
# samples$tcga_id

# Get codes specififying normal and cancer samples.
# Names with -01 correspond to primary tumour while those -10 (blood) or -11 (normal tissue) are normal.

# Get 4th element of sample ID
samples$normal_cancer_code = sapply(
	strsplit(samples$tcga_id, "-"),
	function(x) x[4])

samples$portion_analyte = sapply(
	strsplit(samples$tcga_id, "-"),
	function(x) x[5])

samples$analyte = substr(samples$portion_analyte, 3, 4)


samples$tcga_patient_id = str_match(samples$tcga_id, "^(TCGA-[0-9]+-[0-9]+)")[, 2]

samples$match_group = paste0(samples$tcga_patient_id, "_", samples$analyte)


sample_pairs = sapply(unique(samples$match_group), function(group) {
	idx = samples$match_group == group

	primary_bam = samples$tcga_id[idx & samples$normal_cancer_code == "01A"][1]  # primary tumor
	blood_bam = samples$tcga_id[idx & samples$normal_cancer_code == "10A"][1]  # blood

	file_pair = c(primary_bam,  blood_bam)

	return(file_pair)
})

sample_pairs = t(sample_pairs)
colnames(sample_pairs) = c("tumor", "blood")
sample_pairs = data.frame(sample_pairs)

write.csv(sample_pairs, file="data/TCGA/sample_pairs.csv")

sum(complete.cases(sample_pairs))


# Some examples
samples[samples$tcga_patient_id == "TCGA-04-1331", ]
samples[samples$tcga_patient_id == "TCGA-13-0757", ]
samples[samples$tcga_patient_id == "TCGA-04-1335", ]
samples[samples$tcga_patient_id == "TCGA-10-0937", ]
samples[samples$tcga_patient_id == "TCGA-10-0931", ]

samples[samples$match_group == "TCGA-10-0936_W", ]
