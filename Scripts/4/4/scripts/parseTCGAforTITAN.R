rm(list=ls())

library(stringr)
library(TCGAutils)
library(stringr)
library(dplyr)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")


# bam_set = "original"  # Original set of BAM files downloaded by Danish
bam_set = "extended"  # Extended BAM files


if (bam_set == "original") {
	# TCGA whole exome data downloaded by Danish Memon, stored at
	# /data/memon01/tcga/TCGA-OV
	# including the 'manifest.txt' file describing available .bam files

	samples = read.table(
		# "data/TCGA/manifest.txt",  # hg38
		"data/TCGA/gdc_manifest.2019-08-01.txt",  # hg19
		stringsAsFactors=FALSE,
		header=TRUE)

	samples$path = paste0("/data/memon01/tcga/TCGA-OV-LEGACY/", samples$id, "/", samples$filename)

	# write(samples$path, "data/TCGA/sample_paths.txt")  # for

} else if (bam_set == "extended") {

	samples0 = read.table(
		# "data/TCGA/manifest.txt",  # hg38
		"data/TCGA/gdc_manifest.2019-08-01.txt",  # hg19
		stringsAsFactors=FALSE,
		header=TRUE)

	samples0$path = paste0("/data/memon01/tcga/TCGA-OV-LEGACY/", samples0$id, "/", samples0$filename)

	# ADD
	samples1 = read.table(
		"data/TCGA/gdc_manifest.2019-09-13.add.txt",
		stringsAsFactors=FALSE,
		header=TRUE)

	samples1$path = paste0("/data/memon01/tcga/TCGA-OV-LEGACY-ADD/", samples1$id, "/", samples1$filename)

	# WGS
	samples2 = read.table(
		"data/TCGA/gdc_manifest.wgs.txt",
		stringsAsFactors=FALSE,
		header=TRUE)

	samples2$path = paste0("/data/memon01/tcga/TCGA-OV-LEGACY-WITH-WGS/", samples2$id, "/", samples2$filename)

	samples = rbind(samples0, samples1, samples2)
	# samples = samples2

} else {
	stop("Unrecognized bam_set value")
}


# Load tests for NCBI ('1') or USCF ('chr1') format
# From running testChrFormat.sh on list of sample paths, on mmlab cluster
# samples_info = read.csv("data/TCGA/bamfile_info.csv")
samples_info = read.csv("data/TCGA/bamfile_info_add_wgs.csv")


samples_info$reference = str_match(samples_info$header, "AS:(.*?)\\s+")[, 2]
table(samples_info$reference)


# Inner join, only include paths contained in both tables
# This ensures that bam files have been downloaded (confirmed in samples_info), and that
# they are specificed in the provided samples table
samples = merge(samples, samples_info, by="path")


# # Exclude USCF chromosome coodinates
# message("Excluding samples with UCSF chromosome coordinates, n=", sum(samples$chr_count > 0))
# samples = samples[samples$chr_count == 0, ]

# Also excluding non GRCh37 reference genomes...
idx = samples$chr_count == 0 & samples$reference %in% c("GRCh37", "GRCh37-lite")
message("Excluding samples, n=", sum(!idx))
samples = samples[idx, ]


samples = samples[!is.na(samples$path), ]

# Get TCGA barcodes from file names
samples = merge(samples,
	rename(filenameToBarcode(samples$filename, legacy=TRUE),
		filename=file_name,
		tcga_id=aliquots.submitter_id),
	by="filename")


table(table(samples$tcga_id))


# Get codes specififying normal and cancer samples.
# Names with -01 correspond to primary tumour while those -10 (blood) or -11 (normal tissue) are normal.

# Get 4th element of sample ID
samples$normal_cancer_code = sapply(
	strsplit(samples$tcga_id, "-"),
	function(x) substr(x[4], 1, 2)  # Excluding vial part, eg 01A
)

table(samples$normal_cancer_code)

samples$portion_analyte = sapply(
	strsplit(samples$tcga_id, "-"),
	function(x) x[5])

samples$analyte = substr(samples$portion_analyte, 3, 4)


# samples$tcga_patient_id = str_match(samples$tcga_id, "^(TCGA-[0-9]+-[0-9]+)")[, 2]
samples$tcga_patient_id = str_match(samples$tcga_id, "^(TCGA-[:alnum:]+-[:alnum:]+)")[, 2]  # allows letter barcodes

samples$match_group = paste0(samples$tcga_patient_id, "_", samples$analyte)


samples = samples[order(samples$normal_cancer_code), ]  # blood preference over normal

sample_pairs = sapply(unique(samples$match_group), function(group) {
	idx = samples$match_group == group

	primary_bam = samples$tcga_id[idx & samples$normal_cancer_code == "01"][1]  # primary tumor

	# normal_bam = samples$tcga_id[idx & samples$normal_cancer_code == "10"][1]  # blood only

	normal_bam = samples$tcga_id[idx & (samples$normal_cancer_code == "10" | samples$normal_cancer_code == "11")][1]  # blood or normal tissue

	file_pair = c(primary_bam,  normal_bam)

	return(file_pair)
})

sample_pairs = t(sample_pairs)
colnames(sample_pairs) = c("tumor", "normal")
sample_pairs = data.frame(sample_pairs)



# Use only complete pairs
sum(complete.cases(sample_pairs))


sample_pairs_compl = sample_pairs[complete.cases(sample_pairs), ]

included_tcga_ids = as.character(unlist(sample_pairs_compl))

samples_included = samples[match(included_tcga_ids, samples$tcga_id), ]  # picks first if duplicate IDs

sample_pairs_compl$tumor_path = samples$path[match(sample_pairs_compl$tumor, samples$tcga_id)]
sample_pairs_compl$normal_path = samples$path[match(sample_pairs_compl$normal, samples$tcga_id)]

sample_pairs_compl$tumor_ref = samples$reference[match(sample_pairs_compl$tumor, samples$tcga_id)]
sample_pairs_compl$blood_ref = samples$reference[match(sample_pairs_compl$normal, samples$tcga_id)]


if (bam_set == "original") {
	write.csv(sample_pairs, file="data/TCGA/sample_pairs.csv")  # also contains unmatched TCGA barcodes, for manual inspection

	write(
		c(sample_pairs_compl$tumor_path, sample_pairs_compl$normal_path),
		file="data/TCGA/sample_pairs_paths.txt")

	write(
		dirname(c(sample_pairs_compl$tumor_path, sample_pairs_compl$normal_path)),
		file="data/TCGA/sample_pairs_dirs.txt")

	yaml_file = "config/tcga_samples.yaml"
} else if (bam_set == "extended") {
	write.csv(sample_pairs_compl, file="data/TCGA/sample_pairs_compl_extended.csv", row.names=FALSE, quote=FALSE)

	yaml_file = "config/tcga_samples_extend.yaml"
} else {
	stop("bam_set invalid.")
}


# Write as .yaml sample file for use as input to TITAN
# yaml_file variable is set above
# -------------------------------------

yaml_indent = "  "

# File paths of all available samples
write("samples:",
	file=yaml_file)

write(
	paste0(yaml_indent, samples_included$tcga_id, ": ", samples_included$path),
	append=TRUE,
	file=yaml_file)

# Tumor-normal pairs
write("pairings:",
	append=TRUE,
	file=yaml_file)
write(
	paste0(yaml_indent, sample_pairs_compl$tumor, ": ", sample_pairs_compl$normal),
	append=TRUE,
	file=yaml_file)

# # Some examples
samples[samples$tcga_id == "TCGA-61-2610-02A-01W-1092-09", ]

# # Some examples
# samples[samples$tcga_patient_id == "TCGA-04-1331", ]
# samples[samples$tcga_patient_id == "TCGA-13-0757", ]
# samples[samples$tcga_patient_id == "TCGA-04-1335", ]
# samples[samples$tcga_patient_id == "TCGA-10-0937", ]
# samples[samples$tcga_patient_id == "TCGA-10-0931", ]

# samples[samples$match_group == "TCGA-10-0936_W", ]


# sort(samples$tcga_id)
