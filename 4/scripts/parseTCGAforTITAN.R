rm(list=ls())

library(stringr)
library(TCGAutils)
library(stringr)

setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

# TCGA whole exome data downloaded by Danish Memon, stored at
# /data/memon01/tcga/TCGA-OV
# including the 'manifest.txt' file describing available .bam files

samples = read.table(
	# "data/TCGA/manifest.txt",  # hg38
	"data/TCGA/gdc_manifest.2019-08-01.txt",  # hg19
	stringsAsFactors=FALSE,
	header=TRUE)


samples$path = paste0("/data/memon01/tcga/TCGA-OV-LEGACY/", samples$id, "/", samples$filename)

write(samples$path, "data/TCGA/sample_paths.txt")  # for

# Load tests for NCBI ('1') or USCF ('chr1') format
# From running testChrFormat.sh on list of sample paths, on mmlab cluster
samples_chrom_format = read.csv("data/TCGA/bamfile_chr_counts.csv")

stopifnot(samples_chrom_format$path == samples$path)

samples$chrom_counts = samples_chrom_format$chr_count


# Exclude USCF chromosome coodinates
message("Excluding samples with UCSF chromosome coordinates, n=", sum(samples$chrom_counts > 0))
samples = samples[samples$chrom_counts == 0, ]


# Get TCGA barcodes from file names
samples$tcga_id = filenameToBarcode(samples$filename, legacy=TRUE)[, 3]

table(table(samples$tcga_id))


# Get codes specififying normal and cancer samples.
# Names with -01 correspond to primary tumour while those -10 (blood) or -11 (normal tissue) are normal.

# Get 4th element of sample ID
samples$normal_cancer_code = sapply(
	strsplit(samples$tcga_id, "-"),
	function(x) x[4])

table(samples$normal_cancer_code)

samples$portion_analyte = sapply(
	strsplit(samples$tcga_id, "-"),
	function(x) x[5])

samples$analyte = substr(samples$portion_analyte, 3, 4)


# samples$tcga_patient_id = str_match(samples$tcga_id, "^(TCGA-[0-9]+-[0-9]+)")[, 2]
samples$tcga_patient_id = str_match(samples$tcga_id, "^(TCGA-[:alnum:]+-[:alnum:]+)")[, 2]  # allows letter barcodes

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

write.csv(sample_pairs, file="data/TCGA/sample_pairs.csv")  # also contains unmatched TCGA barcodes, for manual inspection


# Use only complete pairs
sum(complete.cases(sample_pairs))


sample_pairs_compl = sample_pairs[complete.cases(sample_pairs), ]

included_tcga_ids = as.character(unlist(sample_pairs_compl))

samples_included = samples[match(included_tcga_ids, samples$tcga_id), ]  # picks first if duplicate IDs

sample_pairs_compl$tumor_path = samples$path[match(sample_pairs_compl$tumor, samples$tcga_id)]
sample_pairs_compl$blood_path = samples$path[match(sample_pairs_compl$blood, samples$tcga_id)]

write.csv(sample_pairs_compl, file="data/TCGA/sample_pairs_compl.csv", row.names=FALSE, quote=FALSE)


write(c(sample_pairs_compl$tumor_path, sample_pairs_compl$blood_path), file="data/TCGA/sample_pairs_paths.txt")

# Write as .yaml sample file for use as input to TITAN
# -------------------------------------
yaml_file = "config/tcga_samples.yaml"

yaml_indent = "  "

# File paths of all available samples
write("samples:",
	file=yaml_file)

write(
	# paste0(yaml_indent, samples$tcga_id, ": ", "/data/memon01/tcga/TCGA-OV-LEGACY/", samples$id, "/", samples$filename),
	paste0(yaml_indent, samples_included$tcga_id, ": ", "/data/memon01/tcga/TCGA-OV-LEGACY/", samples_included$id, "/", samples_included$filename),
	append=TRUE,
	file=yaml_file)

# Tumor-normal pairs
write("pairings:",
	append=TRUE,
	file=yaml_file)
write(
	paste0(yaml_indent, sample_pairs_compl$tumor, ": ", sample_pairs_compl$blood),
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
