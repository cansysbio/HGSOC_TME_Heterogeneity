loadOptClust = function(file_path, sep="\t") {
	require(stringr)
	require(plyr)
	opt_clust = read.table(file_path, sep=sep, header=TRUE)

	# Remove erroneous comparison, from samples.yaml configuration
	# Due to input error when typing in the configuration file
	opt_clust = opt_clust[opt_clust$barcode != "RG19-N", ]

	# Extract patient IDs
	opt_clust$patient_id = str_extract(opt_clust$barcode, "^[:alpha:]+[0-9]+")
	opt_clust$patient_id = factor(opt_clust$patient_id)

	# Correct column name for mislabelled tumors
	opt_clust$tumor =opt_clust$barcode
	opt_clust$tumor = revalue(opt_clust$tumor, c("RG13T122"="RG13T12", "RG5T14Rep"="RG5T14"))

	return(opt_clust)
}

# Loads and parses Alejandro's sample annotation file
loadSampleAnnot = function(file_path) {
	require(data.table)

	# Column annotation
	annot = fread(file_path)

	# Parse
	annot$tumor = sapply(strsplit(annot$wes_label, "S"), function(x) x[1])
	annot$descr = sapply(strsplit(annot$wes_label, "S"), function(x) x[2])

	return(annot)
}

# Load and reformats Alejandro's RNA-based tumor purity scores from the RNA deconvolution method
file_path = "data/alejandro/OVCTp_log2exp_loess_norm_estimate_score.txt"
loadConvScores = function(file_path) {
	conv_scores = fread(file_path, skip=2)

	# string matrix
	temp = t(conv_scores)

	colnames(temp) = temp[1, ]  # names

	# Skip description entries
	temp = temp[c(-1, -2), ]

	labels = rownames(temp)

	# Convert to numbers
	temp = apply(temp, 2, as.numeric)

	# df = as.data.frame(temp, stringsAsFactors=FALSE)
	df = as.data.frame(temp)

	df$exp_label = labels

	return(df)
}


loadCNA = function(file_path) {
	cna = fread(file_path)

	# Correct column name for mislabelled tumor
	colnames(cna)[colnames(cna) == "RG13T122"] = "RG13T12"

	return(cna)
}

# Get correlation test statistics (r and p) from list of cor.test() outputs, allowing for NA entries		.
getCorStats = function(tests) {
	# Format output
	cor_pvals = sapply(tests, function(x) {
		out = tryCatch(x$p.value,
			error=function(err) {
				return(NA)
		})

		return(out)
	})

	cor_coeff = sapply(tests, function(x) {
		out = tryCatch(x$estimate,
			error=function(err) {
				return(NA)
		})

		return(out)
	})

	gene_symbol = names(tests)

	return(data.frame(gene_symbol=gene_symbol, r=cor_coeff, p=cor_pvals))
}

loadPatientColors = function(opt_clust) {
	require(data.table)
	require(RColorBrewer)

	# Get previous color scheme
	color_tab = fread("data/p_ks_cols.txt")
	color_tab$patient_id = paste0("RG", color_tab$case_number)

	patient_colors = color_tab$case[match(levels(opt_clust$patient_id), color_tab$patient_id)]

	# Fill in missing colors
	patient_colors[is.na(patient_colors)] = brewer.pal(8, "Set2")[-2][1:sum(is.na(patient_colors))]

	return(patient_colors)
}