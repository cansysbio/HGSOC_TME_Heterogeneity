# Test data structure for various dependencies of the segment visualization
testSegments = function(segs) {
	# Any x or y chromosomes?
	# ----------------------------------
	non_numeric_chr = lapply(segs, function(s) {
		idx = !s$Chromosome %in% 1:22
		unique(s$Chromosome[idx])
	})

	non_numeric_chr = unique(unlist(non_numeric_chr))

	if (length(non_numeric_chr) > 0) {
		warning("Non numeric chromosome: ", non_numeric_chr)
	}
}
