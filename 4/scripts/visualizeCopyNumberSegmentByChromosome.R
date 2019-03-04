rm(list=ls())

library(RColorBrewer)
library(viridis)
library(dendextend)
library(data.table)

library(RCircos)  # for UCSC HG19 data
data(UCSC.HG19.Human.CytoBandIdeogram)

# Just in time compilation
library(compiler)
enableJIT(3)

setwd("/Users/koplev01/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/tests.R")

# Load segments from TitanCNA analysis
# --------------------------------

load("data/segs.RData", verbose=TRUE)

# Load TitanCNA results table
opt_clust = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")


# Load colorscheme, based on the levels of factor in opt_clust$patient_id
# -----------------------
patient_colors = loadPatientColors(opt_clust)

# Load CNA count matrix
# -------------------------------------
cna = fread("data/cna.csv")

cna_mat = data.matrix(cna[, c(-1, -2)])


# Genome-wide plot of CNA segments
# -----------------------------------------

# Get chromosome lengths
chr_length = sapply(unique(UCSC.HG19.Human.CytoBandIdeogram$Chromosome), function(chr) {
	max(UCSC.HG19.Human.CytoBandIdeogram$ChromEnd[
		UCSC.HG19.Human.CytoBandIdeogram$Chromosome == chr]
	)
})
names(chr_length) = unique(UCSC.HG19.Human.CytoBandIdeogram$Chromosome)

chr_offset = c(0, cumsum(as.numeric(chr_length)))

# Run various tests for copy-number segments
testSegments(segs)


# Hierarical clustering in gene space
dmat = dist(t(cna_mat), method="euclidean")
hc = hclust(dmat)

idx = hc$order


segs_order = segs[idx]
opt_clust_order = opt_clust[idx, ]

pdf(file="plots/copyNumberSegment.pdf", height=5, width=8)
plot(0, 0,
	type="n",
	xlim=c(0, sum(chr_length)),  # all genomic coordinates
	ylim=c(0, length(segs)),
	bty="n",
	xaxt="n", xlab="Chromosome",
	yaxt="n", ylab="Tumour"
)

poly_height = 0.8


cp_colors = c(
	rev(brewer.pal(5, "Blues")[c(4, 5)]), # 0, 1
	rgb(220, 220, 220, maxColorValue=255),
	brewer.pal(8, "Reds")[c(-1, -2)]
)


for (i in 1:length(segs_order)) {
	chr = segs_order[[i]]$Chromosome
	start = segs_order[[i]][["Start_Position.bp."]]
	end = segs_order[[i]][["End_Position.bp."]]
	cp_num = segs_order[[i]][["Copy_Number"]]

	for (j in 1:nrow(segs_order[[i]])) {
		offset = chr_offset[as.integer(chr[j])]
		genome_start = offset + start[j]
		genome_end = offset + end[j]

		polygon(
			c(genome_start, genome_end, genome_end, genome_start),
			c(i, i, i + poly_height, i + poly_height),
			border=NA,
			col=cp_colors[cp_num[j] + 1]
		)
	}

}

# Chromosome separators
for (chr_end in chr_offset) {
	abline(v=chr_end,
		col="black"
	)
}

par(xpd=TRUE)  # allow plotting outside box

# Plot chromosome numbers
for (chr in 1:(length(chr_offset) - 1)) {
	mid_point = mean(chr_offset[c(chr, chr + 1)])
	text(mid_point, 1, labels=chr, pos=1)  # below
}

# Copy number colors legend
legend("topright", legend=0:8, col=cp_colors, pch=15, xpd=TRUE)

# plot patient indicators
pt_col = patient_colors[as.integer(opt_clust_order$patient_id)]

points(
	rep(-20000000, nrow(opt_clust_order)),
	1:nrow(opt_clust_order) + poly_height / 2,
	pch=15,
	cex=0.8,
	col=pt_col)

# Labels
text(rep(-22000000, nrow(opt_clust_order)),
	1:nrow(opt_clust_order) + poly_height / 2,
	labels=opt_clust_order$barcode,
	cex=0.4,
	pos=2
)

dev.off()


# Plot dendrogram for hierarchical clustering
# -----------------------------------------------------

pdf(file="plots/copyNumberSegmentClustering.pdf", height=5, width=4)

dend = as.dendrogram(hc)

# Color labels by patient IDs
pts_col = patient_colors[as.integer(opt_clust$patient_id[idx])]

dend = dendextend::set(dend, "labels_col", pts_col)
dend = dendextend::set(dend, "labels_cex", 0.5)


par(mfrow=c(1, 2))
plot(dend, xlab="Gene copy-number distance",
	horiz=TRUE)

barplot(opt_clust$purity[idx],
	names.arg=opt_clust$barcode[idx],
	las=2,
	space=0,
	col=brewer.pal(9, "Pastel1")[4],
	xlab="Tumour purity",
	horiz=TRUE)
dev.off()
