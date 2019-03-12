# log ratio test taking into account the compositional nature of the copy-number signature fractions.
rm(list=ls())

library(data.table)
library(compositions)
# library(RColorBrewer)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../")
#setwd("~/GoogleDrive/projects/cambridge/ovarianCancerHeterogeneityChemo/repo/HGSOC_TME_Heterogeneity/4")

source("lib/parse.R")
source("lib/plots.R")

titan = loadOptClust("data/titanCNA_optimalClusterSolutions.csv", sep=",")

patient_colors = loadPatientColors(titan)

cn_sig = readRDS("data/miller_cnsignatures.rds")

# Test if tumor labels are aligned
stopifnot(all(titan$barcode == colnames(cn_sig)))

# One vs all log ratios
log_ratio = sapply(1:nrow(cn_sig), function(k) {
	log10(cn_sig[k, ] / apply(cn_sig[-k, ], 2, sum))
})

log_ratio[is.infinite(log_ratio)] = NA
log_ratio = t(log_ratio)

#pdf("plots/copyNumberSignaturesCompositional.pdf", height=2.8)
par(mfrow=c(1, 3))
for (k in c(1, 6, 4)) {
  scatterPlot(log_ratio[k, ], titan$purity,
              xlab=paste0("log s", k, " / s"),
              ylab="Purity (TitanCNA)"
  )
}
#dev.off()

ilr_transformation <- data.frame(ilr_1=as.vector(ilr(acomp(t(cn_sig)), V = c(1, rep(-1, 6)))),
                                 ilr_4=as.vector(ilr(acomp(t(cn_sig)), V = c(rep(-1, 3), 1, rep(-1, 3)))),
                                 ilr_6=as.vector(ilr(acomp(t(cn_sig)), V = c(rep(-1, 5), 1, -1))))
head(ilr_transformation)

pdf("plots/copyNumberSignaturesCompositional2.pdf", height=2.8)
par(mfrow=c(1, 3))
for (k in c(1,4,6)) {
  scatterPlot(ilr_transformation[,paste0('ilr_', k)], titan$purity,
              xlab=paste0("ilr with binary partition\nof (", k, "|",
                          paste0(which(!(1:7) %in% k),collapse=','), ")"),
              ylab="Purity (TitanCNA)"
  )
}
dev.off()

#' Explanation: the binary partition indicates that we are comparing the exposition to
#' a given signature (1,4,6) to the geometric mean of all others by their log-ratio
#' (i.e. log( (exposure to 4/) / (geometric mean of exposure to others) )