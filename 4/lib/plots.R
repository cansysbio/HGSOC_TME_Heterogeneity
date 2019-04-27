scatterPlot = function(x, y,
	alternative="two.sided", corner="topleft", col="black", method="pearson", pch=16, ...) {
	cor_test = cor.test(x, y, alternative=alternative, method=method)

	plot(x, y,
		pch=pch,
		col=col,
		...)

	legend(corner, legend=c(
			paste0("r=", format(cor_test$estimate, digits=3)),
			paste0("P=", format(cor_test$p.value, digits=3))
		),
		bty="n"
	)

	fit = lm(y~x)
	abline(fit, col=col)
}


scatterPlotJit = function(x, y,
	jit_amount=1,
	alternative="two.sided", corner="topleft", col="black", method="pearson", pch=16, ...) {
	cor_test = cor.test(x, y, alternative=alternative, method=method)

	plot(jitter(x, amount=jit_amount), jitter(y, amount=jit_amount),
		#pch=pch,
		pch=21,
		bg=col,
		...)

	legend(corner, legend=c(
			paste0("r=", format(cor_test$estimate, digits=3)),
			paste0("P=", format(cor_test$p.value, digits=3))
		),
		bty="n"
	)

	fit = lm(y~x)
	abline(fit, col="black")
}


# Scatter plot with confidence intervals on x-axis
scatterPlotCI = function(x, y, x_ci,
	alternative="two.sided", corner="topright", col="black", method="pearson", pch=21, bg="white", ...) {

	cor_test = cor.test(x, y, alternative=alternative, method=method)

	plot(x, y,
		# col=col,
		# xlim=range(c(0, x), na.rm=TRUE),
		xlim=c(0, max(x, na.rm=TRUE) * 1.1),
		type="n",
		...)

	# CI segments
	segments(
		x - x_ci,
		y,
		x + x_ci,
		y
	)

	points(x, y, pch=pch, col=col, bg=bg)

	legend(corner, legend=c(
			paste0("r=", format(cor_test$estimate, digits=3)),
			paste0("P=", format(cor_test$p.value, digits=3))
		),
		bty="n"
	)

	fit = lm(y~x)
	abline(fit, col=bg)
}
