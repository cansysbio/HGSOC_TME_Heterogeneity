scatterPlot = function(x, y,
	alternative="two.sided", corner="topleft", col="black", method="pearson", pch=16, ...) {
	require(lmtest)

	# Print test for normality
	print(shapiro.test(x))
	print(shapiro.test(y))

	n = sum(!is.na(x) & !is.na(y))

	cor_test = cor.test(x, y, alternative=alternative, method=method)

	plot(x, y,
		pch=pch,
		col=col,
		main=method,
		...)

	legend(corner, legend=c(
			paste0("r=", format(cor_test$estimate, digits=3)),
			paste0("P=", format(cor_test$p.value, digits=3)),
			paste0("n=", n)
		),
		bty="n"
	)

	# Print test
	fit = lm(y~x)

	# Breusch-Pagan test of heteroscedasity
	print(bptest(fit))

	abline(fit, col=col)
}


scatterPlotJit = function(x, y,
	jit_amount=1,
	alternative="two.sided", corner="topleft", col="black", method="pearson", pch=16, ...) {
	cor_test = cor.test(x, y, alternative=alternative, method=method)

	require(lmtest)

	# Print test for normality
	print(shapiro.test(x))
	print(shapiro.test(y))

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

	# Breusch-Pagan test of heteroscedasity
	print(bptest(fit))

	abline(fit, col="black")
}


# Scatter plot with confidence intervals on x-axis
scatterPlotCI = function(x, y, x_ci,
	alternative="two.sided", corner="topright", col="black", method="pearson", pch=21, bg="white", ...) {

	require(lmtest)

	# Print test for normality
	print(shapiro.test(x))
	print(shapiro.test(y))

	n = sum(!is.na(x) & !is.na(y))

	cor_test = cor.test(x, y, alternative=alternative, method=method)

	plot(x, y,
		# col=col,
		# xlim=range(c(0, x), na.rm=TRUE),
		xlim=c(0, max(x, na.rm=TRUE) * 1.1),
		type="n",
		main=method,
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
			paste0("P=", format(cor_test$p.value, digits=3)),
			paste0("n=", n)
		),
		bty="n"
	)

	fit = lm(y~x)

	# Breusch-Pagan test of heteroscedasity
	print(bptest(fit))

	abline(fit, col=bg)
}
