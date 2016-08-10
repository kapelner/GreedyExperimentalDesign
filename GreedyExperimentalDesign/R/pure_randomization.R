
pure_randomization = function(n, r){
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}
	indicT_0 = c(rep(n / 2, 0), rep(n / 2, 1))
	indicTs = matrix(NA, nrow = n, ncol = r)
	for (j in 1 : r){
		indicTs[, j] = sample(indicT_0)
	}
	indicTs
}