setwd("C:\\Users\\kapelner\\workspace\\GreedyExperimentalDesign")
library(roxygen2)
roxygenise("GreedyExperimentalDesign", clean = TRUE)

options(java.parameters = "-Xmx1000m")
library(GreedyExperimentalDesign)

#for debugging purposes only
generate_stdzied_design_matrix = function(n = 50, p = 3, covariate_dist = "iid_std_normal"){
	if (covariate_dist == "iid_std_uniform"){
		X = matrix(runif(n * p), nrow = n, ncol = p)	
	} else if (covariate_dist == "iid_std_normal"){
		X = matrix(rnorm(n * p), nrow = n, ncol = p)
	}
	#now standardize the matrix to make things easier later
	apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})	
}


X = generate_stdzied_design_matrix()
ged = initGreedyExperimentalDesignObject(X)
ged$java_obj
greedySearchCurrentProgress(ged)
ged

startGreedySearch(ged)




