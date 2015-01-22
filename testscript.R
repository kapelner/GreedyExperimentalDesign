setwd("C:\\Users\\kapelner\\workspace\\GreedyExperimentalDesign")
library(roxygen2)
roxygenise("GreedyExperimentalDesign", clean = TRUE)

options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

#for debugging purposes only
generate_stdzied_design_matrix = function(n = 50, p = 1, covariate_dist = "iid_std_normal"){
	if (covariate_dist == "iid_std_uniform"){
		X = matrix(runif(n * p), nrow = n, ncol = p)	
	} else if (covariate_dist == "iid_std_normal"){
		X = matrix(rnorm(n * p), nrow = n, ncol = p)
	}
	#now standardize the matrix to make things easier later
	apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})	
}


X = generate_stdzied_design_matrix(n = 200, p = 40)
ged = initGreedyExperimentalDesignObject(X, max_designs = 1000, num_cores = 3, objective = "abs_sum_diff")
ged

startGreedySearch(ged)
data = plot_obj_val_order_statistic(ged, order_stat = 200)
vals = data$val_order_stats
plot(1:length(vals), log(log(vals)))
plot(1:length(vals), vals^(-15))