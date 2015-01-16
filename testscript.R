setwd("C:\\Users\\kapelner\\workspace\\GreedyExperimentalDesign")
library(roxygen2)
roxygenise("GreedyExperimentalDesign", clean = TRUE)

options(java.parameters = "-Xmx1000m")
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

max_designs = 10
ged = initGreedyExperimentalDesignObject(X, max_designs = max_designs, num_cores = 5)

startGreedySearch(ged)
ged
#stopGreedySearch(ged)

plot(ged)

res = resultsGreedySearch(ged, max_vectors = 100)
hist(res$obj_vals, br = max_designs / 10)

obj_vals_rand_order = res$obj_vals[order(rnorm(max_designs))]
obj_vals_rand_order_mins = array(NA, max_designs)
for (d in 1 : max_designs){
	obj_vals_rand_order_mins[d] = min(obj_vals_rand_order[1 : d])
}
plot(1 : max_designs, obj_vals_rand_order_mins)


