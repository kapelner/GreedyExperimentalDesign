
options(java.parameters = "-Xmx10000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 10, p = 1)
oed = initOptimalExperimentalDesignObject(X, num_cores = 4, objective = "abs_sum_diff", wait = TRUE)
res = resultsOptimalSearch(oed)
res

X = generate_stdzied_design_matrix(n = 20, p = 1)
oed = initOptimalExperimentalDesignObject(X, num_cores = 4, objective = "abs_sum_diff", wait = TRUE)
res = resultsOptimalSearch(oed)
res

X = generate_stdzied_design_matrix(n = 22, p = 1)
oed = initOptimalExperimentalDesignObject(X, num_cores = 4, objective = "abs_sum_diff", wait = TRUE)
res = resultsOptimalSearch(oed)
res

X = generate_stdzied_design_matrix(n = 24, p = 1)
oed = initOptimalExperimentalDesignObject(X, num_cores = 4, objective = "abs_sum_diff", wait = TRUE)
res = resultsOptimalSearch(oed)
res

#ps = c(1, 2, 5, 10, 20, 40, 100, 200, 400, 1000, 2000)
ps = c(1, 2, 5, 10, 20)
opt_obj_vals = array(NA, length(ps))
greedy_obj_vals = array(NA, length(ps))

NUM_CORES = 4
n = 20
for (i in 1 : length(ps)){
	X = generate_stdzied_design_matrix(n = n, p = ps[i])
	
	ged = initGreedyExperimentalDesignObject(X, max_designs = 1, num_cores = NUM_CORES, wait = TRUE)
	greedy_obj_vals[i] = resultsGreedySearch(ged, max_vectors = 0)$obj_vals[1]
	
	oed = initOptimalExperimentalDesignObject(X, num_cores = NUM_CORES, objective = "abs_sum_diff", wait = TRUE)
	opt_obj_vals[i] = resultsOptimalSearch(oed)$obj_val
}

plot(log(ps), log(greedy_obj_vals), ylim = c(min(log(opt_obj_vals)), max(log(greedy_obj_vals))),
		type = "o", col = "blue", 
		ylab = "log obj function", xlab = "log p", main = paste("greedy switch vs optimal for n =", n))
points(log(ps), log(opt_obj_vals), type = "o", col = "green")
