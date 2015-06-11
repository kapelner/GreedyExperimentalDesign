
options(java.parameters = "-Xmx10000m")
library(GreedyExperimentalDesign)
#
#X = generate_stdzied_design_matrix(n = 10, p = 1)
#oed = initOptimalExperimentalDesignObject(X, num_cores = 4, objective = "abs_sum_diff", wait = TRUE)
#res = resultsOptimalSearch(oed)
#res
#
#X = generate_stdzied_design_matrix(n = 20, p = 1)
#oed = initOptimalExperimentalDesignObject(X, num_cores = 4, objective = "abs_sum_diff", wait = TRUE)
#res = resultsOptimalSearch(oed)
#res
#
#X = generate_stdzied_design_matrix(n = 22, p = 1)
#oed = initOptimalExperimentalDesignObject(X, num_cores = 4, objective = "abs_sum_diff", wait = TRUE)
#res = resultsOptimalSearch(oed)
#res
#
#X = generate_stdzied_design_matrix(n = 26, p = 1)
#oed = initOptimalExperimentalDesignObject(X, num_cores = 4, objective = "abs_sum_diff", wait = TRUE)
#res = resultsOptimalSearch(oed)
#res

NUM_CORES = 3
n = 26
ps = c(1, 2, 5, 10, 20)
rs = c(1, 5, 10, 100, 1000, 10000)
opt_obj_vals = array(NA, length(ps))
greedy_obj_vals = list() 
for (r in 1 : length(rs)){
	greedy_obj_vals[[r]] = array(NA, length(ps))
}


for (i in 1 : length(ps)){
	X = generate_stdzied_design_matrix(n = n, p = ps[i])
	
	for (r in 1 : length(rs)){
		ged = initGreedyExperimentalDesignObject(X, max_designs = rs[r], num_cores = NUM_CORES, wait = TRUE)
		greedy_obj_vals[[r]][i] = resultsGreedySearch(ged, max_vectors = 0)$obj_vals[1]
	}
	
	oed = initOptimalExperimentalDesignObject(X, num_cores = NUM_CORES, objective = "abs_sum_diff", wait = TRUE)
	opt_obj_vals[i] = resultsOptimalSearch(oed)$obj_val
}

log_greedy_obj_vals = log(greedy_obj_vals) / log(10)
log_opt_obj_vals = log(opt_obj_vals) / log(10)
log_ps = log(ps) / log(10)

plot(log_ps, log(greedy_obj_vals[[1]]) / log(10), ylim = c(min(log_opt_obj_vals), max(log(greedy_obj_vals[[1]]) / log(10))),
		type = "o", col = "blue", 
		ylab = "log10 obj function", xlab = "log10 p", main = paste("greedy switch vs optimal for n =", n))
points(log_ps, log_opt_obj_vals, type = "o", col = "green")
for (r in 1 : length(rs)){
	points(log_ps, log(greedy_obj_vals[[r]]) / log(10), type = "o", col = "blue")
}



#look at jsut optimal
rep = 4
ns = seq(from = 6, to = 26, by = 2)
#opt_obj_vals_b = matrix(NA, nrow = 0, ncol = 2)
for (i in 1 : length(ns)){
	for (r in 1 : rep){
		X = generate_stdzied_design_matrix(n = ns[i], p = 1)
		
		oed = initOptimalExperimentalDesignObject(X, num_cores = NUM_CORES, objective = "abs_sum_diff", wait = TRUE)
		opt_obj_vals_b = rbind(opt_obj_vals_b, c(ns[i], resultsOptimalSearch(oed)$obj_val))
	}
}

plot(opt_obj_vals_b[, 1], log(opt_obj_vals_b[, 2]) / log(10),
		ylab = "log10 obj function", xlab = "n", main = paste("greedy switch vs optimal for n =", n))

colnames(opt_obj_vals_b) = c("n", "obj_val")
y = log(opt_obj_vals_b[, 2]) / log(10)
x = opt_obj_vals_b[, 1]
mod = lm(y ~ x)
abline(mod)
summary(mod)
