
options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

NUM_CORES = 3
MINN = 4
MAXN = 26
ns = seq(from = MINN, to = MAXN, by = 2)
ps = c(1) #c(1, 2, 5, 10, 20)
rs = c(1) #c(1, 5, 10, 100, 1000, 10000)
opt_obj_vals = array(NA, length(ns))
greedy_obj_vals = array(NA, length(ns))

for (n_i in 1 : length(ns)){
	for (i in 1 : length(ps)){
		X = generate_stdzied_design_matrix(n = ns[n_i], p = ps[i])
		
		ged = initGreedyExperimentalDesignObject(X, max_designs = 100, num_cores = NUM_CORES, wait = TRUE)
		greedy_obj_vals[n_i] = resultsGreedySearch(ged, max_vectors = 0)$obj_vals[1]

		
		oed = initOptimalExperimentalDesignObject(X, num_cores = NUM_CORES, objective = "abs_sum_diff", wait = TRUE)
		opt_obj_vals[n_i] = resultsOptimalSearch(oed)$obj_val
	}	
}


log_greedy_obj_vals = log(greedy_obj_vals) / log(10)
log_opt_obj_vals = log(opt_obj_vals) / log(10)

plot(ns, log_greedy_obj_vals, ylim = c(min(log_opt_obj_vals), max(log_greedy_obj_vals)),
		type = "o", col = "blue", main = paste("greedy switch vs optimal for n =", n))
points(ns, log_opt_obj_vals, type = "o", col = "green")




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
