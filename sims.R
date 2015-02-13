options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

num_diff_datasets = 10
num_reps_per_dataset = 10
ns = c(50, 100, 200, 400, 1000)

###do sims for p = 1
ps = c(1, 2, 5, 10) #, 5, 10, 50
max_iterss = c(1, 2, Inf)

all_results = data.frame(matrix(NA, nrow = 0, ncol = 5))

for (i in 1 : length(ns)){
	for (max_iters in max_iterss){
		n = ns[i]
		
		n_col = rep(n, num_reps_per_dataset)
		for (p in ps){
			cat ("n =", n, "p =", p, "max_iters =", max_iters, "\n")
			p_col = rep(p, num_reps_per_dataset)
			mi_col = rep(max_iters, num_reps_per_dataset)
			for (j in 1 : num_diff_datasets){
				X = generate_stdzied_design_matrix(n, p)
				ged = initGreedyExperimentalDesignObject(X, max_designs = num_reps_per_dataset, num_cores = 4, wait = TRUE, max_iters = max_iters)
				startGreedySearch(ged)
				res = resultsGreedySearch(ged, max_vectors = 0)
#				print(res$obj_vals)
				all_results = rbind(all_results, cbind(n_col, p_col, rep(j, num_diff_datasets), mi_col, res$obj_vals))
			}	
			write.csv(all_results, "all_results.csv")
		}
	}
}
colnames(all_results) = c("n", "p", "dataset", "max_iters", "val")

x = all_results$n
y = all_results$val
ln_x = log(x)
ln_y = log(y)
par(mfrow = c(2, 1))
plot(x, y)
plot(ln_x, ln_y)
mod = lm(ln_y ~ ln_x)
lm_res = coef(summary(mod))
mod = lm(ln_y ~ poly(ln_x, 4))
lm_res = coef(summary(mod))