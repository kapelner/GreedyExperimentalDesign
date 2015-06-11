options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

ns = c(50, 100, 200, 400)
num_reps = 50
karp_obj_vals = matrix(NA, nrow = 0, ncol = 2)
colnames(karp_obj_vals) = c("n", "obj_val")
karp_obj_vals_tab = matrix(NA, nrow = length(ns), ncol = num_reps)

for (i in 1 : length(ns)){
	for (j in 1 : num_reps){
		X = generate_stdzied_design_matrix(n = ns[i], p = 1)	
		ked = initKarpExperimentalDesignObject(X, wait = TRUE)
		obj_val = compute_objective_val(X, resultsKarpSearch(ked)$indicT)
		karp_obj_vals = rbind(karp_obj_vals, c(ns[i], obj_val))	
		karp_obj_vals_tab[i, j] = obj_val
	}
}

karp_obj_vals = as.data.frame(karp_obj_vals)
plot(log(karp_obj_vals$n), log(karp_obj_vals$obj_val), ylab = "log obj val", xlab = "log n")
mod = lm(log(obj_val) ~ log(n), data = karp_obj_vals)
summary(mod)
abline(mod)

###COMPARE to greedy

NUM_CORES = 4
rs = c(1, 5, 10, 30, 100, 500) #, 1000, 10000
karp_obj_vals2 = array(NA, length(ns))
greedy_obj_vals = list() 
for (r in 1 : length(rs)){
	greedy_obj_vals[[r]] = array(NA, length(ns))
}


for (i in 1 : length(ns)){
	X = generate_stdzied_design_matrix(n = ns[i], p = 1)
	ked = initKarpExperimentalDesignObject(X, wait = TRUE)
	karp_obj_vals2[i] = compute_objective_val(X, resultsKarpSearch(ked)$indicT)
	
	for (r in 1 : length(rs)){
		ged = initGreedyExperimentalDesignObject(X, max_designs = rs[r], num_cores = NUM_CORES, wait = TRUE)
		greedy_obj_vals[[r]][i] = resultsGreedySearch(ged, max_vectors = 0)$obj_vals[1]
	}
}

log_karp_obj_vals = log(karp_obj_vals2) / log(10)

plot(ns, log(greedy_obj_vals[[1]]) / log(10), ylim = c(min(log_karp_obj_vals, log(greedy_obj_vals[[length(rs)]]) / log(10)), max(log_karp_obj_vals, log(greedy_obj_vals[[1]]) / log(10))),
		type = "o", col = "blue", 
		ylab = "log10 obj function", xlab = "n", main = paste("greedy switch vs karp for values of n"))
points(ns, log_karp_obj_vals, type = "o", col = "chartreuse3")
for (r in 1 : length(rs)){
	points(ns, log(greedy_obj_vals[[r]]) / log(10), type = "o", col = ifelse(r == length(rs), "red", "blue"))
}



#options(java.parameters = "-Xmx3000m")
#library(GreedyExperimentalDesign)
#X = t(t(c(1, 3, 5, 7, 10, 20, 30, 35, 40, 45, 50, 55, 80, 100)))
#ked = initKarpExperimentalDesignObject(X, wait = TRUE)
#resultsKarpSearch(ked)
#oed = initOptimalExperimentalDesignObject(X, objective = "abs_sum_diff", wait = TRUE)
#resultsOptimalSearch(oed)
#ged = initGreedyExperimentalDesignObject(X, objective = "abs_sum_diff", max_designs = 10, wait = TRUE)
#resultsGreedySearch(ged)