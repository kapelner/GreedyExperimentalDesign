options(java.parameters = "-Xmx5000m")
library(GreedyExperimentalDesign)

NUM_CORES = 3

ns = c(50, 100, 200, 400, 1000)
ps = c(1, 5, 10, 20, 30)
rs = c(1, 10, 50, 100, 300)
DUPS = 25 #number of repeats per cell in the grid

all_results = data.frame(matrix(NA, nrow = 0, ncol = 5))

for (n in ns){
	for (p in ps){
		for (r in rs){
			for (d in 1 : DUPS){
				cat ("n =", n, "p =", p, "r =", r, "d =", d, "\n")
				#generate data independently for each run with the n and p characteristics
				X = generate_stdzied_design_matrix(n, p)
				#generate the greedy search object with r starting points and begin the search
				ged = initGreedyExperimentalDesignObject(X, max_designs = r, num_cores = NUM_CORES, wait = TRUE)
				startGreedySearch(ged)
				#get back the results
				res = resultsGreedySearch(ged, max_vectors = 0)
				#record the minimum balance found
				all_results = rbind(all_results, c(n, p, r, d, res$obj_vals[1]))				
			}					
		}
		#record the results in a CSV iteratively
		colnames(all_results) = c("n", "p", "r", "d", "val")
		write.csv(all_results, "all_results.csv")			
	}
}

### now analyze

X = read.csv("all_results.csv")
invp = X$p^-1

#model time
mod = lm(log(val) ~ log(n), data = X)
summary(mod)


mod = lm(log(val) ~ log(n) * invp, data = X)
summary(mod)

mod1 = lm(log(val) ~ log(p) + log(n) * invp, data = X)
summary(mod1)

mod2 = lm(log(val) ~ log(r) * invp + log(n) * invp, data = X)
summary(mod2)

anova(mod1, mod2)