options(java.parameters = "-Xmx4000m")
#install.packages("GreedyExperimentalDesign_1.0.tar.gz", repos=NULL, type= "source", lib="~")
library(GreedyExperimentalDesign)

NUM_CORES = 1

n = 100
p = 1
r = 1
ns = c(100, 200, 500, 1000, 2000, 4000, 6000, 8000, 10000)
ps = c(10)
rs = c(1)
DUPS = 2 #number of repeats per cell in the grid

all_results = data.frame(matrix(NA, nrow = 0, ncol = 5))

for (n in ns){
	for (p in ps){
		if (p >= n){
			next
		}
		for (r in rs){
			for (d in 1 : DUPS){
				cat ("n =", n, "p =", p, "r =", r, "d =", d, "\n")
				#generate data independently for each run with the n and p characteristics
				X = generate_stdzied_design_matrix(n, p)
				#generate the greedy search object with r starting points and begin the search
				ged = initGreedyExperimentalDesignObject(X, max_designs = r, num_cores = NUM_CORES, wait = TRUE)
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

library(data.table)
X = fread("all_results.csv")
dim(X)
tail(X)
X$invp = X$p^-1

#model time
mod = lm(log(val) ~ log(n), data = X[X$p == 1, ])
summary(mod)

mod = lm(log(val) ~ log(n), data = X[X$p == 5, ])
summary(mod)

mod = lm(log(val)[1:8] ~ log(n)[1:8], data = X[X$p == 10, ])
summary(mod)
plot(log(X$n), log(X$val))
abline(mod, col = "red")

mod = lm(log(val) ~ log(n), data = X[X$p == 20, ])
summary(mod)

mod = lm(log(val) ~ log(n), data = X[X$p == 30, ])
summary(mod)

mod = lm(log(val) ~ 0 + log(n) + log(n) : as.factor(p))
summary(mod)

Xp = X[X$p == 1, ]
attach(Xp)
plot(log(Xp$n), log(Xp$val))



mod = lm(log(val) ~ log(n) * log(p), data = X)
summary(mod)

mod = lm(log(val) ~ log(n) : log(p), data = X)
summary(mod)

mod4 = lm(log(val) ~ log(n) * as.factor(p) - as.factor(p) - log(n), data = X)
summary(mod4)

mod4 = lm(log(val) ~ log(n) * as.factor(invp) - invp, data = X)
summary(mod4)

mod5 = lm(log(val) ~ log(n) * invp, data = X)
summary(mod5)

anova(mod5, mod4)

mod4 = lm(log(val) ~ r * log(n) * as.factor(p) - as.factor(p) - log(n), data = X)
summary(mod4)

mod5 = lm(log(val) ~ r * as.factor(p) + log(n) * as.factor(p) - as.factor(p) - log(n), data = X)
summary(mod5)

mod4 = lm(log(val) ~ log(n) * as.factor(p), data = X)
summary(mod4)

mod1 = lm(log(val) ~ log(p) + log(n) * invp, data = X)
summary(mod1)

mod2 = lm(log(val) ~ log(r) * log(p) + log(n) * invp, data = X)
summary(mod2)

anova(mod1, mod2)

mod8 = lm(log(val) ~ log(r) + log(p) * log(n), data = X)
summary(mod8)

Xr = X[X$r == 3, ]
Xp = X[X$p == 1, ]

mod = lm(log(val) ~ log(n) * invp, data = Xpr)
summary(mod)


Xpr = X[X$p == 20 & X$r == 1, ]
hist(Xpr[Xpr$n == 100, ]$val, br = 50)
mod = lm(log(val) ~ log(n), data = Xpr)
summary(mod)
plot(log(Xpr$n), log(Xpr$val))
abline(mod)