
options(java.parameters = "-Xmx6000m")
library(GreedyExperimentalDesign)
library(xtable)

#tab 1 results for 2n=100
n = 100
p = 1
r = 20
X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, diagnostics = TRUE)
res = resultsGreedySearch(rd, max_vectors = NULL)
res$obj_vals

initial_obj_vals = array(NA, r)
for (j in 1 : r){
  initial_obj_vals[j] = compute_objective_val(X, res$starting_indicTs[, j])
}
initial_obj_vals

print(xtable(cbind(initial_obj_vals, res$num_iters, res$obj_vals), digits = c(0,2,0,7)), include.rownames = FALSE)

#now demonstrate optimal at 2n=28
n = 28
p = 1
r = 5
X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, diagnostics = TRUE)
res = resultsGreedySearch(rd, max_vectors = NULL)
res$obj_vals

initial_obj_vals = array(NA, r)
for (j in 1 : r){
  initial_obj_vals[j] = compute_objective_val(X, res$starting_indicTs[, j])
}

greedy_res = cbind(initial_obj_vals, res$num_iters, res$obj_vals)

NUM_CORES = 3
oed = initOptimalExperimentalDesignObject(X, num_cores = NUM_CORES, objective = "abs_sum_diff", wait = TRUE)
opt_res = resultsOptimalSearch(oed)$obj_val

print(xtable(rbind(greedy_res, c(NA, NA, opt_res)), digits = c(0,2,0,9)), include.rownames = FALSE)

#fig 1
log_ns = seq(1, 2.5, by = 0.25)
sim_res = matrix(NA, nrow = length(log_ns), ncol = 3)
p = 1
r = 1000

for (i in 1 : length(log_ns)){
  n = 2 * round(10^(log_ns[i]) / 2) * 2 ## paper 2n = 100
  cat("sim for n = ", n, "\n")
  X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
  rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE)
  res = resultsGreedySearch(rd, max_vectors = NULL)
  sim_res[i, 1] = log(mean(res$obj_vals)) / log(10)
  
  X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rexp)
  rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE)
  res = resultsGreedySearch(rd, max_vectors = NULL)
  sim_res[i, 2] = log(mean(res$obj_vals)) / log(10)
  
  X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = runif)
  rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE)
  res = resultsGreedySearch(rd, max_vectors = NULL)
  sim_res[i, 3] = log(mean(res$obj_vals)) / log(10)
}

plot(log_ns, sim_res[, 1], col = "blue", type = "l", 
     xlab = "n",
     ylab = "log10(b)",
     xaxt = "n",
     ylim = c(min(sim_res), max(sim_res)), lwd = 2)
points(log_ns, sim_res[, 2], col = "red", type = "l", lty = 2, lwd = 3)
points(log_ns, sim_res[, 3], col = "darkgreen", type = "l", lty = 3, lwd = 4)
axis(1, at = log_ns[c(1,3,5,7)], labels = round(10^log_ns)[c(1,3,5,7)])

#fig 2
log_ns = seq(1, 2.5, by = 0.25)
ps = c(1, 2, 5, 10, 40)
sim_res = matrix(NA, nrow = length(log_ns), ncol = length(ps))
switches_res = matrix(NA, nrow = length(log_ns), ncol = length(ps))
r = 1000

for (i in 1 : length(log_ns)){
  n = 2 * round(10^(log_ns[i]) / 2) * 2 #2n in the paper
  for (j in 1 : length(ps)){
    p = ps[j]
    X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
    rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, diagnostics = TRUE)
    res = resultsGreedySearch(rd, max_vectors = NULL)
    sim_res[i, j] = log(mean(res$obj_vals)) / log(10)
    switches_res[i, j] = mean(res$num_iters)
  }
}

#fig 2a
par(mar = c(5,5,1,5))
plot(log_ns, log_ns, type = "n", 
     xlab = "n",
     xaxt = "n",
     ylab = "log10(b)",
     ylim = c(min(sim_res), max(sim_res)))
for (j in 1: length(ps)){
  points(log_ns, sim_res[, j], type = "l")
}
axis(side = 4, at = sim_res[length(log_ns), ], labels = ps)
axis(1, at = log_ns[c(1,3,5,7)], labels = round(10^log_ns)[c(1,3,5,7)])
mtext("p", side = 4, padj = 4)

#fig 2b
plot(log_ns, log_ns, type = "n", 
     xlab = "n",
     xaxt = "n",
     ylab = "average # of switches",
     ylim = c(min(switches_res), max(switches_res)))
for (j in 1: length(ps)){
  points(log_ns, switches_res[, j], type = "l")
}
axis(side = 4, at = switches_res[length(log_ns), ], labels = ps)
axis(1, at = log_ns[c(1,3,5,7)], labels = round(10^log_ns)[c(1,3,5,7)])
mtext("p", side = 4, padj = 4)

#table 3
log_ns = seq(1, 2.5, by = 0.25)
r = 100 # for speed
Xregr = data.frame(matrix(NA, nrow = 0, ncol = 3))
for (i in 1 : length(log_ns)){
  n = 2 * round(10^(log_ns[i]) / 2) * 2 #2n is the parameterization in the paper but not the R package
  for (j in 1 : length(ps)){
    p = ps[j]
    X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
    rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, diagnostics = TRUE)
    res = resultsGreedySearch(rd, max_vectors = NULL)
    Xregr = rbind(Xregr, cbind(log(res$obj_vals), log(n / 2), p^-1)) #n/2 because 2n is how it is parameterized in the paper
  }
}

colnames(Xregr) = c("logb", "logn", "invp")
mod = lm(logb ~ logn * invp, Xregr)
summary(mod)
library(stargazer)
stargazer(mod)

#fig 3
log_ns = seq(1, 2.5, by = 0.25)
ps = c(1, 2, 5, 10, 40)
entropy_res = matrix(NA, nrow = length(log_ns), ncol = length(ps))
norm_se_res = matrix(NA, nrow = length(log_ns), ncol = length(ps))
max_eigenval_min_one_over_n = matrix(NA, nrow = length(log_ns), ncol = length(ps))
r = 5000

for (i in 1 : length(log_ns)){
  n = 2 * round(10^(log_ns[i]) / 2)
  for (j in 1 : length(ps)){
    p = ps[j]
    X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
    rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE)
    res = resultsGreedySearch(rd, max_vectors = NULL)
    designs = res$ending_indicTs
    metr = compute_randomization_metrics(designs)
    entropy_res[i, j] = metr$rand_entropy_metric
    norm_se_res[i, j] = metr$rand_norm_se_metric
    max_eigenval_min_one_over_n[i, j] = (metr$max_eigenval - 1) / n
  }
}

#need metrics for complete randomization
cr_entropy = array(NA, length(log_ns))
cr_norm_se = array(NA, length(log_ns))
cr_max_eigenval_minus_one_over_n = array(NA, length(log_ns))
for (i in 1 : length(log_ns)){
  n = 2 * round(10^(log_ns[i]) / 2)
  X = generate_stdzied_design_matrix(n = n, p = 1)
  rd = initRerandomizationExperimentalDesignObject(X, r, wait = TRUE)
  designs = sapply(.jcall(rd$java_obj, "[[I", "getEndingIndicTs"), .jevalArray)
  metr = compute_randomization_metrics(designs)
  cr_entropy[i] = metr$rand_entropy_metric
  cr_norm_se[i] = metr$rand_norm_se_metric
  cr_max_eigenval_minus_one_over_n[i] = (metr$max_eigenval - 1) / n
}


par(mar = c(5,5,1,1))

#fig 3a
plot(log_ns, log_ns, type = "n", 
     xlab = "n",
     xaxt = "n",
     ylab = "E",
     ylim = c(-0.02, 1))
for (j in 1: length(ps)){
  points(log_ns, entropy_res[, j], type = "l")
}
points(log_ns, cr_entropy, type = "l", lty = 2, col = "gray")
axis(1, at = log_ns[c(1,3,5,7)], labels = round(10^log_ns)[c(1,3,5,7)])
label_locs = entropy_res[1, ] - 0.03 #needs to be adjusted each sim
label_locs[1] = label_locs[1] + 0.06 #needs to be adjusted each sim
text(x = log_ns[1], y = label_locs, labels = ps)

#fig 3b
plot(log_ns, log_ns, type = "n", 
     xlab = "n",
     ylab = "D",
     xaxt = "n",
     ylim = c(0, 1.02))
for (j in 1: length(ps)){
  points(log_ns, norm_se_res[, j], type = "l")
}
points(log_ns, cr_norm_se, type = "l", lty = 2, col = "gray")
axis(1, at = log_ns[c(1,3,5,7)], labels = round(10^log_ns)[c(1,3,5,7)])
label_locs = norm_se_res[1, ] + 0.03 #needs to be adjusted each sim
label_locs[1] = norm_se_res[1] - 0.03 #needs to be adjusted each sim
text(x = log_ns[1], y = label_locs, labels = ps)


#fig 3c
plot(log_ns, log_ns, type = "n", 
     xlab = "n",
     ylab = "(maximum eigenvalue - 1) / n",
     xaxt = "n",
     ylim = c(0, max(max_eigenval_min_one_over_n) + 0.02))
for (j in 1: length(ps)){
  points(log_ns, max_eigenval_min_one_over_n[, j], type = "l")
}
points(log_ns, cr_max_eigenval_minus_one_over_n, type = "l", lty = 2, col = "gray")
axis(1, at = log_ns[c(1,3,5,7)], labels = round(10^log_ns)[c(1,3,5,7)])
label_locs = max_eigenval_min_one_over_n[1, ] + 0.025 #needs to be adjusted each sim
label_locs[1] = max_eigenval_min_one_over_n[1] - 0.03 #needs to be adjusted each sim
text(x = log_ns[1], y = label_locs, labels = ps)


#Fig 4
options(java.parameters = "-Xmx20000m")
library(GreedyExperimentalDesign)
n = 100
p = 1
r = 5e6
X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
rd = initGreedyExperimentalDesignObject(X, r, wait = FALSE)

plotted = FALSE

while (TRUE){
	res = resultsGreedySearch(rd)$obj_vals_unordered
	l = length(res)
	cat("num done: ", l, " i.e. ", l / r * 100, "%\n")
	
	xs = seq(1, l, r / 100)
	ys = array(NA, length(xs))
	for (i in 1 : length(xs)){
		ys[i] = (min(res[1 : xs[i]]))
	}
	plot(log(xs), log(ys), type = "l")
	
	if (l == r){
		break
	}
	Sys.sleep(5)
}







####Gurobi testing stuff
options(java.parameters = "-Xmx4000m")
NUM_CORES = 3
library(GreedyExperimentalDesign)
.jaddClassPath("/gurobi752/win64/lib/gurobi.jar")


#sample a covariate matrix
n = 50
p = 1
X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
r = 100

#numerical optimization
gnoed = initGurobiNumericalOptimizationExperimentalDesignObject(X, time_limit_min = 0.5, num_cores = NUM_CORES)
indicT = resultsGurobiNumericalOptimizeExperimentalDesign(gnoed)$indicT
compute_objective_val(X, indicT, objective = "mahal_dist")

random_indices = sample(1 : n)
X_randomized = X[random_indices, , drop = FALSE]
gnoed = initGurobiNumericalOptimizationExperimentalDesignObject(X_randomized, time_limit_min = 0.5, num_cores = NUM_CORES)
indicT = resultsGurobiNumericalOptimizeExperimentalDesign(gnoed)$indicT
compute_objective_val(X_randomized, indicT, objective = "mahal_dist")
# compute_objective_val(X, indicT[order(random_indices)], objective = "mahal_dist")

indicTs = gurobi_multiple_designs(X, r, time_limit_min = 0.5, num_cores = NUM_CORES)
for (i in 1 : nrow(indicTs)){
  indicTs[i, ] = 2 * (indicTs[i, ] - 0.5)
}
var_cov_matrix_est = var(indicTs)
abs(var_cov_matrix_est[1:10, 1:10])
eigen(var_cov_matrix_est)$values[1]





#greedy pair switching
rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, objective = "mahal_dist", num_cores = 4)
res = resultsGreedySearch(rd, max_vectors = r)
indicTs = t(res$ending_indicTs)
for (i in 1 : nrow(indicTs)){
  indicTs[i, ] = 2 * (indicTs[i, ] - 0.5)
}
var_cov_matrix_est = var(indicTs)
abs(var_cov_matrix_est[1:10, 1:10])
eigen(var_cov_matrix_est)$values[1]
# res$obj_vals


indicTs = complete_randomization_with_balanced_ns(n, r)
var_cov_matrix_est = var(indicTs)
abs(var_cov_matrix_est[1:10, 1:10])
eigen(var_cov_matrix_est)$values[1]


indicTs = complete_randomization(n, r)
var_cov_matrix_est = var(indicTs)
abs(var_cov_matrix_est[1:10, 1:10])
eigen(var_cov_matrix_est)$values[1]
#optimal solutions
# oed = initOptimalExperimentalDesignObject(X, num_cores = NUM_CORES, objective = "mahal_dist", wait = TRUE)
# opt_res = resultsOptimalSearch(oed)$obj_val
# opt_indicT = resultsOptimalSearch(oed)$indicT
# #calc Mahal
# XT_bar = colMeans(X[opt_indicT == 1, , drop = FALSE])
# XC_bar = colMeans(X[opt_indicT == 0, , drop = FALSE])
# XT_bar_min_XC_bar = XT_bar - XC_bar
# XT_bar_min_XC_bar %*% gnoed$SinvX %*% XT_bar_min_XC_bar



complete_randomization = function(n, r){
  indicTs = matrix(NA, nrow = r, ncol = n)
  
  for (nsim in 1 : r){
    indicTs[nsim, ] = 2 * (rbinom(n, 1, 0.5) - 0.5)
  }
  indicTs
}

complete_randomization_with_balanced_ns = function(n, r){
  indicTs = matrix(NA, nrow = r, ncol = n)
  
  for (nsim in 1 : r){
    indicTs[nsim, ] = sample(c(rep(1, n / 2), rep(-1, n / 2)))
  }
  indicTs
}


