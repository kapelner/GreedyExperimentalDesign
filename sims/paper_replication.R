
options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)
library(xtable)

#tab 1
n = 20
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

#tab 1a
print(xtable(cbind(initial_obj_vals, res$num_iters, res$obj_vals), digits = c(0,2,0,5)), include.rownames = FALSE)

#tab 1b
TOP_N_VECS = 4
demo = matrix(NA, nrow = n, ncol = 0)
demo = cbind(X, demo)
for (tn in 1 : TOP_N_VECS){
  t0 = res$starting_indicTs[, tn]
  demo = cbind(demo, t0)
  switches = res$switches[[tn]] + 1 #java indexes from 0...n-1
  obj_vals_iter = res$obj_val_by_iters[[tn]]
  for (s in 1 : ncol(switches)){
    switch = switches[, s]
    t_new = t0
    t_new[switch[1]] = t0[switch[2]]
    t_new[switch[2]] = t0[switch[1]]
    demo = cbind(demo, t_new)
    t1 = t_new
  }
}

#get delta b's
delta_bs = list()
for (tn in 1 : TOP_N_VECS){
  switches = res$switches[[tn]] + 1 #java indexes from 0...n-1
  obj_vals_iter = res$obj_val_by_iters[[tn]]
  delta_bs[[tn]] = -diff(c(initial_obj_vals[tn], obj_vals_iter))
}
delta_bs

print(xtable(demo, digits = c(0,2,rep(0, 11))), include.rownames = FALSE)

#no relationship --- fig not in paper
plot(initial_obj_vals, res$obj_vals)


#fig 1
log_ns = seq(1, 2.5, by = 0.25)
sim_res = matrix(NA, nrow = length(log_ns), ncol = 3)
p = 1
r = 1000

for (i in 1 : length(log_ns)){
  n = 2 * round(10^(log_ns[i]) / 2)
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
     xlab = "log10(n)",
     ylab = "log10(b)",
     ylim = c(min(sim_res), max(sim_res)))
points(log_ns, sim_res[, 2], col = "red", type = "l", lty = 2)
points(log_ns, sim_res[, 3], col = "darkgreen", type = "l", lty = 3)


#fig 2
log_ns = seq(1, 2.5, by = 0.25)
ps = c(1, 2, 5, 10, 40)
sim_res = matrix(NA, nrow = length(log_ns), ncol = length(ps))
switches_res = matrix(NA, nrow = length(log_ns), ncol = length(ps))
r = 100

for (i in 1 : length(log_ns)){
  n = 2 * round(10^(log_ns[i]) / 2)
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
     xlab = "log10(n)",
     ylab = "log10(b)",
     ylim = c(min(sim_res), max(sim_res)))
for (j in 1: length(ps)){
  points(log_ns, sim_res[, j], type = "l")
}
axis(side = 4, at = sim_res[length(log_ns), ], labels = ps)
mtext("p", side = 4, padj = 4)

#fig 2b
plot(log_ns, log_ns, type = "n", 
     xlab = "log10(n)",
     ylab = "average # of switches",
     ylim = c(min(switches_res), max(switches_res)))
for (j in 1: length(ps)){
  points(log_ns, switches_res[, j], type = "l")
}
axis(side = 4, at = switches_res[length(log_ns), ], labels = ps)
mtext("p", side = 4, padj = 4)


#fig 3
log_ns = seq(1, 2.5, by = 0.25)
ps = c(1, 2, 5, 10, 40)
entropy_res = matrix(NA, nrow = length(log_ns), ncol = length(ps))
norm_se_res = matrix(NA, nrow = length(log_ns), ncol = length(ps))
r = 1000

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
  }
}

#need metrics for complete randomization


cr_entropy = array(NA, length(log_ns))
cr_norm_se = array(NA, length(log_ns))
for (i in 1 : length(log_ns)){
  n = 2 * round(10^(log_ns[i]) / 2)
  X = generate_stdzied_design_matrix(n = n, p = 1)
  rd = initRerandomizationExperimentalDesignObject(X, r, wait = TRUE)
  designs = sapply(.jcall(rd$java_obj, "[[I", "getEndingIndicTs"), .jevalArray)
  designs = res$ending_indicTs
  metr = compute_randomization_metrics(designs)
  cr_entropy[i] = metr$rand_entropy_metric
  cr_norm_se[i] = metr$rand_norm_se_metric
}


par(mar = c(5,5,1,1))

#fig 3a
plot(log_ns, log_ns, type = "n", 
     xlab = "log10(n)",
     ylab = "E",
     ylim = c(0, 1))
for (j in 1: length(ps)){
  points(log_ns, entropy_res[, j], type = "l")
}
entropy_res[, 1]
points(log_ns, cr_entropy, type = "l", lty = 2, col = "gray")
# axis(side = 2, at = entropy_res[1, ], labels = ps, padj = 8)
# mtext("p", side = 2, padj = 8)

?axis
#fig 3b
plot(log_ns, log_ns, type = "n", 
     xlab = "log10(n)",
     ylab = "D",
     ylim = c(0, 1))
for (j in 1: length(ps)){
  points(log_ns, norm_se_res[, j], type = "l")
}
points(log_ns, cr_norm_se, type = "l", lty = 2, col = "gray")
# axis(side = 4, at = norm_se_res[length(log_ns), ], labels = ps)
# mtext("p", side = 4, padj = 4)
