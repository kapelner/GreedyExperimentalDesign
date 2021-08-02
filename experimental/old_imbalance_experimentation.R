options(java.parameters = "-Xmx2g")
# devtools::install_github("mkuhn/dict")
pacman::p_load(GreedyExperimentalDesign, tidyverse, lobstr, data.table, magrittr, dict, microbenchmark, Rcpp, RcppArmadillo, pkgbuild, nbpMatching)
# pkgbuild::find_rtools(TRUE)


n = 60
p = 1
X = matrix(rnorm(n * p), ncol = p)
X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})

WBCRD = complete_randomization_with_forced_balanced(n, 10000)
mean(colSums(abs(t(X) %*% t(WBCRD) / (n / 2))))


bin_exp_des = binaryMatchExperimentalDesignSearch(X)
# bin_exp_des$indices_pairs
X1 = X[bin_exp_des$indicies_pairs[, 1], , drop = FALSE]
X2 = X[bin_exp_des$indicies_pairs[, 2], , drop = FALSE]
# x_pairs = cbind(, )
colavg1 = colMeans(X1)
colavg2 = colMeans(X2)
sum(abs(colavg1 - colavg2)) #balance for an original optimal nonbipartite matching set

pair_diffs = X1 - X2


nvec = 100
greedy_exp_des = initGreedyExperimentalDesignObject(pair_diffs, objective = "abs_sum_diff", max_designs = nvec, wait = TRUE, diagnostics = TRUE)
greedy_exp_res = resultsGreedySearch(greedy_exp_des, max_vectors = nvec, form = "pos_one_min_one")
W = greedy_exp_res$ending_indicTs

res = array(NA, nvec)
for (i in 1 : nvec){
  X1a = rbind(X1[W[i, ] == 1, , drop = FALSE], X2[W[i, ] == -1, , drop = FALSE])
  X2a = rbind(X2[W[i, ] == 1, , drop = FALSE], X1[W[i, ] == -1, , drop = FALSE])
  colavg1 = colMeans(X1a)
  colavg2 = colMeans(X2a)
  res[i] = sum(abs(colavg1 - colavg2) / sdX)
}
mean(res)
summary(unlist(lapply(greedy_exp_res$switches, nrow)))



