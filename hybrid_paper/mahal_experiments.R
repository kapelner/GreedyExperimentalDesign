options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table, ggpubr, microbenchmark, checkmate)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign

set.seed(1)

n = 100
p = 2

X = matrix(NA, nrow = n, ncol = p)
X[, 1] = rnorm(n)
X[, 2] = 0.5 * X[, 1] + 0.5 * rnorm(n)
w = rep(c(1, -1), n / 2)

Sigma_X_inv = solve(var(X))

d_M_X = 1/n^2 * t(w) %*% X %*% Sigma_X_inv %*% t(X) %*% w


Delta = X[seq(from = 1, to = n, by = 2), ] - X[seq(from = 2, to = n, by = 2), ]
Sigma_Delta_inv = solve(var(Delta))

w_Delta = rep(1, n / 2)
d_M_Delta = 1 / (n / 2)^2 * t(w_Delta) %*% Delta %*% Sigma_Delta_inv %*% t(Delta) %*% w_Delta

d_M_X
d_M_Delta

t(w) %*% X
t(w_Delta) %*% Delta

Sigma_X_inv
Sigma_Delta_inv

b = binaryMatchExperimentalDesignSearch(X)
all.equal(sort(c(b$indicies_pairs)), 1 : n)

b$indicies_pairs

Delta = X[b$indicies_pairs[, 1], ] - X[b$indicies_pairs[, 2], ]
Sigma_Delta_inv = solve(var(Delta))

