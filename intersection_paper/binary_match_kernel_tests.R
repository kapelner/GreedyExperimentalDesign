options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table, ggpubr, microbenchmark, checkmate)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign



n = 100
p = 1

X = matrix(rnorm(n * p), nrow = n)

b = binaryMatchExperimentalDesignSearch(X)
all.equal(sort(c(b$indicies_pairs)), 1 : n)
P = b$indicies_pairs - 1
P[order(apply(P, 1, min)), ]

d = initBinaryMatchFollowedByGreedyExperimentalDesignSearch(X, max_designs = 100, wait = TRUE)
d = initGreedyExperimentalDesignSearch(X, max_designs = 100, wait = TRUE)

# ged_res = resultsGreedySearch(d$greedy_design, 100, "zero_one")
d_res = resultsBinaryMatchThenGreedySearch(d)

Pd = matrix(NA, nrow = n / 2, ncol = 2)
for (i in 1 : (n / 2)){
  Pd[i, ] = d_res$indicTs[40, d$binary_match_design$indicies_pairs[i, ]]
}
Pd
rowSums(Pd)
d$binary_match_design$indicies_pairs
