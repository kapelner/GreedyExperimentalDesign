options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, doParallel, tidyverse, magrittr, data.table, r2r, checkmate, rlist)

nC = 1
n = 48
x = 1 : n
nR = 1e7
p = 20
PM_match_structure = computeBinaryMatchStructure(as.matrix(data.frame(x = x)), mahal_match = TRUE)
PM_match_structure = computeBinaryMatchStructure(matrix(rnorm(n * p), ncol = p, nrow = n), mahal_match = TRUE)
PM_match_structure
PM_des = initBinaryMatchExperimentalDesignSearch(PM_match_structure, num_cores = nC, wait = TRUE, max_designs = min(nR, 2^(n/2)))
PM_des
w_PM = resultsBinaryMatchSearch(PM_des, form = "pos_one_min_one")
# w_PM

dim(w_PM)
nrow(w_PM)
nrow(unique(data.table(w_PM)))
all(rowSums(w_PM) == 0)
pair_sums = array(NA, n / 2)
for (m in 1 : nrow(PM_match_structure$indicies_pairs)){
  pair = PM_match_structure$indicies_pairs[m, ]
  pair_sums[m] = all(rowSums(w_PM[, pair]) == 0)
}
all(pair_sums)
