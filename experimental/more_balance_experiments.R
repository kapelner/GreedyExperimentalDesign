options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table, viridis, RColorBrewer, ggsci, lmtest, sandwich)


ns = c(32)
p = 1
nR = 1000
nC = 10
objective = "abs_sum_diff"

# res = data.table(imbalance = numeric(), n = numeric(), design = character())
for (n in ns){
  # X = matrix(runif(n * p, -sqrt(3), sqrt(3)), nrow = n)
  X = matrix(rnorm(n * p), nrow = n)
  
  bmfged = binaryMatchFollowedByGreedyExperimentalDesignSearch(X, objective = objective, max_designs = nR, wait = TRUE, num_cores = nC)
  bmfged_res = resultsBinaryMatchThenGreedySearch(bmfged, compute_obj_vals = TRUE)
  res = rbind(res, data.table(imbalance = bmfged_res$obj_vals, n = n, design = "matching_then_greedy"))
  
  ged = initGreedyExperimentalDesignObject(X, objective = objective, wait = TRUE, max_designs = nR, num_cores = nC)
  ged_res = resultsGreedySearch(ged)
  res = rbind(res, data.table(imbalance = ged_res$obj_vals, n = n, design = "greedy"))
}

avg_summary = res[, .(avg_imbalance = mean(imbalance)), by = c("n", "design")]
avg_summary

mod = lm(imbalance ~ design * n, res)
summary(mod)

ggplot(res) + 
  aes(x = imbalance) + 
  geom_density() + 
  facet_grid(design ~ n) + 
  scale_x_log10() + 
  geom_vline(data = avg_summary, aes(xintercept = avg_imbalance), col = "green")

