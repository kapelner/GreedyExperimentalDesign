options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, doParallel, tidyverse, magrittr, scales, data.table, r2r, checkmate, rlist)
#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.3/ GreedyExperimentalDesign

# alloc = optimize_asymmetric_treatment_assignment(
#   c_total_max = 100, c_treatment = 2, c_control = 1
# )
# 
# X = matrix(rnorm(alloc$n), ncol = 1)
# 
# ged = initGreedyExperimentalDesignObject(X, alloc$nT, objective = "abs_sum_diff", max_designs = 2, diagnostics = TRUE)
# res_ged = resultsGreedySearch(ged)
# res_ged


set.seed(1984)
ns = seq(from = 100, to = 1300, by = 300)
ps = c(1, 2, 5, 10)

X = matrix(rnorm(max(ns) * max(ps)), ncol = max(ps))
max_designs = 100
res = data.table(n = numeric(), p = numeric(), avg_obj_val = numeric())
for (i_n in 1 : length(ns)){
  n = ns[i_n]
  
  alloc = optimize_asymmetric_treatment_assignment(
    n = n, c_treatment = 2, c_control = 1
  )
  
  for (i_j in 1 : length(ps)){
    p = ps[i_j]
    cat("n", n, "p", p, "\n")
    
    ged = initGreedyExperimentalDesignObject(
            X[1 : alloc$n, 1 : p, drop = FALSE], 
            nT = alloc$nT, 
            objective = "abs_sum_diff", 
            max_designs = max_designs, 
            diagnostics = TRUE, 
            num_cores = 10,
            wait = TRUE
          )
    res_ged = resultsGreedySearch(ged)
    avg_obj_val = mean(res_ged$obj_vals)
    
    res = rbind(res, data.table(n = n, p = p, avg_obj_val = avg_obj_val))
  }
}
res[, p := as.factor(p)]
ggplot(res) +
  geom_point(aes(x = n, y = avg_obj_val, color = p)) + 
  scale_x_log10() + 
  scale_y_log10()
summary(lm(log(avg_obj_val) ~ log(n), res[p == 1]))
summary(lm(log(avg_obj_val) ~ log(n), res[p == 2]))
summary(lm(log(avg_obj_val) ~ log(n), res[p == 5]))
summary(lm(log(avg_obj_val) ~ log(n), res[p == 10]))

