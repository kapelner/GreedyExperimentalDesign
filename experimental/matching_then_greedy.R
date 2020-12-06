options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign

ns = c(32, 100, 1000, 2000)
ps = c(1, 2, 5, 10)
nR = 100
nX = 5
nC = 10
objective = "abs_sum_diff"
res = data.frame(imbalance = numeric(), method = character(), n = numeric(), p = numeric(), nx = numeric())

for (n in ns){
  for (p in ps){
    for (nx in 1 : nX){
      X = matrix(runif(n * p), nrow = n)
      # X = data.matrix(MASS::Pima.tr[1 : n, 1 : p])
      #stdize
      # X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
      bcrds = complete_randomization_with_forced_balanced(n, nR)
      bcrd_obj_vals = apply(bcrds, 1, function(w){compute_objective_val(X, w, objective = objective)})
      
      bmeds = binaryMatchExperimentalDesignSearch(X)
      bmeds_obj_vals = resultsBinaryMatchSearch(bmeds, num_vectors = nR, objective = objective)$obj_vals[1 : nR]
      # dim(bmeds_res$indicTs)
      
      bmfged = binaryMatchFollowedByGreedyExperimentalDesignSearch(X, objective = objective, max_designs = nR, wait = TRUE, num_cores = nC)
      # bmfged$binary_match_design
      # bmfged$greedy_design
      bmfged_obj_vals = resultsBinaryMatchThenGreedySearch(bmfged, compute_obj_vals = TRUE)$obj_vals
      # dim(bmfged_res$indicTs)
      
      ged = initGreedyExperimentalDesignObject(X, objective = objective, max_designs = nR, wait = TRUE, num_cores = nC)
      ged_res = resultsGreedySearch(ged, max_vectors = nR)
      # dim(ged_res$ending_indicTs)
      
      #now let's compare balance
      res = rbind(res,
                  data.frame(imbalance = bcrd_obj_vals, method = "bcrd", n = n, p = p, nx = nx),
                  data.frame(imbalance = ged_res$obj_vals, method = "greedy", n = n, p = p, nx = nx),
                  data.frame(imbalance = bmeds_obj_vals, method = "pair-matching", n = n, p = p, nx = nx),
                  data.frame(imbalance = bmfged_obj_vals, method = "pair-matching-then-greedy", n = n, p = p, nx = nx)
      )
      avg_summary = res %>% group_by(method, n, p) %>% summarize(med = median(log10(imbalance)))
      print(ggplot(res) +
              geom_density(aes(x = log10(imbalance), col = method, fill = method), position = "identity", alpha = 0.2) +
              facet_wrap(n ~ p, 
                         scales = "free_x",
                         ncol = length(ps),
                         labeller = function(labs){label_both(labs, multi_line = FALSE)}) +
              # scale_x_log10() + 
              xlab("log10 abs sum diff imbalance in x") +
              geom_vline(data = avg_summary, aes(xintercept = med, col = method), lwd = 2, alpha = 0.5, linetype = "dashed"))
      
      save(res, file = "matching_then_greedy_balance_sims_large_n_2000.RData")
    }
    
    
  }
}

res %>% group_by(method, n, p) %>% summarize(avg = mean(log10(imbalance)))

res = data.table(res)
res = res[imbalance != 0]
res[method == "pair-matching-then-greedy" & n == 500, imbalance]
sort(res[method == "pair-matching-then-greedy" & n == 500, imbalance])

res[, .(avg_log10_imbalance = mean(log10(imbalance))), by = c("n", "p", "method")]

t.test(res[n == 1000 & method == "greedy", log10(imbalance)],
       res[n == 1000 & method == "pair-matching-then-greedy", log10(imbalance)])


summary(lm(imbalance ~ log(n) * as.factor(p), res))
