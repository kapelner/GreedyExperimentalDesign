options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, doParallel, tidyverse, magrittr, data.table, viridis, RColorBrewer, ggsci, lmtest, sandwich)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign

cl = makeCluster(nC)
registerDoParallel(cl)
clusterExport(cl, list(
  "exp_settings", "all_betas_and_correlations", "Xalls", "betaT", "sigma_e", "nR"
), envir = environment())

res = data.table(n = factor(), p = factor(), nx = factor(), neps = factor(), nsim = factor(), model = factor(), design = factor(), squared_error = numeric())


res = foreach(n_setting = 1 : nrow(exp_settings), .inorder = FALSE, .combine = rbind, .packages = "GreedyExperimentalDesign") %dopar% {
  inner_res = data.frame()
  
  #get this run's settings
  p = exp_settings$p[n_setting]
  nx = exp_settings$nx[n_setting]
  neps = exp_settings$neps[n_setting]
  
  #get this run's data
  X = Xalls[[nx]][, 1 : p]
  epsilons = epsilon_alls[[neps]]
  
  ##ADD BCRD HERE
  cat("bcrds...")
  bcrds = complete_randomization_with_forced_balanced(n, nR)
  
  cat("bmeds...")
  bmeds = binaryMatchExperimentalDesignSearch(X)
  bmeds_res = resultsBinaryMatchSearch(bmeds, num_vectors = min(nR, 2^(n/2)), objective = objective)
  # dim(bmeds_res$indicTs)
  
  cat("bmfged...")
  bmfged = binaryMatchFollowedByGreedyExperimentalDesignSearch(X, objective = objective, max_designs = nR, wait = TRUE)
  # bmfged$binary_match_design
  # bmfged$greedy_design
  bmfged_res = resultsBinaryMatchThenGreedySearch(bmfged)
  # dim(bmfged_res$indicTs)
  
  #rerand
  cat("rerands...")
  bcrds_rerand = complete_randomization_with_forced_balanced(n, nR * 1 / rerand_threshold)
  bcrds_rerand_balances = apply(bcrds_rerand, 1, function(w){compute_objective_val(X, w)})
  rerands = bcrds_rerand[order(bcrds_rerand_balances)[1 : nR], ]
  
  #matching then rerand
  cat("matching_then_rerands...")
  bmeds_rerand = resultsBinaryMatchSearch(bmeds, num_vectors = nR * 1 / rerand_threshold, objective = objective)$indicTs
  bmeds_rerand_balances = apply(bmeds_rerand, 1, function(w){compute_objective_val(X, w)})
  matching_then_rerands = bmeds_rerand[order(bmeds_rerand_balances)[1 : nR], ]
  
  cat("greedy...")
  ged = initGreedyExperimentalDesignObject(X, objective = objective, wait = TRUE)
  ged_res = resultsGreedySearch(ged, max_vectors = nR)
  # dim(ged_res$ending_indicTs)
  
  
  for (nsim in 1 : nR){
    for (model_name in names(all_betas_and_correlations)){
      params = all_betas_and_correlations[[model_name]]
      y = params$beta_x1 * X[, 1] + 
        params$beta_x2 * X[, 2] + 
        params$beta_x1sq * X[, 1]^2 + 
        params$beta_x2sq * X[, 2]^2 + 
        params$beta_x1_x2 * X[, 1] * X[, 2] + 
        epsilons
      
      adjust_y_for_w = function(y, w){
        y + betaT * w + params$beta_x1_T * w * X[, 1] + params$beta_x2_T * w * X[, 2]
      }
      
      #bcrd's
      w = bcrds[nsim, ]
      y_bcrd = adjust_y_for_w(y, w)
      betaT_hat_cbrd = mean(y_bcrd[w == 1]) - mean(y_bcrd[w == 0])
      #rerand
      w = rerands[nsim, ]
      y_rerand = adjust_y_for_w(y, w)
      betaT_hat_rerand = mean(y_rerand[w == 1]) - mean(y_rerand[w == 0])        
      #matching
      w = bmeds_res$indicTs[nsim,]
      y_matching = adjust_y_for_w(y, w)
      betaT_hat_matching = mean(y_matching[w == 1]) - mean(y_matching[w == 0])
      #matching then greedy
      w = bmfged_res$indicTs[nsim,]
      y_matching_then_greedy = adjust_y_for_w(y, w)
      betaT_hat_matching_then_greedy = mean(y_matching_then_greedy[w == 1]) - mean(y_matching_then_greedy[w == 0])
      #matching then rerand
      w = matching_then_rerands[nsim,]
      y_matching_then_rerand = adjust_y_for_w(y, w)
      betaT_hat_matching_then_rerand = mean(y_matching_then_rerand[w == 1]) - mean(y_matching_then_rerand[w == 0])
      
      #just greedy
      w = ged_res$ending_indicTs[nsim,]
      y_greedy = adjust_y_for_w(y, w)
      betaT_greedy = mean(y_greedy[w == 1]) - mean(y_greedy[w == 0])
      
      inner_res = rbind(
        inner_res,
        data.frame(design = "bcrd",                 squared_error = (betaT_hat_cbrd - betaT)^2,                 n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
        data.frame(design = "rerand",               squared_error = (betaT_hat_rerand - betaT)^2,               n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
        data.frame(design = "matching",             squared_error = (betaT_hat_matching - betaT)^2,             n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
        data.frame(design = "matching_then_greedy", squared_error = (betaT_hat_matching_then_greedy - betaT)^2, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
        data.frame(design = "matching_then_rerand", squared_error = (betaT_hat_matching_then_rerand - betaT)^2, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
        data.frame(design = "greedy",               squared_error = (betaT_greedy - betaT)^2,                   n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name)
      )
    }
  }
  
  save(inner_res, file = paste0(filename, "_", n_setting, ".RData"))
  inner_res
}
stopCluster(cl)
save(res, file = paste0(filename, ".RData"))



