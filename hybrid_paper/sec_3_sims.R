options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, doParallel, tidyverse, magrittr, data.table, viridis, RColorBrewer, ggsci, lmtest, sandwich)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign

source("sim_params.R")
nC = 6
cl = makeCluster(nC)
registerDoParallel(cl)


res = data.table(n = factor(), p = factor(), nx = factor(), neps = factor(), nsim = factor(), model = factor(), design = factor(), squared_error = numeric())
set.seed(1)

Xalls = list()
for (nx in 1 : nX){
  Xalls[[nx]] = matrix(runif(n * max(ps), -sqrt(3), sqrt(3)), nrow = n)
  # Xtemp = matrix(rnorm(n * p, 0, 1), nrow = n)
  # Xtemp = data.matrix(MASS::Pima.tr[1 : n, 1 : p])
  # Xalls[[nx]] = apply(Xtemp, 2, function(xj){(xj - mean(xj)) / sd(xj)}) #stdize
  # rm(Xtemp)
}

epsilon_alls = list()
for (neps in 1 : nEPS){
  epsilon_alls[[neps]] = rnorm(n, 0, sigma_e)
}


#create experiment setting matrix
# exp_settings = data.frame(matrix(NA, nrow = length(ps) * nX * nEPS, ncol = 3))
# names(exp_settings) = c("p", "nx", "neps")
# exp_settings$neps = rep(1 : nEPS, times = length(ps) * nX)
# exp_settings$nx   = rep(1 : nX,   each  = nEPS)
# exp_settings$p    = rep(ps,       each  = nEPS * nX)

exp_settings = data.frame(matrix(NA, nrow = 150, ncol = 3))
names(exp_settings) = c("p", "nx", "neps")
exp_settings$neps = rep(1, 150)
exp_settings$nx =   rep(1, 150)
exp_settings$p =    rep(ps, each = 50)

clusterExport(cl, list(
    "exp_settings", "all_betas_and_correlations", "Xalls", "betaT", "sigma_e", "nR", "nR0"
  ), envir = environment())

res = foreach(n_setting = 1 : nrow(exp_settings), .inorder = FALSE, .combine = rbind, .packages = "GreedyExperimentalDesign") %dopar% {
  inner_res = data.frame()

  #get this run's settings
  p = exp_settings$p[n_setting]
  nx = exp_settings$nx[n_setting]
  neps = exp_settings$neps[n_setting]

  #get this run's data
  X = Xalls[[nx]][, 1 : p]
  epsilons = epsilon_alls[[neps]]
  nr = ifelse(nx == 1 & neps == 1, nR0, nR)

  ##ADD BCRD HERE
  cat("bcrds...")
  bcrds = complete_randomization_with_forced_balanced(n, nr)
  
  cat("bmeds...")
  bmeds = binaryMatchExperimentalDesignSearch(X)
  bmeds_res = resultsBinaryMatchSearch(bmeds, num_vectors = min(nr, 2^(n/2)), objective = objective)
  # dim(bmeds_res$indicTs)
  
  cat("bmfged...")
  bmfged = binaryMatchFollowedByGreedyExperimentalDesignSearch(X, objective = objective, max_designs = nr, wait = TRUE)
  # bmfged$binary_match_design
  # bmfged$greedy_design
  bmfged_res = resultsBinaryMatchThenGreedySearch(bmfged)
  # dim(bmfged_res$indicTs)
  
  #rerand
  cat("rerands...")
  bcrds_rerand = complete_randomization_with_forced_balanced(n, nr * 1 / rerand_threshold)
  bcrds_rerand_balances = apply(bcrds_rerand, 1, function(w){compute_objective_val(X, w)})
  rerands = bcrds_rerand[order(bcrds_rerand_balances)[1 : nr], ]
  
  #matching then rerand
  cat("matching_then_rerands...")
  bmeds_rerand = resultsBinaryMatchSearch(bmeds, num_vectors = nr * 1 / rerand_threshold, objective = objective)$indicTs
  bmeds_rerand_balances = apply(bmeds_rerand, 1, function(w){compute_objective_val(X, w)})
  matching_then_rerands = bmeds_rerand[order(bmeds_rerand_balances)[1 : nr], ]
  
  cat("greedy...")
  ged = initGreedyExperimentalDesignObject(X, objective = objective, wait = TRUE)
  ged_res = resultsGreedySearch(ged, max_vectors = nr)
  # dim(ged_res$ending_indicTs)
  
  
  for (nsim in 1 : nr){
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
  
  save(inner_res, file = paste0("matching_then_greedy_squared_error_sims_unif_", n_setting, ".RData"))
  inner_res
}
stopCluster(cl)
save(res, file = "matching_then_greedy_squared_error_sims_unif_all.RData")



# another idea: get ideal error which is...
Nsim_ideal = 10000000
bcrds = complete_randomization_with_forced_balanced(n, Nsim_ideal, form = "pos_one_min_one")
betaThat_ideals = array(NA, Nsim_ideal)
for (nsim_ideal in 1 : Nsim_ideal){
  w = bcrds[nsim_ideal, ]
  y = w * betaT + epsilon_alls[[1]]
  betaThat_ideals[nsim_ideal] = sum(y * w) / n
}
ideal_sq_err = mean((betaThat_ideals - betaT)^2)
log10(ideal_sq_err)
