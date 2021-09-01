

cl = makeCluster(nC)
registerDoParallel(cl)
# registerDoSEQ() #### for debugging only
clusterExport(cl, list(
  "exp_settings", "all_betas_and_correlations", "Xalls", "betaT", "nR"
), envir = environment())

res = data.table(n = factor(), p = factor(), nx = factor(), neps = factor(), nsim = factor(), model = factor(), design = factor(), estimator = factor(), squared_error = numeric())


res = foreach(n_setting = 1 : nrow(exp_settings), .inorder = FALSE, .combine = rbind, .packages = "GreedyExperimentalDesign") %dopar% {
# for (n_setting in 1 : nrow(exp_settings)){
  cat("n_setting: ", n_setting, "/", nrow(exp_settings), " ", paste(exp_settings[n_setting, ], collapse = "-"), "\n")
  
  inner_res = data.frame()
  
  #get this run's settings
  covariate_distribution = exp_settings$covariate_distribution[n_setting]
  p = exp_settings$p[n_setting]
  nx = exp_settings$nx[n_setting]
  neps = exp_settings$neps[n_setting]
  
  #get this run's data
  X = Xalls[[covariate_distribution]][[nx]][, 1 : p, drop = FALSE]
  Xscale = scale(X) #this is the way it goes into the design
  epsilons = epsilon_alls[[neps]]
  
  cat("  bcrds...")
  bcrds = complete_randomization_with_forced_balanced(n, nR)
  cat("done\n")
  #rerand
  cat("  rerands...")
  bcrds_rerand = complete_randomization_with_forced_balanced(n, nR * 1 / rerand_threshold)
  bcrds_rerand_balances = apply(bcrds_rerand, 1, function(w){compute_objective_val(Xscale, w)})
  rerands = bcrds_rerand[order(bcrds_rerand_balances)[1 : nR], ]
  cat("done\n")
  
  cat("  bmeds...")
  bmeds = initBinaryMatchExperimentalDesignSearch(Xscale)
  bmeds_res = resultsBinaryMatchSearch(bmeds, num_vectors = min(nR, 2^(n/2)))
  cat("done\n")
  # dim(bmeds_res$indicTs)
  
  #matching then rerand
  # cat("matching_then_rerands...")
  # bmeds_rerand = resultsBinaryMatchSearch(bmeds, num_vectors = nR * 1 / rerand_threshold, objective = objective)$indicTs
  # bmeds_rerand_balances = apply(bmeds_rerand, 1, function(w){compute_objective_val(Xscale, w)})
  # matching_then_rerands = bmeds_rerand[order(bmeds_rerand_balances)[1 : nR], ]
  
  Ks = list()
  for (kernel_name in kernel_names){
    if (kernel_name == "mahalanobis"){
      Ks[[kernel_name]] = Xscale %*% solve(var(Xscale)) %*% t(Xscale)
    } else {
      Ks[[kernel_name]] = matrix(NA, n, n)
      for (i in 1 : n){
        for (j in 1 : n){
          xi = Xscale[i, , drop = FALSE]
          xj = Xscale[j, , drop = FALSE]
          #https://people.eecs.berkeley.edu/~jordan/kernels/0521813972c01_p03-24.pdf
          if (kernel_name == "poly_2"){
            Ks[[kernel_name]][i, j] = (1 + xi %*% t(xj) / 2)^2
          } else if (kernel_name == "exponential"){
            Ks[[kernel_name]][i, j] = exp(xi %*% t(xj))
          } else if (kernel_name == "gaussian"){
            Ks[[kernel_name]][i, j] = exp(-sum((xi - xj)^2))
          } else if (kernel_name == "laplacian"){
            Ks[[kernel_name]][i, j] = exp(-sqrt(sum((xi - xj)^2)))
          } else if (kernel_name == "inv_mult_quad"){
            #http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/#inverse_multiquadric
            Ks[[kernel_name]][i, j] = 1 / sqrt(sum((xi - xj)^2) + 1)
          }  
        }
      }       
    }
  }
    
  bmfged_ress = list()
  ged_ress = list()
  for (kernel_objective_function in kernel_objective_functions){
    if (kernel_objective_function %in% kernel_names){
      K = Ks[[kernel_objective_function]]
    } else if (kernel_objective_function == "mahalanobis_plus_gaussian"){
      K = Ks[["mahalanobis"]] + Ks[["gaussian"]]
    } else if (kernel_objective_function == "mahalanobis_plus_exponential"){
      K = Ks[["mahalanobis"]] + Ks[["exponential"]]
    } else if (kernel_objective_function == "mahalanobis_plus_laplacian"){
      K = Ks[["mahalanobis"]] + Ks[["laplacian"]]
    } else if (kernel_objective_function == "mahalanobis_plus_inv_mult_quad"){
      K = Ks[["mahalanobis"]] + Ks[["inv_mult_quad"]]
    } else if (kernel_objective_function == "mahalanobis_plus_gaussian_plus_inv_mult_quad"){
      K = Ks[["mahalanobis"]] + Ks[["gaussian"]] + Ks[["inv_mult_quad"]]
    }
    
    cat("  greedy ", kernel_objective_function, "...")
    ged = initGreedyExperimentalDesignObject(Xscale, Kgram = K, objective = "kernel", max_designs = nR, wait = TRUE)
    ged_ress[[kernel_objective_function]] = resultsGreedySearch(ged, max_vectors = nR)
    cat("done\n")

    cat("  bmfged ", kernel_objective_function, "...")
    bmfged = initBinaryMatchFollowedByGreedyExperimentalDesignSearch(Xscale, Kgram = K, objective = "kernel", max_designs = nR, wait = TRUE)
    bmfged_ress[[kernel_objective_function]] = resultsBinaryMatchThenGreedySearch(bmfged)
    cat("done\n")
    
    # dim(ged_res$ending_indicTs)
    # bmfged$binary_match_design
    # bmfged$greedy_design
    # dim(bmfged_res$indicTs)  
  }
  
  for (nsim in 1 : nR){
    if (nsim %% 100 == 0 || nsim == 1){
      cat("    recording results for ", nsim, "/", nR, "\n")
    }
    
    for (model_name in names(all_betas_and_correlations)){
      params = all_betas_and_correlations[[model_name]]
      y = 
        ifelse(is.null(params$beta_x1), 0, params$beta_x1) * X[, 1] + 
        ifelse(is.null(params$beta_x2), 0, params$beta_x2) * X[, 2] + 
        ifelse(is.null(params$beta_x1sq), 0, params$beta_x1sq) * X[, 1]^2 + 
        ifelse(is.null(params$beta_x2sq), 0, params$beta_x2sq) * X[, 2]^2 + 
        ifelse(is.null(params$beta_x1_x2), 0, params$beta_x1_x2) * X[, 1] * X[, 2] + 
        sin(ifelse(is.null(params$beta_x1_sin), 0, params$beta_x1_sin) * X[, 1]) +
        sin(ifelse(is.null(params$beta_x_2_sin), 0, params$beta_x_2_sin) * X[, 2]) +
        sin(ifelse(is.null(params$beta_x1_x_2_sin), 0, params$beta_x1_x_2_sin) * X[, 1] * X[, 2]) +
        ifelse(is.null(params$beta_x1_exp), 0, params$beta_x1_exp) * exp(X[, 1]) + 
        ifelse(is.null(params$beta_x2_exp), 0, params$beta_x2_exp) * exp(X[, 2]) + 
        ifelse(is.null(params$beta_x1_x_2_exp), 0, params$beta_x1_x_2_exp) * exp(X[, 1] * X[, 2]) + 
        epsilons
      
      add_effect_of_w_to_y = function(y, w){
        y + betaT * w + 
          ifelse(is.null(params$beta_x1_T), 0, params$beta_x1_T) * w * X[, 1] + 
          ifelse(is.null(params$beta_x2_T), 0, params$beta_x2_T) * w * X[, 2]
      }
      
      #bcrd's
      w = bcrds[nsim, ]
      y_bcrd = add_effect_of_w_to_y(y, w)
      betaT_hat_bcrd = mean(y_bcrd[w == 1]) - mean(y_bcrd[w == 0])
      betaT_hat_bcrd_ols = lm.fit(cbind(1, w, X), y_bcrd)$coefficients[2]
      
      #rerand
      w = rerands[nsim, ]
      y_rerand = add_effect_of_w_to_y(y, w)
      betaT_hat_rerand = mean(y_rerand[w == 1]) - mean(y_rerand[w == 0])
      betaT_hat_rerand_ols = lm.fit(cbind(1, w, X), y_rerand)$coefficients[2]
      
      #matching
      w = bmeds_res$indicTs[nsim,]
      y_matching = add_effect_of_w_to_y(y, w)
      betaT_hat_matching = mean(y_matching[w == 1]) - mean(y_matching[w == 0])
      pairs = bmeds$indicies_pairs
      iCs = which(w == 0)
      pairs_to_flip = intersect(iCs, pairs[, 1])
      flip_encoding = rep(1, n / 2)
      flip_encoding[match(pairs_to_flip, pairs[, 1])] = -1
      XDelta = (X[pairs[, 1], ] - X[pairs[, 2], ]) * flip_encoding
      yDelta = (y_matching[pairs[, 1]] - y_matching[pairs[, 2]]) * flip_encoding
      betaT_hat_matching_ols = lm.fit(cbind(1, XDelta), yDelta)$coefficients[1]
      
      #matching then rerand
      # w = matching_then_rerands[nsim,]
      # y_matching_then_rerand = add_effect_of_w_to_y(y, w)
      # betaT_hat_matching_then_rerand = mean(y_matching_then_rerand[w == 1]) - mean(y_matching_then_rerand[w == 0])
      
      inner_res = rbind(
        inner_res,
        data.frame(design = "B",   kernel = "",    squared_error = (betaT_hat_bcrd - betaT)^2,         estimator = "naive", covariate_distribution = covariate_distribution, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
        data.frame(design = "R",   kernel = "",    squared_error = (betaT_hat_rerand - betaT)^2,       estimator = "naive", covariate_distribution = covariate_distribution, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
        data.frame(design = "M",   kernel = "",    squared_error = (betaT_hat_matching - betaT)^2,     estimator = "naive", covariate_distribution = covariate_distribution, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
        data.frame(design = "B",   kernel = "",    squared_error = (betaT_hat_bcrd_ols - betaT)^2,     estimator = "ols",   covariate_distribution = covariate_distribution, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
        data.frame(design = "R",   kernel = "",    squared_error = (betaT_hat_rerand_ols - betaT)^2,   estimator = "ols",   covariate_distribution = covariate_distribution, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
        data.frame(design = "M",   kernel = "",    squared_error = (betaT_hat_matching_ols - betaT)^2, estimator = "ols",   covariate_distribution = covariate_distribution, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name)
      )
      
      
      for (kernel_objective_function in kernel_objective_functions){
        #just greedy
        w = ged_ress[[kernel_objective_function]]$ending_indicTs[nsim, ]
        y_greedy = add_effect_of_w_to_y(y, w)
        betaT_greedy = mean(y_greedy[w == 1]) - mean(y_greedy[w == 0])
        betaT_greedy_ols = lm.fit(cbind(1, w, X), y_greedy)$coefficients[2]  
        
        #matching then greedy
        w = bmfged_ress[[kernel_objective_function]]$indicTs[nsim, ]
        y_matching_then_greedy = add_effect_of_w_to_y(y, w)
        betaT_hat_matching_then_greedy = mean(y_matching_then_greedy[w == 1]) - mean(y_matching_then_greedy[w == 0])
        pairs = bmfged$binary_match_design$indicies_pairs
        #we need the first col of pairs to be all T indicies and the second col of pairs to be all C indicies
        #to do so, we flip the indicies that correspond to C's in the first col
        iCs = which(w == 0)
        pairs_to_flip = intersect(iCs, pairs[, 1])
        flip_encoding = rep(1, n / 2)
        flip_encoding[match(pairs_to_flip, pairs[, 1])] = -1
        XDelta = (X[pairs[, 1], ] - X[pairs[, 2], ]) * flip_encoding
        yDelta = (y_matching_then_greedy[pairs[, 1]] - y_matching_then_greedy[pairs[, 2]]) * flip_encoding
        betaT_hat_matching_then_greedy_ols = lm.fit(cbind(1, XDelta), yDelta)$coefficients[1]
        
        inner_res = rbind(
          inner_res,
          data.frame(design = "G",    kernel = kernel_objective_function,  squared_error = (betaT_greedy - betaT)^2,                       estimator = "naive", covariate_distribution = covariate_distribution, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
          data.frame(design = "MG",   kernel = kernel_objective_function,  squared_error = (betaT_hat_matching_then_greedy - betaT)^2,     estimator = "naive", covariate_distribution = covariate_distribution, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
          data.frame(design = "G",    kernel = kernel_objective_function,  squared_error = (betaT_greedy_ols - betaT)^2,                   estimator = "ols",   covariate_distribution = covariate_distribution, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name),
          data.frame(design = "MG",   kernel = kernel_objective_function,  squared_error = (betaT_hat_matching_then_greedy_ols - betaT)^2, estimator = "ols",   covariate_distribution = covariate_distribution, n = n, p = p, nx = nx, neps = neps, nsim = nsim, model = model_name)
        )
      }
    }
  }
  
  save(inner_res, file = paste0(filename, "_", n_setting, ".RData"))
  inner_res
}
stopCluster(cl)
save(res, file = paste0(filename, ".RData"))



