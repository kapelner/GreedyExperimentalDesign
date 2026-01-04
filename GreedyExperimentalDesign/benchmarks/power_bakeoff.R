#!/usr/bin/env Rscript
# Power bake-off across design search methods on a fixed X (n=100, p=5).
#
# Usage: Rscript GreedyExperimentalDesign/benchmarks/power_bakeoff.R
# Adjust the parameters below to control effect sizes, sample size, and reps.

options(java.parameters = "-Xmx10g")

suppressPackageStartupMessages(library(GreedyExperimentalDesign))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(survival))
if (!requireNamespace("fastLogisticRegressionWrap", quietly = TRUE)) {
  devtools::install_github("kapelner/fastLogisticRegressionWrap")
}
suppressPackageStartupMessages(library(fastLogisticRegressionWrap))
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Package \"Rcpp\" is required for fast OLS.")
}
if (!requireNamespace("RcppEigen", quietly = TRUE)) {
  stop("Package \"RcppEigen\" is required for fast OLS.")
}
if (!exists("fast_ols_cpp", mode = "function")) {
  Rcpp::cppFunction(code = '
    #include <RcppEigen.h>
    Rcpp::List fast_ols_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
      Eigen::MatrixXd XtX = X.transpose() * X;
      Eigen::VectorXd Xty = X.transpose() * y;
      Eigen::VectorXd beta = XtX.ldlt().solve(Xty);
      return Rcpp::List::create(
        Rcpp::Named("b") = beta,
        Rcpp::Named("XtX") = XtX
      );
    }', depends = "RcppEigen", plugins = "cpp11")
}

SEED = 1984
N = 100
P = 2
NUM_CORES = 1
NUM_DESIGNS = 201
INIT_NUM_DESIGNS_BANK_MULTIPLE = 20
NUM_DESIGNS_BANK = NUM_DESIGNS * INIT_NUM_DESIGNS_BANK_MULTIPLE
NUM_OUTCOME_REPS = 2000
ALPHA = 0.05
OBJECTIVE = "mahal_dist"
MULTI_KERNEL_NAMES = c(
  "mahalanobis",
  "gaussian",
  "laplacian",
  "exponential",
  "inv_mult_quad",
  "poly_2",
  "poly_3"
)
USE_SAFE_INVERSE = TRUE

BETA_T_VALUES = c(0.5, 0)
BETA_T_CONT_SCALE = 0.5
BETA_T_BIN_SCALE = 1
BETA_T_SURV_SCALE = 0.5
BETA_X = c(0.3, -0.2, 0.15, 0.05, -0.1)

ALPHA_CONT = 0
SIGMA_Y = 1

ALPHA_BIN = -0.2

WEIBULL_SHAPE = 1.5
WEIBULL_LAMBDA = 0.05
CENSOR_RATE = 0

OUT_PATH = "power_bakeoff_results.csv"

set.seed(SEED)
X = generate_stdzied_design_matrix(n = N, p = P, covariate_gen = rnorm)

if (length(BETA_X) != P) {
  BETA_X = rep_len(BETA_X, P)
}
if (is.null(colnames(X))) {
  colnames(X) = paste0("X", seq_len(ncol(X)))
}

normalize_indicTs = function(indicTs, n) {
  if (is.null(indicTs)) {
    return(NULL)
  }
  indicTs = as.matrix(indicTs)
  if (ncol(indicTs) == n) {
    # ok
  } else if (nrow(indicTs) == n) {
    indicTs = t(indicTs)
  } else {
    stop("Unexpected indicTs dimensions: ", nrow(indicTs), "x", ncol(indicTs))
  }
  if (length(indicTs) > 0 && min(indicTs) < 0) {
    indicTs = (indicTs + 1) / 2
  }
  storage.mode(indicTs) = "integer"
  indicTs
}

take_designs = function(indicTs, num_designs) {
  if (is.null(indicTs)) {
    return(NULL)
  }
  if (nrow(indicTs) <= num_designs) {
    return(indicTs)
  }
  indicTs[seq_len(num_designs), , drop = FALSE]
}

take_best_by_obj_vals = function(indicTs, obj_vals, frac) {
  if (is.null(indicTs) || is.null(obj_vals) || length(obj_vals) == 0) {
    return(indicTs)
  }
  indicTs = as.matrix(indicTs)
  n_keep = min(length(obj_vals), nrow(indicTs))
  if (n_keep == 0) {
    return(indicTs)
  }
  n_take = max(1L, floor(n_keep * frac))
  order_idx = order(obj_vals[seq_len(n_keep)])
  indicTs[order_idx[seq_len(n_take)], , drop = FALSE]
}

safe_beta = function(expr) {
  tryCatch(expr, error = function(e) NA_real_)
}

simulate_continuous = function(w) {
  eta = as.numeric(X %*% BETA_X) + BETA_T_CONT * w + ALPHA_CONT
  eta + rnorm(N, sd = SIGMA_Y)
}

simulate_binary = function(w) {
  eta = as.numeric(X %*% BETA_X) + BETA_T_BIN * w + ALPHA_BIN
  prob = plogis(eta)
  rbinom(N, size = 1, prob = prob)
}

simulate_survival = function(w) {
  eta = as.numeric(X %*% BETA_X) + BETA_T_SURV * w
  u = runif(N)
  time = (-log(u) / (WEIBULL_LAMBDA * exp(eta)))^(1 / WEIBULL_SHAPE)
  status = rep(1L, N)
  if (CENSOR_RATE > 0) {
    censor_time = rexp(N, rate = CENSOR_RATE)
    status = as.integer(time <= censor_time)
    time = pmin(time, censor_time)
  }
  list(time = time, status = status)
}

beta_hat_diff = function(y, w) {
  mean(y[w == 1]) - mean(y[w == 0])
}

beta_hat_continuous_reg = function(y, w) {
  Xmm = cbind(1, w, X)
  fit = fast_ols_cpp(Xmm, y)
  if (is.null(fit$b) || length(fit$b) < 2) {
    return(NA_real_)
  }
  fit$b[[2]]
}

beta_hat_binary_reg = function(y, w) {
  Xmm = cbind(Intercept = 1, w = w, X)
  colnames(Xmm) = make.unique(colnames(Xmm))
  fit = fastLogisticRegressionWrap::fast_logistic_regression(
    Xmm = Xmm,
    ybin = y,
    do_inference_on_var = "none"
  )
  coef = fit$coefficients[["w"]]
  if (is.null(coef) || !is.finite(coef)) {
    return(NA_real_)
  }
  coef
}

beta_hat_survival = function(outcome, w) {
  mean(outcome$time[w == 1]) - mean(outcome$time[w == 0])
}

beta_hat_survival_reg = function(outcome, w) {
  n_obs = length(outcome$time)
  fit = suppressWarnings(survival::coxph.fit(
    x = cbind(w, X),
    y = cbind(outcome$time, outcome$status),
    strata = rep(1L, n_obs),
    offset = rep(0, n_obs),
    init = NULL,
    control = survival::coxph.control(),
    weights = rep(1, n_obs),
    method = "breslow",
    rownames = as.character(seq_len(n_obs)),
    resid = FALSE
  ))
  if (is.null(fit$coefficients) || length(fit$coefficients) < 1) {
    return(NA_real_)
  }
  fit$coefficients[[1]]
}

estimate_power = function(indicTs, method, outcome, simulate_fun, beta_hat_fun) {
  indicTs = take_designs(normalize_indicTs(indicTs, N), NUM_DESIGNS_BANK)
  if (is.null(indicTs) || nrow(indicTs) == 0) {
    return(data.frame(
      method = method,
      outcome = outcome,
      power = NA_real_,
      n_eval = 0,
      n_success = 0,
      n_total = NUM_OUTCOME_REPS,
      n_fail = NUM_OUTCOME_REPS,
      stringsAsFactors = FALSE
    ))
  }
  n_designs = nrow(indicTs)
  rejects = rep(NA, NUM_OUTCOME_REPS)
  n_fail = 0
  n_eval = 0
  for (j in seq_len(NUM_OUTCOME_REPS)) {
    obs_idx = sample.int(n_designs, 1)
    w_obs = as.numeric(indicTs[obs_idx, ])
    other_idx = setdiff(seq_len(n_designs), obs_idx)
    if (length(other_idx) < NUM_DESIGNS) {
      n_fail = n_fail + 1
      next
    }
    test_idx = sample(other_idx, NUM_DESIGNS, replace = FALSE)
    outcome_data = simulate_fun(w_obs)
    beta_obs = safe_beta(beta_hat_fun(outcome_data, w_obs))
    if (!is.finite(beta_obs)) {
      n_fail = n_fail + 1
      next
    }
    betas = vapply(test_idx, function(i) {
      safe_beta(beta_hat_fun(outcome_data, as.numeric(indicTs[i, ])))
    }, numeric(1))
    betas = betas[is.finite(betas)]
    if (length(betas) < 2) {
      n_fail = n_fail + 1
      next
    }
    qs = stats::quantile(betas, probs = c(ALPHA / 2, 1 - ALPHA / 2), type = 1, names = FALSE)
    rejects[j] = (beta_obs < qs[1]) || (beta_obs > qs[2])
    n_eval = n_eval + 1
  }
  n_success = sum(rejects, na.rm = TRUE)
  data.frame(
    method = method,
    outcome = outcome,
    power = if (n_eval > 0) n_success / n_eval else NA_real_,
    n_eval = n_eval,
    n_success = n_success,
    n_total = NUM_OUTCOME_REPS,
    n_fail = n_fail,
    stringsAsFactors = FALSE
  )
}

get_designs_greedy = function() {
  obj = initGreedyExperimentalDesignObject(
    X,
    max_designs = NUM_DESIGNS_BANK,
    objective = OBJECTIVE,
    num_cores = NUM_CORES,
    start = TRUE,
    wait = TRUE,
    seed = SEED,
    verbose = FALSE,
    use_safe_inverse = USE_SAFE_INVERSE
  )
  res = resultsGreedySearch(obj, max_vectors = NUM_DESIGNS_BANK, form = "one_zero")
  res$ending_indicTs
}

get_designs_multiple_kernel = function() {
  obj = initGreedyMultipleKernelExperimentalDesignObject(
    X,
    max_designs = NUM_DESIGNS_BANK,
    kernel_pre_num_designs = NUM_DESIGNS_BANK,
    kernel_names = MULTI_KERNEL_NAMES,
    objective = "added_pct_reduction",
    num_cores = NUM_CORES,
    start = TRUE,
    wait = TRUE,
    seed = SEED,
    diagnostics = FALSE,
    verbose = FALSE,
    use_safe_inverse = USE_SAFE_INVERSE
  )
  res = resultsMultipleKernelGreedySearch(obj, max_vectors = NUM_DESIGNS_BANK, form = "one_zero")
  res$ending_indicTs
}

get_designs_binary_match = function() {
  bms = computeBinaryMatchStructure(X, use_safe_inverse = USE_SAFE_INVERSE)
  obj = initBinaryMatchExperimentalDesignSearchObject(
    bms,
    max_designs = NUM_DESIGNS_BANK,
    num_cores = NUM_CORES,
    start = TRUE,
    wait = TRUE,
    seed = SEED,
    verbose = FALSE
  )
  resultsBinaryMatchSearch(obj, form = "one_zero")
}

get_designs_binary_match_then_greedy = function() {
  obj = initBinaryMatchFollowedByGreedyExperimentalDesignSearchObject(
    X,
    max_designs = NUM_DESIGNS_BANK,
    objective = OBJECTIVE,
    num_cores = NUM_CORES,
    start = TRUE,
    wait = TRUE,
    seed = SEED,
    diagnostics = FALSE,
    verbose = FALSE,
    use_safe_inverse = USE_SAFE_INVERSE
  )
  res = resultsBinaryMatchThenGreedySearch(obj, num_vectors = NUM_DESIGNS_BANK, form = "one_zero")
  res$indicTs
}

get_designs_binary_match_then_rerand = function() {
  obj = initBinaryMatchFollowedByRerandomizationDesignSearchObject(
    X,
    max_designs = NUM_DESIGNS_BANK * 10,
    objective = OBJECTIVE,
    obj_val_cutoff_to_include = Inf,
    num_cores = NUM_CORES,
    start = TRUE,
    wait = TRUE,
    seed = SEED,
    verbose = FALSE,
    use_safe_inverse = USE_SAFE_INVERSE
  )
  res = resultsBinaryMatchThenRerandomizationSearch(
    obj,
    num_vectors = NUM_DESIGNS_BANK * 10,
    compute_obj_vals = TRUE,
    form = "one_zero",
    use_safe_inverse = USE_SAFE_INVERSE
  )
  take_best_by_obj_vals(res$indicTs, res$obj_vals, 0.1)
}

get_designs_rerandomization = function() {
  obj = initRerandomizationExperimentalDesignObject(
    X,
    max_designs = NUM_DESIGNS_BANK * 10,
    objective = OBJECTIVE,
    obj_val_cutoff_to_include = Inf,
    num_cores = NUM_CORES,
    start = TRUE,
    wait = TRUE,
    seed = SEED,
    verbose = FALSE,
    use_safe_inverse = USE_SAFE_INVERSE
  )
  res = resultsRerandomizationSearch(obj, include_assignments = TRUE, form = "one_zero")
  take_best_by_obj_vals(res$ending_indicTs, res$obj_vals, 0.1)
}

get_designs_karp = function() {
  if (P != 1) {
    return(NULL)
  }
  obj = initKarpExperimentalDesignObject(
    X,
    balanced = TRUE,
    start = TRUE,
    wait = TRUE,
    verbose = FALSE
  )
  res = resultsKarpSearch(obj)
  matrix(res$indicT, nrow = 1)
}

get_designs_optimal = function() {
  if (N > 30) {
    return(NULL)
  }
  obj = initOptimalExperimentalDesignObject(
    X,
    objective = OBJECTIVE,
    num_cores = NUM_CORES,
    start = TRUE,
    wait = TRUE,
    verbose = FALSE,
    use_safe_inverse = USE_SAFE_INVERSE
  )
  res = resultsOptimalSearch(obj, num_vectors = NUM_DESIGNS_BANK, form = "one_zero")
  res$indicTs
}

get_designs_gurobi = function() {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    message("  Gurobi package not available; skipping.")
    return(NULL)
  }
  message("  Running Gurobi pool search.")
  obj = initGurobiNumericalOptimizationExperimentalDesignObject(
    X,
    r = NUM_DESIGNS_BANK,
    pool_solutions = 2 * NUM_DESIGNS_BANK,
    pool_search_mode = 2,
    pool_gap = 0.2,
    objective = OBJECTIVE,
    num_cores = NUM_CORES,
    initial_time_limit_sec = 20,
    restart_time_limit_sec = 5,
    allow_restarts = FALSE,
    verbose = FALSE,
    use_safe_inverse = USE_SAFE_INVERSE
  )
  res = resultsGurobiNumericalOptimizeSearch(obj)
  indicTs = res$indicTs
  if (is.null(dim(indicTs))) {
    indicTs = matrix(indicTs, nrow = 1)
  }
  gurobi_diag = list()
  if (!is.null(res$obj_vals) && length(res$obj_vals) > 0) {
    obj_summary = summary(res$obj_vals)
    gurobi_diag$obj_vals = c(
      min = obj_summary[[1]],
      median = obj_summary[[3]],
      mean = obj_summary[[4]],
      max = obj_summary[[6]]
    )
  } else {
    gurobi_diag$obj_vals = NULL
  }
  if (nrow(indicTs) > 1) {
    unique_indicTs = unique(indicTs)
    n_unique = nrow(unique_indicTs)
    gurobi_diag$n_unique = n_unique
    gurobi_diag$n_total = nrow(indicTs)
    if (n_unique > 1) {
      dists = stats::dist(unique_indicTs, method = "manhattan")
      if (length(dists) > 0) {
        dist_summary = summary(dists)
        gurobi_diag$hamming = c(
          min = dist_summary[[1]],
          median = dist_summary[[3]],
          mean = dist_summary[[4]],
          max = dist_summary[[6]]
        )
      }
    }
  } else {
    gurobi_diag$n_unique = 1L
    gurobi_diag$n_total = nrow(indicTs)
  }
  if (nrow(indicTs) < NUM_DESIGNS_BANK) {
    indicTs = indicTs[sample(seq_len(nrow(indicTs)), NUM_DESIGNS_BANK, replace = TRUE), , drop = FALSE]
  } else if (nrow(indicTs) > NUM_DESIGNS_BANK) {
    indicTs = indicTs[seq_len(NUM_DESIGNS_BANK), , drop = FALSE]
  }
  attr(indicTs, "gurobi_diag") = gurobi_diag
  indicTs
}

get_designs_complete_randomization = function() {
  complete_randomization_with_forced_balanced(
    n = N,
    r = NUM_DESIGNS_BANK,
    form = "one_zero",
    seed = SEED
  )
}

methods = list(
  complete_randomization = get_designs_complete_randomization,
  greedy = get_designs_greedy,
  multiple_kernel = get_designs_multiple_kernel,
  binary_match = get_designs_binary_match,
  binary_match_then_greedy = get_designs_binary_match_then_greedy,
  binary_match_then_rerandomization = get_designs_binary_match_then_rerand,
  rerandomization = get_designs_rerandomization,
  karp = get_designs_karp,
  optimal = get_designs_optimal,
  gurobi = get_designs_gurobi
)

method_abbrev = c(
  binary_match = "BM",
  binary_match_then_greedy = "BMG",
  binary_match_then_rerandomization = "BMR",
  complete_randomization = "BCRD",
  greedy = "G",
  gurobi = "GUR",
  multiple_kernel = "GMK",
  rerandomization = "R"
)
method_order = vapply(names(methods), function(m) {
  if (m %in% names(method_abbrev)) {
    method_abbrev[[m]]
  } else {
    m
  }
}, character(1))
pval_method_order = method_order
progress_base_cols = c(
  "outcome",
  "method",
  "beta_T",
  "regression",
  "power",
  "pval_vs_alpha",
  "n_eval",
  "n_fail",
  "n",
  "p",
  "num_designs",
  "seed"
)
progress_cols = c(progress_base_cols, pval_method_order)
if (file.exists(OUT_PATH)) {
  invisible(file.remove(OUT_PATH))
}

outcomes = list(
  continuous = list(simulate = simulate_continuous, beta_hat = list(nonreg = beta_hat_diff, reg = beta_hat_continuous_reg)),
  binary = list(simulate = simulate_binary, beta_hat = list(nonreg = beta_hat_diff, reg = beta_hat_binary_reg)),
  survival = list(simulate = simulate_survival, beta_hat = list(nonreg = beta_hat_survival, reg = beta_hat_survival_reg))
)

prop_z_pval = function(x1, n1, x2, n2) {
  if (is.na(x1) || is.na(n1) || is.na(x2) || is.na(n2) || n1 == 0 || n2 == 0) {
    return(NA_real_)
  }
  p_pool = (x1 + x2) / (n1 + n2)
  if (!is.finite(p_pool) || p_pool <= 0 || p_pool >= 1) {
    return(NA_real_)
  }
  se = sqrt(p_pool * (1 - p_pool) * (1 / n1 + 1 / n2))
  if (!is.finite(se) || se <= 0) {
    return(NA_real_)
  }
  z = (x1 / n1 - x2 / n2) / se
  2 * stats::pnorm(-abs(z))
}

prop_one_pval = function(x, n, p0) {
  if (is.na(x) || is.na(n) || n == 0 || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
    return(NA_real_)
  }
  p_hat = x / n
  se = sqrt(p0 * (1 - p0) / n)
  if (!is.finite(se) || se <= 0) {
    return(NA_real_)
  }
  z = (p_hat - p0) / se
  2 * stats::pnorm(-abs(z))
}

designs = list()
for (method_name in names(methods)) {
  cat("Generating designs for:", method_name, "\n")
  indicTs = tryCatch(methods[[method_name]](), error = function(e) {
    cat("  Skipping", method_name, "due to error:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(indicTs)) {
    cat("  No designs available for", method_name, "\n")
    next
  }
  if (identical(method_name, "gurobi")) {
    gurobi_diag = attr(indicTs, "gurobi_diag")
    if (!is.null(gurobi_diag)) {
      message("  Gurobi diagnostics:")
      if (!is.null(gurobi_diag$obj_vals)) {
        message(
          "    obj_vals: min ", format(gurobi_diag$obj_vals[["min"]], digits = 6),
          " median ", format(gurobi_diag$obj_vals[["median"]], digits = 6),
          " mean ", format(gurobi_diag$obj_vals[["mean"]], digits = 6),
          " max ", format(gurobi_diag$obj_vals[["max"]], digits = 6)
        )
      } else {
        message("    obj_vals: unavailable")
      }
      if (!is.null(gurobi_diag$n_unique) && !is.null(gurobi_diag$n_total)) {
        message(
          "    unique designs: ",
          gurobi_diag$n_unique,
          " of ",
          gurobi_diag$n_total
        )
      }
      if (!is.null(gurobi_diag$hamming)) {
        message(
          "    Hamming distances: min ", format(gurobi_diag$hamming[["min"]], digits = 6),
          " median ", format(gurobi_diag$hamming[["median"]], digits = 6),
          " mean ", format(gurobi_diag$hamming[["mean"]], digits = 6),
          " max ", format(gurobi_diag$hamming[["max"]], digits = 6)
        )
      }
    }
  }
  designs[[method_name]] = indicTs
}

write_progress_row = function(res_row) {
  row = as.data.frame(res_row, stringsAsFactors = FALSE)
  row$method = ifelse(
    row$method %in% names(method_abbrev),
    method_abbrev[row$method],
    row$method
  )
  row$power = round(row$power, 4)
  row$pval_vs_alpha = ifelse(
    row$beta_T == 0,
    formatC(prop_one_pval(row$n_success, row$n_eval, ALPHA), format = "f", digits = 4),
    ""
  )
  row$n = N
  row$p = P
  row$num_designs = NUM_DESIGNS
  row$seed = SEED
  for (col in pval_method_order) {
    row[[col]] = ""
  }
  row = row[, progress_cols, drop = FALSE]
  if (!file.exists(OUT_PATH)) {
    write.csv(row, OUT_PATH, row.names = FALSE)
  } else {
    write.table(row, OUT_PATH, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

results = list()
for (regression_flag in c(TRUE, FALSE)) {
  regression_label = if (regression_flag) "yes" else "no"
  for (beta_T_base in BETA_T_VALUES) {
    BETA_T_CONT = beta_T_base * BETA_T_CONT_SCALE
    BETA_T_BIN = beta_T_base * BETA_T_BIN_SCALE
    BETA_T_SURV = beta_T_base * BETA_T_SURV_SCALE
    for (outcome_name in names(outcomes)) {
      beta_hat_fun = if (regression_flag) {
        outcomes[[outcome_name]]$beta_hat$reg
      } else {
        outcomes[[outcome_name]]$beta_hat$nonreg
      }      
      group_results = list()
      for (method_name in names(designs)) {
        indicTs = designs[[method_name]]
        cat(
          "  Evaluating outcome:", outcome_name,
          "beta_T =", beta_T_base,
          "method =", method_name,
          "regression =", regression_label, "\n"
        )
        res = estimate_power(
          indicTs,
          method_name,
          outcome_name,
          outcomes[[outcome_name]]$simulate,
          beta_hat_fun
        )
        res$beta_T = ifelse(
          outcome_name == "continuous",
          BETA_T_CONT,
          ifelse(outcome_name == "survival", BETA_T_SURV, BETA_T_BIN)
        )
        res$regression = regression_label
        res$pval_vs_alpha = if (res$beta_T == 0) {
          prop_one_pval(res$n_success, res$n_eval, ALPHA)
        } else {
          NA_real_
        }
        cat(
          "    Power =", formatC(res$power, format = "f", digits = 4),
          ifelse(res$beta_T == 0,
            paste("pval_vs_alpha =", formatC(res$pval_vs_alpha, format = "f", digits = 4)),
            ""
          ),
          "n_eval =", res$n_eval,
          "n_fail =", res$n_fail, "\n"
        )
        write_progress_row(res)
        results[[length(results) + 1]] = res
        group_results[[length(group_results) + 1]] = res
      }
      if (beta_T_base != 0 && length(group_results) > 0) {
        group_df = do.call(rbind, group_results)
        group_df$method = ifelse(
          group_df$method %in% names(method_abbrev),
          method_abbrev[group_df$method],
          group_df$method
        )
        method_levels = method_order[method_order %in% group_df$method]
        pval_method_levels = method_levels
        if (length(pval_method_levels) >= 2) {
          group_df$method = factor(group_df$method, levels = method_levels)
          group_df = group_df[order(group_df$method), , drop = FALSE]
          pval_mat = matrix(NA_real_, nrow = length(pval_method_levels), ncol = length(pval_method_levels),
            dimnames = list(pval_method_levels, pval_method_levels))
          for (i in seq_len(nrow(group_df))) {
            method_i = as.character(group_df$method[i])
            for (j in seq_len(nrow(group_df))) {
              method_j = as.character(group_df$method[j])
              if (method_i %in% pval_method_levels &&
                method_j %in% pval_method_levels &&
                match(method_i, pval_method_levels) < match(method_j, pval_method_levels)) {
                pval_mat[method_i, method_j] = prop_z_pval(
                  group_df$n_success[i],
                  group_df$n_eval[i],
                  group_df$n_success[j],
                  group_df$n_eval[j]
                )
              }
            }
          }
          pval_mat = round(pval_mat, 2)
          cat(
            "P-value comparison matrix (beta_T =", beta_T_base,
            "outcome =", outcome_name,
            "regression =", regression_label, ")\n"
          )
          print(pval_mat, na.print = "")
        }
      }
    }
  }
}

results_df = do.call(rbind, results)
results_df$n = N
results_df$p = P
results_df$num_designs = NUM_DESIGNS
results_df$num_outcome_reps = NUM_OUTCOME_REPS
results_df$seed = SEED
results_df = results_df[, c("outcome", "method", "beta_T", "regression",
  setdiff(names(results_df), c("outcome", "method", "beta_T", "regression")))]
if ("method" %in% names(results_df)) {
  results_df$method = ifelse(
    results_df$method %in% names(method_abbrev),
    method_abbrev[results_df$method],
    results_df$method
  )
}
method_levels_all = method_order[method_order %in% results_df$method]
pval_method_levels = method_levels_all
results_df$method = factor(results_df$method, levels = method_levels_all)
results_df = results_df[order(results_df$outcome, results_df$method, results_df$regression, results_df$beta_T), , drop = FALSE]
results_df$method = as.character(results_df$method)
results_df$power = round(results_df$power, 4)
results_df$pval_vs_alpha = ifelse(
  results_df$beta_T == 0,
  formatC(
    round(mapply(prop_one_pval, results_df$n_success, results_df$n_eval, MoreArgs = list(p0 = ALPHA)), 4),
    format = "f",
    digits = 4
  ),
  ""
)

group_keys = unique(results_df[, c("outcome", "beta_T", "regression")])
pval_blocks = list()

for (g in seq_len(nrow(group_keys))) {
  outcome_name = group_keys$outcome[g]
  beta_T = group_keys$beta_T[g]
  regression_label = group_keys$regression[g]
  df_out = results_df[
    results_df$outcome == outcome_name &
      results_df$beta_T == beta_T &
      results_df$regression == regression_label,
    ,
    drop = FALSE
  ]
  df_out$method = factor(df_out$method, levels = method_levels_all)
  df_out = df_out[order(df_out$method), , drop = FALSE]
  df_out$method = as.character(df_out$method)
  df_out$pval_vs_alpha = if (beta_T == 0) {
    formatC(
      round(mapply(prop_one_pval, df_out$n_success, df_out$n_eval, MoreArgs = list(p0 = ALPHA)), 4),
      format = "f",
      digits = 4
    )
  } else {
    ""
  }
  if (beta_T == 0) {
    pval_df = as.data.frame(
      matrix(NA_real_, nrow = nrow(df_out), ncol = length(pval_method_levels),
        dimnames = list(NULL, pval_method_levels)),
      stringsAsFactors = FALSE
    )
    df_out = cbind(df_out, pval_df)
  } else {
    pval_mat = matrix(NA_real_, nrow = length(pval_method_levels), ncol = length(pval_method_levels),
      dimnames = list(pval_method_levels, pval_method_levels))
    for (i in seq_len(nrow(df_out))) {
      method_i = df_out$method[i]
      for (j in seq_len(nrow(df_out))) {
        method_j = df_out$method[j]
        if (method_i %in% pval_method_levels &&
          method_j %in% pval_method_levels &&
          match(method_i, pval_method_levels) < match(method_j, pval_method_levels)) {
          pval_mat[method_i, method_j] = prop_z_pval(
            df_out$n_success[i],
            df_out$n_eval[i],
            df_out$n_success[j],
            df_out$n_eval[j]
          )
        }
      }
    }
    pval_mat = round(pval_mat, 2)
    pval_df = as.data.frame(pval_mat, stringsAsFactors = FALSE)
    pval_rows = matrix(NA_real_, nrow = nrow(df_out), ncol = length(pval_method_levels),
      dimnames = list(NULL, pval_method_levels))
    for (i in seq_len(nrow(df_out))) {
      if (df_out$method[i] %in% pval_method_levels) {
        pval_rows[i, ] = pval_mat[df_out$method[i], ]
      }
    }
    df_out = cbind(df_out, as.data.frame(pval_rows, stringsAsFactors = FALSE))
  }
  pval_cols = setdiff(names(df_out), names(results_df))
  if (length(pval_cols) > 0) {
    df_out[pval_cols] = lapply(df_out[pval_cols], function(x) {
      ifelse(is.na(x), "", formatC(as.numeric(x), format = "f", digits = 2))
    })
  }
  pval_blocks[[paste(outcome_name, beta_T, regression_label, sep = "_")]] = df_out
}

results_df = do.call(rbind, pval_blocks)
results_df = results_df[, setdiff(names(results_df), c("n_success", "n_total", "num_outcome_reps")), drop = FALSE]
results_df = results_df[, c(
  "outcome",
  "method",
  "beta_T",
  "regression",
  "power",
  "pval_vs_alpha",
  setdiff(names(results_df), c("outcome", "method", "beta_T", "regression", "power", "pval_vs_alpha"))
), drop = FALSE]
names(results_df) = ifelse(
  names(results_df) %in% names(method_abbrev),
  method_abbrev[names(results_df)],
  names(results_df)
)
pval_cols = setdiff(
  names(results_df),
  c(
    "outcome",
    "method",
    "beta_T",
    "regression",
    "power",
    "pval_vs_alpha",
    "n_eval",
    "n_fail",
    "n",
    "p",
    "num_designs",
    "seed"
  )
)
if (length(pval_cols) > 0) {
  results_df[pval_cols] = lapply(results_df[pval_cols], function(x) {
    x = ifelse(is.na(x) | x == "NA", "", x)
    ifelse(x == "", "", formatC(as.numeric(x), format = "f", digits = 2))
  })
}

write.csv(results_df, OUT_PATH, row.names = FALSE)
print_df = results_df
pval_cols = setdiff(
  names(print_df),
  c(
    "outcome",
    "method",
    "beta_T",
    "regression",
    "power",
    "pval_vs_alpha",
    "n_eval",
    "n_fail",
    "n",
    "p",
    "num_designs",
    "seed"
  )
)
if (length(pval_cols) > 0) {
  print_df[pval_cols] = lapply(print_df[pval_cols], function(x) {
    x = ifelse(is.na(x) | x == "NA", "", x)
    ifelse(x == "", "", formatC(as.numeric(x), format = "f", digits = 2))
  })
}
print_df$power = formatC(print_df$power, format = "f", digits = 4)
print_df$pval_vs_alpha = ifelse(
  print_df$pval_vs_alpha == "" | is.na(print_df$pval_vs_alpha),
  "",
  formatC(as.numeric(print_df$pval_vs_alpha), format = "f", digits = 4)
)
print_df$beta_T = formatC(print_df$beta_T, format = "f", digits = 2)
int_cols = c("n_eval", "n_fail", "n", "p", "num_designs", "num_outcome_reps", "seed")
int_cols = intersect(int_cols, names(print_df))
print_df[int_cols] = lapply(print_df[int_cols], function(x) formatC(x, format = "d"))
print(print_df, row.names = FALSE, right = TRUE, na.print = "")
cat("Saved results to:", OUT_PATH, "\n")
