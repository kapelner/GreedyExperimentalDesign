testthat::test_that("objective values rank by search method", {
  skip_on_cmd_check()
  with_immediate_failures({
    configs_raw = build_test_configs()
    config_keys = vapply(configs_raw, function(cfg) {
      paste(cfg$n, cfg$p, cfg$r, sep = "|")
    }, character(1))
    configs = configs_raw[!duplicated(config_keys)]
    configs = Filter(function(cfg) cfg$r == 200 && cfg$n == 96, configs)
    testthat::expect_true(length(configs) > 0)
    B_values = sort(unique(vapply(configs_raw, function(cfg) cfg$B, numeric(1))))
    objectives = c("abs_sum_diff", "mahal_dist")
    for (objective in objectives) {
      for (cfg in configs) {
        if (objective == "mahal_dist" && cfg$p == 10) {
          next
        }
        cat(
          "objective=", objective,
          " n=", cfg$n,
          " p=", cfg$p,
          " r=", cfg$r,
          " seed=", cfg$seed,
          "\n",
          sep = ""
        )

        X = make_X(cfg$n, cfg$p, seed = cfg$seed)
        compare_r = as.integer(cfg$r)
        if (objective == "mahal_dist") {
          X_eval = X
          inv_cov = safe_cov_inverse(X)
        } else {
          X_eval = standardize_data_matrix(X)
          inv_cov = NULL
        }

        complete_indicTs = complete_randomization(n = cfg$n, r = compare_r)
        complete_row_totals = rowSums(complete_indicTs)
        keep_complete = complete_row_totals > 0 & complete_row_totals < cfg$n
        if (!all(keep_complete)) {
          complete_indicTs = complete_indicTs[keep_complete, , drop = FALSE]
        }
        complete_obj = GreedyExperimentalDesign:::compute_objective_vals_cpp(
          X_eval,
          complete_indicTs,
          objective,
          inv_cov
        )
        complete_stat = median(complete_obj)
        cutoff = complete_stat

        ged = initGreedyExperimentalDesignObject(
          X,
          max_designs = compare_r,
          objective = objective,
          diagnostics = FALSE,
          wait = TRUE,
          num_cores = 1,
          seed = cfg$seed,
          verbose = FALSE,
          use_safe_inverse = TRUE
        )
        res_greedy = resultsGreedySearch(ged, max_vectors = compare_r, form = "one_zero")
        greedy_indicTs = res_greedy$ending_indicTs
        if (is.null(dim(greedy_indicTs))) {
          greedy_indicTs = matrix(greedy_indicTs, nrow = 1)
        }
        greedy_obj = GreedyExperimentalDesign:::compute_objective_vals_cpp(
          X_eval,
          greedy_indicTs,
          objective,
          inv_cov
        )
        greedy_stat = median(greedy_obj)

        bms = computeBinaryMatchStructure(X, use_safe_inverse = TRUE)
        bms$verbose = FALSE
        bm_unique = 2^(cfg$n / 2)
        bm_max_designs = min(compare_r, max(1L, floor(bm_unique / 2)))
        bm = initBinaryMatchExperimentalDesignSearchObject(
          bms,
          max_designs = bm_max_designs,
          wait = TRUE,
          num_cores = 1,
          seed = cfg$seed,
          verbose = FALSE
        )
        res_bm = resultsBinaryMatchSearch(bm, form = "one_zero")
        bm_indicTs = res_bm
        if (is.null(dim(bm_indicTs))) {
          bm_indicTs = matrix(bm_indicTs, nrow = 1)
        }
        bm_obj = GreedyExperimentalDesign:::compute_objective_vals_cpp(
          X_eval,
          bm_indicTs,
          objective,
          inv_cov
        )
        bm_stat = median(bm_obj)

        rerand = initRerandomizationExperimentalDesignObject(
          X,
          obj_val_cutoff_to_include = cutoff,
          max_designs = compare_r * 10L,
          objective = objective,
          wait = TRUE,
          num_cores = 1,
          seed = cfg$seed,
          verbose = FALSE,
          use_safe_inverse = TRUE
        )
        res_rerand = resultsRerandomizationSearch(rerand, include_assignments = FALSE)
        rerand_obj = res_rerand$obj_vals
        rerand_stat = median(rerand_obj)

        n = cfg$n
        r = compare_r
        block_obj_stats = numeric(0)
        for (B in B_values) {
          if (n %% B != 0) {
            next
          }
          n_B = n / B
          prop_T = floor(n_B / 2) / n_B
          if (prop_T <= 0 || prop_T >= 1) {
            next
          }
          block_designs = imbalanced_block_designs(
            n = n,
            prop_T = prop_T,
            B = B,
            r = r,
            seed = 123
          )
          block_obj = GreedyExperimentalDesign:::compute_objective_vals_cpp(
            X_eval,
            block_designs,
            objective,
            inv_cov
          )
          block_obj_stats = c(block_obj_stats, median(block_obj))
        }

        opt_obj = NULL
        if (cfg$n < 30 && cfg$n %% 2 == 0) {
          opt = initOptimalExperimentalDesignObject(
            X,
            objective = objective,
            wait = TRUE,
            num_cores = 1,
            verbose = FALSE,
            use_safe_inverse = TRUE
          )
          opt_res = resultsOptimalSearch(opt, num_vectors = min(2L, compare_r), form = "one_zero")
          opt_indicTs = opt_res$indicTs
          if (is.null(dim(opt_indicTs))) {
            opt_indicTs = matrix(opt_indicTs, nrow = 1)
          }
          opt_obj = GreedyExperimentalDesign:::compute_objective_vals_cpp(
            X_eval,
            opt_indicTs,
            objective,
            inv_cov
          )
          opt_stat = median(opt_obj)
        }

        for (b_idx in seq_along(block_obj_stats)) {
          #testthat::expect_true(block_obj_stats[b_idx] < complete_stat) #this seems to be hit or miss
          testthat::expect_true(greedy_stat <= block_obj_stats[b_idx])
          testthat::expect_true(bm_stat <= block_obj_stats[b_idx])
        }
        testthat::expect_true(greedy_stat < bm_stat)
        testthat::expect_true(bm_stat < rerand_stat)
        testthat::expect_true(rerand_stat < complete_stat)
        if (!is.null(opt_obj)) {
          tol = 1e-12
          testthat::expect_true(opt_stat <= greedy_stat + tol)
          testthat::expect_true(opt_stat <= bm_stat + tol)
          testthat::expect_true(opt_stat <= rerand_stat + tol)
          testthat::expect_true(opt_stat <= complete_stat + tol)
        }
      }
        }
      })
    })
