source("helper-package.R")
test_that("core functions work across sizes", {
  skip_on_cmd_check()
  with_immediate_failures({
    configs_raw = build_test_configs()
    config_keys = vapply(configs_raw, function(cfg) {
      paste(cfg$n, cfg$p, cfg$r, cfg$prop_T, sep = "|")
    }, character(1))
    configs = configs_raw[!duplicated(config_keys)]
    B_values = sort(unique(vapply(configs_raw, function(cfg) cfg$B, numeric(1))))
    objectives = c("abs_sum_diff", "mahal_dist")

    for (objective in objectives) {
      for (cfg in configs) {
        if (objective == "mahal_dist" && cfg$p >= 10) { #too slow
          next
        }
        X = make_X(cfg$n, cfg$p, seed = cfg$seed)
        X1 = X[, 1, drop = FALSE]
        max_vectors = max(1L, min(as.integer(cfg$r), 10L))
        cat(
          "objective=", objective,
          " n=", cfg$n,
          " p=", cfg$p,
          " r=", cfg$r,
          " prop_T=", cfg$prop_T,
          " seed=", cfg$seed,
          "\n",
          sep = ""
        )

    Xstd = standardize_data_matrix(X)
    expect_equal(dim(Xstd), dim(X))
    expect_true(all(abs(colMeans(Xstd)) < 1e-10))
    expect_true(all(abs(apply(Xstd, 2, sd) - 1) < 1e-10))

    Xgen = generate_stdzied_design_matrix(n = cfg$n, p = cfg$p, covariate_gen = rnorm)
    expect_equal(dim(Xgen), dim(X))
    expect_true(all(abs(colMeans(Xgen)) < 1e-10))

    w = 1:cfg$n
    w1 = shuffle_cpp_wrap(w, seed = cfg$seed)
    w2 = shuffle_cpp_wrap(w, seed = cfg$seed)
    expect_equal(w1, w2)
    expect_equal(sort(w1), w)
    expect_true(all_elements_same_cpp_wrap(rep(2, 5)))
    expect_false(all_elements_same_cpp_wrap(c(1, 2, 1)))

    Xsmall = matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
    D = compute_distance_matrix_cpp_wrap(Xsmall)
    expect_equal(dim(D), c(3, 3))
    expect_true(is.na(diag(D)[1]))
    expect_equal(D[1, 2], D[2, 1])
    K = compute_gram_matrix(Xsmall, kernel_type = "vanilla")
    expect_equal(dim(K), c(3, 3))
    expect_equal(K, t(K))

    indic = c(rep(1, cfg$n / 2), rep(0, cfg$n / 2))
    obj_val = compute_objective_val(
      X,
      indic,
      objective = objective,
      use_safe_inverse = TRUE
    )
    expect_true(is.numeric(obj_val))

    designs = complete_randomization(n = cfg$n, r = cfg$r)
    metrics = compute_randomization_metrics(t(designs))
    expect_true(is.list(metrics))
    expect_equal(dim(metrics$p_hat_ijs), c(cfg$n, cfg$n))
    expect_true(is.numeric(metrics$rand_entropy_metric))

    mat = complete_randomization(n = cfg$n, r = cfg$r)
    expect_binary_matrix(mat, nrow_expected = cfg$r, ncol_expected = cfg$n)

    forced = complete_randomization_with_forced_balanced(n = cfg$n, r = cfg$r, seed = cfg$seed)
    expect_binary_matrix(forced, nrow_expected = cfg$r, ncol_expected = cfg$n)
    expect_equal(rowSums(forced), rep(cfg$n / 2, cfg$r))

    imbal = imbalanced_complete_randomization(n = cfg$n, prop_T = cfg$prop_T, r = cfg$r, seed = cfg$seed)
    expect_binary_matrix(imbal, nrow_expected = cfg$r, ncol_expected = cfg$n)
    n_T = as.integer(round(cfg$n * cfg$prop_T))
    expect_equal(rowSums(imbal), rep(n_T, cfg$r))

    for (B in B_values) {
      if (cfg$n %% B != 0) {
        next
      }
      n_B = cfg$n / B
      if (abs(n_B * cfg$prop_T - round(n_B * cfg$prop_T)) > 1e-8) {
        next
      }
      blocks = imbalanced_block_designs(
        n = cfg$n,
        prop_T = cfg$prop_T,
        B = B,
        r = cfg$r,
        seed = cfg$seed
      )
      expect_binary_matrix(blocks, nrow_expected = cfg$r, ncol_expected = cfg$n)
      expect_equal(rowSums(blocks), rep(n_T, cfg$r))

      dummy_block = rep(c(1, 0), length.out = cfg$n / B)
      block_designs = generate_block_design_cpp_wrap(B = B, nR = cfg$r, dummy_block = dummy_block)
      expect_equal(dim(block_designs), c(cfg$r, cfg$n))
      expect_equal(rowSums(block_designs), rep(sum(dummy_block) * B, cfg$r))

      SigmaW = gen_var_cov_matrix_block_designs(n = cfg$n, prop_T = cfg$prop_T, B = B)
      expect_equal(dim(SigmaW), c(cfg$n, cfg$n))
      expect_equal(SigmaW, t(SigmaW))
    }

    pairs = matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
    Wpm = gen_pm_designs_cpp_wrap(pairs, n = 2, r = cfg$r)
    expect_equal(dim(Wpm), c(cfg$r, 4))
    expect_true(all(Wpm %in% c(-1, 1)))

    pCs = c(0.1, 0.2, 0.3, 0.4)
    pTs = c(0.9, 0.8, 0.7, 0.6)
    Wbin = matrix(c(1, 0, 1, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE)
    Y = create_all_ys_cpp_wrap(pCs, pTs, Wbin, two_n = 4, nY = 2)
    expect_equal(dim(Y), c(2, 4))
    expect_true(is.logical(Y))

    res1 = optimize_asymmetric_treatment_assignment(n = cfg$n)
    expect_equal(res1$nT + res1$nC, cfg$n)
    res2 = optimize_asymmetric_treatment_assignment(n = cfg$n, c_treatment = 2, c_control = 1)
    expect_equal(res2$nT + res2$nC, cfg$n)
    expect_true(res2$nT > res2$nC)

    if (cfg$n >= 4) {
      Xhad = X[1:4, , drop = FALSE]
      W = hadamardExperimentalDesign(Xhad, strict = TRUE, form = "one_zero")
      expect_true(is.matrix(W))
      expect_true(all(W %in% c(0, 1)))
      expect_true(all(rowSums(W) == nrow(Xhad) / 2))
    }

    r_cur = max(2, cfg$r)
    Wcur = matrix(sample(c(-1, 1), r_cur * cfg$n, replace = TRUE), nrow = r_cur, ncol = cfg$n)
    res_cur = greedy_orthogonalization_curation(Wcur, Rmin = 2, verbose = FALSE)
    expect_true(is.list(res_cur))
    expect_equal(ncol(res_cur$Wsorted), cfg$n)
    R0 = min(max(2, floor(r_cur / 2)), r_cur)
    Wcur2 = greedy_orthogonalization_curation2(Wcur, R0 = R0, verbose = FALSE)
    expect_true(is.matrix(Wcur2))
    expect_equal(ncol(Wcur2), cfg$n)

    ged = initGreedyExperimentalDesignObject(
      X,
      max_designs = cfg$r,
      objective = objective,
      diagnostics = TRUE,
      wait = TRUE,
      num_cores = 1,
      seed = cfg$seed,
      verbose = FALSE,
      use_safe_inverse = TRUE
    )
    expect_output(print(ged))
    expect_output(summary(ged))
    res_greedy = resultsGreedySearch(ged, max_vectors = max_vectors)
    expect_binary_matrix(res_greedy$ending_indicTs, nrow_expected = max_vectors, ncol_expected = cfg$n)
    expect_equal(rowSums(res_greedy$ending_indicTs), rep(cfg$n / 2, max_vectors))
    with_plot_device(plot_obj_val_by_iter(res_greedy))
    with_plot_device(plot_obj_val_order_statistic(ged, order_stat = 1, skip_every = 1))
    with_plot_device(plot(ged))

    mk = suppressWarnings(initGreedyMultipleKernelExperimentalDesignObject(
      X,
      max_designs = cfg$r,
      kernel_pre_num_designs = cfg$r,
      kernel_names = c("mahalanobis", "gaussian"),
      wait = TRUE,
      num_cores = 1,
      seed = cfg$seed,
      diagnostics = TRUE,
      verbose = FALSE,
      use_safe_inverse = TRUE
    ))
    expect_output(print(mk))
    expect_output(summary(mk))
    res_mk = resultsMultipleKernelGreedySearch(mk, max_vectors = max_vectors, form = "one_zero")
    expect_binary_matrix(res_mk$ending_indicTs, nrow_expected = max_vectors, ncol_expected = cfg$n)
    expect_equal(rowSums(res_mk$ending_indicTs), rep(cfg$n / 2, max_vectors))
    with_plot_device(plot(mk))

    bms = computeBinaryMatchStructure(X, use_safe_inverse = TRUE)
    bms$verbose = FALSE
    bm_unique = 2^(cfg$n / 2)
    bm_max_designs = min(cfg$r, max(1L, floor(bm_unique / 2)))
    bm = initBinaryMatchExperimentalDesignSearchObject(
      bms,
      max_designs = bm_max_designs,
      wait = TRUE,
      num_cores = 1,
      seed = cfg$seed,
      verbose = FALSE
    )
    expect_output(print(bm))
    expect_output(summary(bm))
    res_bm = resultsBinaryMatchSearch(bm, form = "one_zero")
    expect_binary_matrix(res_bm, nrow_expected = bm_max_designs, ncol_expected = cfg$n)
    expect_equal(rowSums(res_bm), rep(cfg$n / 2, bm_max_designs))

    if (cfg$n %% 4 == 0) {
      bm_greedy = initBinaryMatchFollowedByGreedyExperimentalDesignSearchObject(
        X,
        diff_method = TRUE,
        max_designs = cfg$r,
        wait = TRUE,
        num_cores = 1,
        seed = cfg$seed,
        objective = objective,
        diagnostics = TRUE,
        verbose = FALSE
      )
      expect_output(print(bm_greedy))
      expect_output(summary(bm_greedy))
      res_bm_greedy = resultsBinaryMatchThenGreedySearch(
        bm_greedy,
        num_vectors = max_vectors,
        compute_obj_vals = TRUE,
        use_safe_inverse = TRUE
      )
      expect_binary_matrix(res_bm_greedy$indicTs, nrow_expected = max_vectors, ncol_expected = cfg$n)
      expect_equal(rowSums(res_bm_greedy$indicTs), rep(cfg$n / 2, max_vectors))

      bm_rerand = initBinaryMatchFollowedByRerandomizationDesignSearchObject(
        X,
        max_designs = cfg$r,
        wait = TRUE,
        num_cores = 1,
        seed = cfg$seed,
        objective = objective,
        obj_val_cutoff_to_include = Inf,
        verbose = FALSE
      )
      expect_output(print(bm_rerand))
      expect_output(summary(bm_rerand))
      res_bm_rerand = resultsBinaryMatchThenRerandomizationSearch(
        bm_rerand,
        num_vectors = max_vectors,
        compute_obj_vals = TRUE,
        use_safe_inverse = TRUE
      )
      expect_binary_matrix(res_bm_rerand$indicTs, ncol_expected = cfg$n)
      expect_true(nrow(res_bm_rerand$indicTs) >= max_vectors)
      expect_equal(rowSums(res_bm_rerand$indicTs), rep(cfg$n / 2, nrow(res_bm_rerand$indicTs)))
    } else {
      expect_error(
        initBinaryMatchFollowedByGreedyExperimentalDesignSearchObject(
          X,
          diff_method = TRUE,
          max_designs = cfg$r,
          wait = TRUE,
          num_cores = 1,
          seed = cfg$seed,
          objective = objective,
          diagnostics = TRUE,
          verbose = FALSE
        ),
        "divisible by four"
      )
      expect_error(
        initBinaryMatchFollowedByRerandomizationDesignSearchObject(
          X,
          max_designs = cfg$r,
          wait = TRUE,
          num_cores = 1,
          seed = cfg$seed,
          objective = objective,
          obj_val_cutoff_to_include = Inf,
          verbose = FALSE
        ),
        "divisible by four"
      )
    }

    rd = initRerandomizationExperimentalDesignObject(
      X,
      obj_val_cutoff_to_include = Inf,
      max_designs = cfg$r,
      objective = objective,
      wait = TRUE,
      num_cores = 1,
      seed = cfg$seed,
      verbose = FALSE,
      use_safe_inverse = TRUE
    )
    expect_output(print(rd))
    expect_output(summary(rd))
    res_rd = resultsRerandomizationSearch(rd, include_assignments = TRUE, form = "one_zero")
    expect_equal(dim(res_rd$ending_indicTs), c(cfg$r, cfg$n))
    expect_equal(rowSums(res_rd$ending_indicTs), rep(cfg$n / 2, cfg$r))

    Xopt = X[1:4, , drop = FALSE]
    if (objective == "abs_sum_diff" || cfg$p < nrow(Xopt)) {
      od = initOptimalExperimentalDesignObject(
        Xopt,
        objective = objective,
        wait = TRUE,
        num_cores = 1,
        verbose = FALSE,
        use_safe_inverse = TRUE
      )
      expect_output(print(od))
      expect_output(summary(od))
      res_opt = resultsOptimalSearch(od, num_vectors = 2, form = "one_zero")
      expect_equal(dim(res_opt$indicTs), c(2, 4))
    }

    if (objective == "abs_sum_diff" && cfg$p == 1) {
      kd = initKarpExperimentalDesignObject(X1, wait = TRUE, balanced = TRUE, verbose = FALSE)
      expect_output(print(kd))
      expect_output(summary(kd))
      res_kd = resultsKarpSearch(kd)
      expect_equal(length(res_kd$indicT), cfg$n)
    }
      }
    }
  })
})
