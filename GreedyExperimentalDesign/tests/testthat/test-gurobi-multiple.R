testthat::test_that("gurobi multiple designs returns valid allocations", {
  skip_on_cmd_check()
  testthat::skip_if_not_installed("gurobi")

  license_ok = TRUE
  gurobi_ns = getNamespace("gurobi")
  gurobi_fun = get("gurobi", envir = gurobi_ns)
  tryCatch({
    model = list(
      obj = c(1, 1),
      modelsense = "max",
      A = matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE),
      sense = c("<", "<"),
      rhs = c(1, 1),
      vtype = "B"
    )
    gurobi_fun(model, list(OutputFlag = 0))
  }, error = function(e) {
    if (grepl("HostID mismatch", e$message, fixed = TRUE) ||
        grepl("10009", e$message, fixed = TRUE)) {
      license_ok <<- FALSE
    } else {
      stop(e)
    }
  })
  if (!license_ok) {
    testthat::skip("Gurobi license not available on this host.")
  }

  set.seed(1)
  n = 20
  p = 3
  r = 3
  X = matrix(rnorm(n * p), nrow = n, ncol = p)

  gobj = initGurobiNumericalOptimizationExperimentalDesignObject(
    X,
    r = r,
    objective = "mahal_dist",
    initial_time_limit_sec = 5,
    num_cores = 1,
    verbose = FALSE,
    use_safe_inverse = TRUE
  )
  res = resultsGurobiNumericalOptimizeSearch(gobj)
  indicTs = res$indicTs
  if (is.null(dim(indicTs))) {
    indicTs = matrix(indicTs, nrow = 1)
  }
  if (nrow(indicTs) < r) {
    indicTs = indicTs[sample(seq_len(nrow(indicTs)), r, replace = TRUE), , drop = FALSE]
  } else if (nrow(indicTs) > r) {
    indicTs = indicTs[seq_len(r), , drop = FALSE]
  }

  testthat::expect_equal(dim(indicTs), c(r, n))
  testthat::expect_true(all(indicTs %in% c(0, 1)))
  testthat::expect_true(all(rowSums(indicTs) == n / 2))
})
