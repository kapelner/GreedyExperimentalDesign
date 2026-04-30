test_that("resultsMultipleKernelGreedySearch returns best results in order", {
  skip_on_cmd_check()
  library(GreedyExperimentalDesign)
  
  N = 40
  P = 2
  MAX_DESIGNS = 20
  NUM_CORES = 1
  SEED = 1
  
  X = matrix(rnorm(N * P), nrow = N)
  
  mk = suppressWarnings(initGreedyMultipleKernelExperimentalDesignObject(
    X,
    max_designs = MAX_DESIGNS,
    kernel_pre_num_designs = 5,
    num_cores = NUM_CORES,
    kernel_names = c("mahalanobis", "gaussian"),
    start = TRUE,
    wait = TRUE,
    seed = SEED,
    diagnostics = FALSE,
    verbose = FALSE
  ))
  
  # Return fewer than all designs
  request_vectors = 5
  res = resultsMultipleKernelGreedySearch(mk, max_vectors = request_vectors)
  
  # Check dimensions
  expect_equal(nrow(res$ending_indicTs), request_vectors)
  expect_equal(length(res$obj_vals), request_vectors)
  
  # Check ordering (best to worst, so increasing objective values)
  expect_true(all(diff(res$obj_vals) >= 0))
  
  # Check that they are indeed the best ones
  all_obj_vals = .jcall(mk$java_obj, "[D", "getObjectiveVals")
  expect_equal(res$obj_vals, sort(all_obj_vals)[1:request_vectors])
})
