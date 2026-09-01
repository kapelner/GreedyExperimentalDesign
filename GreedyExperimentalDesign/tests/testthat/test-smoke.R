test_that("the packaged Java backend completes a small greedy search", {
  jar = system.file("java", "GreedyExperimentalDesign.jar", package = "GreedyExperimentalDesign")
  expect_true(nzchar(jar))
  expect_true(file.exists(jar))

  set.seed(1)
  X = matrix(rnorm(24), nrow = 8, ncol = 3)
  search = initGreedyExperimentalDesignObject(
    X,
    max_designs = 2,
    objective = "abs_sum_diff",
    wait = TRUE,
    num_cores = 1,
    seed = 1,
    verbose = FALSE
  )
  result = resultsGreedySearch(search, max_vectors = 2)

  expect_length(result$obj_vals, 2)
  expect_true(all(is.finite(result$obj_vals)))
  expect_equal(dim(result$ending_indicTs), c(2L, 8L))
  expect_true(all(result$ending_indicTs %in% c(0, 1)))
})
