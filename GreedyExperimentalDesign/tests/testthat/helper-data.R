make_X <- function(n, p, seed = 1) {
  set.seed(seed)
  matrix(rnorm(n * p), nrow = n, ncol = p)
}

expect_binary_matrix <- function(mat, nrow_expected = NULL, ncol_expected = NULL) {
  if (!is.null(nrow_expected)) {
    testthat::expect_equal(nrow(mat), nrow_expected)
  }
  if (!is.null(ncol_expected)) {
    testthat::expect_equal(ncol(mat), ncol_expected)
  }
  testthat::expect_true(all(mat %in% c(0, 1)))
}

with_plot_device <- function(expr) {
  file = tempfile(fileext = ".png")
  png(file)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}
