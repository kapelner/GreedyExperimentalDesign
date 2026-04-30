options(java.parameters = c("-Xmx10g", "--enable-native-access=ALL-UNNAMED"))

library(testthat)
library(GreedyExperimentalDesign)

is_check = nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
  nzchar(Sys.getenv("R_CMD_CHECK")) ||
  nzchar(Sys.getenv("_R_CHECK_CRAN_INCOMING_")) ||
  nzchar(Sys.getenv("_R_CHECK_CRAN_INCOMING_REMOTE_"))
if (is_check) {
  message("Skipping testthat tests during R CMD check.")
} else {
  test_check("GreedyExperimentalDesign")
}
