options(java.parameters = c("-Xmx2g", "--enable-native-access=ALL-UNNAMED", "--enable-preview"))

library(testthat)
library(GreedyExperimentalDesign)

test_check("GreedyExperimentalDesign")
