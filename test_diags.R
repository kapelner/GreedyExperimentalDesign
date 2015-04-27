options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 200, p = 40)
ged = initGreedyExperimentalDesignObject(X, max_designs = 10, num_cores = 1, objective = "abs_sum_diff", wait = TRUE, diagnostics = TRUE)
resultsGreedySearch(ged)