options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 20, p = 3)
ged = initGreedyExperimentalDesignObject(X, max_designs = 4, num_cores = 4, objective = "abs_sum_diff", wait = TRUE, diagnostics = FALSE)
res = resultsGreedySearch(ged)
res