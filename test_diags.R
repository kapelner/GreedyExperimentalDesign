options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 100, p = 30)
ged = initGreedyExperimentalDesignObject(X, max_designs = 40, num_cores = 4, objective = "abs_sum_diff", wait = TRUE, diagnostics = TRUE)
res = resultsGreedySearch(ged)
plot_obj_val_by_iter(res)
