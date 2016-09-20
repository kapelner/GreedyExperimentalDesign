options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 100, p = 30)
ged = initGreedyExperimentalDesignObject(X, max_designs = 40, num_cores = 4, objective = "abs_sum_diff", wait = TRUE, diagnostics = TRUE)
res = resultsGreedySearch(ged)
plot_obj_val_by_iter(res)
names(res)

compute_randomization_metrics(res$starting_indicTs)


designs = matrix(rep(c(rep(0, 50), rep(1, 50)), 50), ncol = 50)	
compute_randomization_metrics(same)

