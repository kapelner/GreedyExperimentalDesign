#setwd("C:\\Users\\kapelner\\workspace\\GreedyExperimentalDesign")
#library(roxygen2)
#roxygenise("GreedyExperimentalDesign", clean = TRUE)

options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

#X = generate_stdzied_design_matrix(n = 200, p = 40)
ged = initGreedyExperimentalDesignObject(X, max_designs = 50000, num_cores = 2, objective = "abs_sum_diff", semigreedy = FALSE)
startGreedySearch(ged)
resultsGreedySearch(ged)$obj_vals[1]

ged; resultsGreedySearch(ged)$obj_vals[1]; resultsGreedySearch(ged)$obj_vals[100]

data = plot_obj_val_order_statistic(ged, order_stat = 200)
vals = data$val_order_stats
plot(1:length(vals), log(log(vals)))
plot(1:length(vals), vals^(-15))

res = resultsGreedySearch(ged)