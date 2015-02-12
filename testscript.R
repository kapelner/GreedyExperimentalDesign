setwd("C:\\Users\\kapelner\\workspace\\GreedyExperimentalDesign")
library(roxygen2)
roxygenise("GreedyExperimentalDesign", clean = TRUE)

options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 200, p = 40)
ged = initGreedyExperimentalDesignObject(X, max_designs = 50000, num_cores = 2, objective = "abs_sum_diff", semigreedy = FALSE)
startGreedySearch(ged)
resultsGreedySearch(ged)$obj_vals[1]

ged; resultsGreedySearch(ged)$obj_vals[1]; resultsGreedySearch(ged)$obj_vals[100]

data = plot_obj_val_order_statistic(ged, order_stat = 200)
vals = data$val_order_stats
plot(1:length(vals), log(log(vals)))
plot(1:length(vals), vals^(-15))

res = resultsGreedySearch(ged)


options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 200, p = 40)
ged = initGreedyExperimentalDesignObject(X, max_designs = 100, num_cores = 2, objective = "mahal_dist")
startGreedySearch(ged)

data = plot_obj_val_order_statistic(ged, order_stat = 5, skip_every = 1)

#check to make sure Java and R agree
indicT = resultsGreedySearch(ged)$indicTs[, 1]
xbartminxbarc = colMeans(X[indicT == 1, ]) - colMeans(X[indicT == 0, ])
t(xbartminxbarc) %*% solve(var(X)) %*% t(t(xbartminxbarc))
resultsGreedySearch(ged)$obj_vals[1]

ged = initGreedyExperimentalDesignObject(X, max_designs = 100, num_cores = 2)
startGreedySearch(ged)

#check to make sure Java and R agree
indicT = resultsGreedySearch(ged)$indicTs[, 1]
xbartminxbarc = colMeans(X[indicT == 1, ]) - colMeans(X[indicT == 0, ])
sum(abs(xbartminxbarc))
resultsGreedySearch(ged)$obj_vals[1]

options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 100, p = 1)
ged = initGreedyExperimentalDesignObject(X, max_designs = 50, num_cores = 6, wait = TRUE)
startGreedySearch(ged)
res = resultsGreedySearch(ged)
res

ged = initGreedyExperimentalDesignObject(X, max_designs = 50, num_cores = 6, wait = TRUE, max_iters = 1)
startGreedySearch(ged)
res = resultsGreedySearch(ged)
res
