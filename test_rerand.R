
options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 100, p = 1)
rd = initRerandomizationExperimentalDesignObject(X, 1000, 0.000000001)
res = resultsRerandomizationSearch(rd)
res$obj_vals
hist(res$obj_vals, br = 100)
dim(res$ending_indicTs)