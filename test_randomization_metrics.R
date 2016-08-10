
options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 100, p = 1)
rd = initRerandomizationExperimentalDesignObject(X, 100)

designs = pure_randomization(n = 50, r = 100)
