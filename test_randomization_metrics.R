
options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

X = generate_stdzied_design_matrix(n = 14, p = 1)
rd = initRerandomizationExperimentalDesignObject(X, 10)
res = resultsRerandomizationSearch(rd)
designs = res$ending_indicTs
designs
compute_randomization_metrics(designs)
n = nrow(designs)
r = ncol(designs)

java_obj = .jnew("ExperimentalDesign.RandomizationMetrics")
.jcall(java_obj, "V", "setNandR", as.integer(n), as.integer(r))
for (j in 1 : r){
  .jcall(java_obj, "V", "setDesign", as.integer(r - 1), designs[, j, drop = FALSE]) #java indexes from 0...n-1
}
.jcall(java_obj, "V", "compute")

#harvest the data and return it as a list
p_hat_ijs = sapply(.jcall(java_obj, "[[D", "getPhats"), .jevalArray)
rand_entropy_metric = .jcall(obj$java_obj, "D", "getRandEntropyMetric")
rand_norm_se_metric = .jcall(obj$java_obj, "D", "getRandStdErrMetric")
list(
  p_hat_ijs = p_hat_ijs, 
  rand_entropy_metric = rand_entropy_metric, 
  rand_norm_se_metric = rand_norm_se_metric
)