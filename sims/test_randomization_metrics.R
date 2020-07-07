
options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

n = 100
p = 1
r = 10000
X = generate_stdzied_design_matrix(n = n, p = p)
rd = initRerandomizationExperimentalDesignObject(X, r, wait = TRUE)
designs = .jcall(rd$java_obj, "[[I", "getEndingIndicTs", simplify = TRUE)

metr = compute_randomization_metrics(designs)
metr$rand_entropy_metric
metr$rand_norm_se_metric
java_obj = .jnew("DesignMetrics.RandomizationMetrics")
.jcall(java_obj, "V", "setNandR", as.integer(n), as.integer(r))
for (j in 1 : r){
  .jcall(java_obj, "V", "setDesign", as.integer(j - 1), as.integer(designs[, j])) #java indexes from 0...n-1
}
.jcall(java_obj, "V", "compute")
p_hat_ijs = .jcall(java_obj, "[[D", "getPhats", simplify = TRUE)
phats = c(metr$p_hat_ijs)
phats = phats[phats!=0]
length(phats)

hist(phats, br=100)

s_n = (n-2) / (2 * n - 2)

sum((phats - s_n)^2)


gd = initGreedyExperimentalDesignObject(X, r, wait = TRUE)
res = resultsGreedySearch(gd)
designs = res$ending_indicTs
metr = compute_randomization_metrics(designs)
metr$rand_entropy_metric
metr$rand_norm_se_metric

designs = matrix(as.numeric(rep(sample(c(rep(1, n / 2), rep(0, n / 2))), r)), ncol = r)
metr = compute_randomization_metrics(designs)
metr$rand_entropy_metric
metr$rand_norm_se_metric


