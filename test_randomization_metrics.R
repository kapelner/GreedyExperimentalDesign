
options(java.parameters = "-Xmx3000m")
library(GreedyExperimentalDesign)

n = 100
p = 1
r = 100000
X = generate_stdzied_design_matrix(n = n, p = p)
rd = initRerandomizationExperimentalDesignObject(X, r, wait = TRUE)
designs = sapply(.jcall(rd$java_obj, "[[I", "getEndingIndicTs"), .jevalArray)

metr = compute_randomization_metrics(designs)
metr$rand_entropy_metric
metr$rand_norm_se_metric


gd = initGreedyExperimentalDesignObject(X, r, wait = TRUE)
res = resultsGreedySearch(gd)
designs = res$ending_indicTs
metr = compute_randomization_metrics(designs)
metr$rand_entropy_metric
metr$rand_norm_se_metric

