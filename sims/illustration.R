n = 28
n_T = 14
choose(n, n_T)

X = matrix(rnorm(n), nrow = n)

options(java.parameters = "-Xmx5000m") 
library(GreedyExperimentalDesign)
opt = initOptimalExperimentalDesignObject(X, wait = TRUE)

res_opt = resultsOptimalSearch(opt)
all_obj_vals = res_opt$all_obj_vals
ci_cr = quantile(all_obj_vals, c(0.025, 0.975))
avg_cr = mean(all_obj_vals)

greedy = initGreedyExperimentalDesignObject(X, wait = TRUE)
res_gr = resultsGreedySearch(greedy)
ci_greedy_res = quantile(res_gr$obj_vals, c(0.025, 0.975))
avg_gr = mean(res_gr$obj_vals)

###rerand
Nrerand = 1e3
RerandSize = 1e5
obj_vals_rerand = array(NA, Nrerand)
for (i in 1 : Nrerand){
  obj_vals_rerand[i] = min(sample(all_obj_vals, RerandSize))
}

ci_rerand = quantile(obj_vals_rerand, c(0.025, 0.975))
avg_rerand = mean(obj_vals_rerand)

hist(all_obj_vals, br = 1000)
abline(v = ci_cr[1], col = "red")
abline(v = avg_cr, col = "red")
abline(v = ci_cr[2], col = "red")
abline(v = ci_greedy_res[1], col = "green")
abline(v = avg_gr, col = "green")
abline(v = ci_greedy_res[2], col = "green")
abline(v = ci_rerand[1], col = "orange")
abline(v = avg_rerand, col = "orange")
abline(v = ci_rerand[2], col = "orange")

hist(log10(all_obj_vals), br = 1000)
abline(v = log10(ci_cr[1]), col = "red")
abline(v = log10(avg_cr), col = "red")
abline(v = log10(ci_cr[2]), col = "red")
abline(v = log10(ci_greedy_res[1]), col = "green")
abline(v = log10(avg_gr), col = "green")
abline(v = log10(ci_greedy_res[2]), col = "green")
abline(v = log10(res_opt$opt_obj_val), col = "purple")
abline(v = log10(ci_rerand[1]), col = "orange")
abline(v = log10(avg_rerand), col = "orange")
abline(v = log10(ci_rerand[2]), col = "orange")
