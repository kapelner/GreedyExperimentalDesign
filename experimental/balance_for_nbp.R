options(java.parameters = c("-Xmx8000m"))
pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table ,ggthemes)

ns = c(20, 50, 76, 100, 126, 150, 200, 300)
ps = c(1, 2, 5, 7, 10, 15, 20, 25, 30, 40)
Nsim = 5000

res = data.table(n = factor(), p = factor(), imbalances = numeric())
for (n in ns){
  X = matrix(rnorm(n * max(ps)), nrow = n)
  for (p in ps){
    if (p > n){
      next
    }
    cat("n", n, "p", p, "\n")
    binary_match_design = binaryMatchExperimentalDesignSearch(X[, 1 : p, drop = FALSE])
    W = resultsBinaryMatchSearch(binary_match_design, min(Nsim, 2^(n/2)))$indicTs
    imbalances = apply(W, 1, function(w){compute_objective_val(X[, 1 : p, drop = FALSE], w)})
    res = rbind(res, data.table(n = n, p = p, imbalances = imbalances))
  }
}

avg_summary = res %>% group_by(n, p) %>% summarize(avg = mean(log10(imbalance)))
print(ggplot(res) +
        geom_density(aes(x = log10(imbalance)), position = "identity", fill = "blue", alpha = 0.2) +
        facet_wrap(n ~ p, 
                   # scales = "free_x", 
                   ncol = length(ps),
                   labeller = function(labs){label_both(labs, multi_line = FALSE)}) + 
        theme(axis.line=element_line()) +
        # scale_x_log10() + 
        xlim(-3, 2) + 
        xlab("abs sum diff imbalance in x") +
        geom_vline(data = avg_summary, aes(xintercept = avg), lwd = 1, alpha = 0.5))

save(res, file = "nbm_balances_n_p.RData")

summary(res$imbalances)

res$n = relevel(res$n, as.character(min(ns)))
res$p = relevel(res$p, as.character(min(ps)))

mod = lm(imbalances ~ p * n, res)
summary(mod)
mod = lm(log(imbalances) ~ p * n, res)
summary(mod)

res[, n_num := as.numeric(as.character(n))]
res[, p_num := as.numeric(as.character(p))]
mod = lm(imbalances ~ n_num * p_num, res)
summary(mod)

mod = lm(imbalances ~ log10(n_num) + log10(p_num), res)
summary(mod)
mod = lm(imbalances ~ log10(n_num) * log10(p_num), res)
summary(mod)
mod = lm((imbalances) ~ log10(n_num), res[p == 1])
summary(mod)
mod = lm((imbalances) ~ 0+log(n_num), res[p == 1])
summary(mod)

res[, imbalances_times_n := imbalances * n_num * log(n_num)]
res[, imbalances_times_n := imbalances * n_num]
ggplot(res[p == 1]) + geom_boxplot(aes(x = n, y = log(imbalances_times_n)))

mod = lm(imbalances_times_n ~ n_num, res[p == 1])
summary(mod)

