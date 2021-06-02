options(java.parameters = c("-Xmx10000m"))
pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table, ggthemes, multcompView)

n = 40
p = 1
nC = 1
betaT = 1
sigma_e = 0.1
Nsim = 1000

res = data.table(design = character(), model = character(), n = factor(), p = factor(), imbalance = numeric(), betaThat = numeric())

set.seed(1984)
#switch this to optimal
# optimal_design = initOptimalExperimentalDesignObject(x, wait = TRUE, num_cores = nC)
# w_opt = resultsOptimalSearch(optimal_design, num_vectors = 2)$indicTs[1, ]
# opt_imbalance = compute_objective_val(x, w_opt)
# rm(optimal_design); gc()

for (nsim in 1 : Nsim){
  if (nsim %% 10 == 0){cat("nsim", nsim, "\n"); print(tail(res, 20))}
  
  x = as.matrix(sort(runif(n, 0, 3)))
  epsilons = rnorm(n, 0, sigma_e)
  
  binary_match_design = binaryMatchExperimentalDesignSearch(x)
  w = resultsBinaryMatchSearch(binary_match_design, num_vectors = 1)$indicTs[1, ]
  imbalance = compute_objective_val(x, w)
  
  # res = resultsBinaryMatchSearch(binary_match_design, num_vectors = 1, objective = "abs_sum_diff")
  # w = res$indicTs[1, ]
  # compute_objective_val(x, w)
  
  y_linear = x + w * betaT + epsilons
  betaThat = mean(y_linear[w == 1]) - mean(y_linear[w == 0])
  res = rbind(res, data.table(design = "matching", model = "linear", n = n, p = p, imbalance = imbalance, betaThat = betaThat))
  
  y_nonlinear = x^2 + w * betaT + epsilons
  betaThat = mean(y_nonlinear[w == 1]) - mean(y_nonlinear[w == 0])
  res = rbind(res, data.table(design = "matching", model = "nonlinear", n = n, p = p, imbalance = imbalance, betaThat = betaThat))
  
  # y_linear = x + w_opt * betaT + epsilons
  # betaThat = mean(y_linear[w_opt == 1]) - mean(y_linear[w_opt == 0])
  # res = rbind(res, data.table(design = "optimal", model = "linear", n = n, p = p, imbalance = opt_imbalance, betaThat = betaThat))
  # 
  # y_nonlinear = x^2 + w_opt * betaT + epsilons
  # betaThat = mean(y_nonlinear[w_opt == 1]) - mean(y_nonlinear[w_opt == 0])
  # res = rbind(res, data.table(design = "optimal", model = "nonlinear", n = n, p = p, imbalance = opt_imbalance, betaThat = betaThat))
  
  greedy_design = initGreedyExperimentalDesignObject(x, wait = TRUE, max_designs = 1, num_cores = nC)
  w = resultsGreedySearch(greedy_design, max_vectors = 1)$ending_indicTs[1, ]
  imbalance = compute_objective_val(x, w)
  
  y_linear = x + w * betaT + epsilons
  betaThat = mean(y_linear[w == 1]) - mean(y_linear[w == 0])
  res = rbind(res, data.table(design = "greedy_pair", model = "linear", n = n, p = p, imbalance = imbalance, betaThat = betaThat))
  
  y_nonlinear = x^2 + w * betaT + epsilons
  betaThat = mean(y_nonlinear[w == 1]) - mean(y_nonlinear[w == 0])
  res = rbind(res, data.table(design = "greedy_pair", model = "nonlinear", n = n, p = p, imbalance = imbalance, betaThat = betaThat))

  matching_then_greedy_design = binaryMatchFollowedByGreedyExperimentalDesignSearch(x, max_designs = 1, wait = TRUE, num_cores = nC)
  w = resultsBinaryMatchThenGreedySearch(matching_then_greedy_design)$indicTs[1, ]
  imbalance = compute_objective_val(x, w)
  
  y_linear = x + w * betaT + epsilons
  betaThat = mean(y_linear[w == 1]) - mean(y_linear[w == 0])
  res = rbind(res, data.table(design = "matching_then_greedy_pair", model = "linear", n = n, p = p, imbalance = imbalance, betaThat = betaThat))
  
  y_nonlinear = x^2 + w * betaT + epsilons
  betaThat = mean(y_nonlinear[w == 1]) - mean(y_nonlinear[w == 0])
  res = rbind(res, data.table(design = "matching_then_greedy_pair", model = "nonlinear", n = n, p = p, imbalance = imbalance, betaThat = betaThat))
  
}
res[, squared_error := (betaThat - betaT)^2]
res = res[imbalance != 0]
avg_summary = res %>% 
  group_by(design, model) %>% 
  summarize(avg_log10_imbalance = mean(log10(imbalance)), mse = mean(squared_error))
print(ggplot(res) +
        geom_density(aes(x = squared_error), position = "identity", fill = "blue", alpha = 0.2) +
        facet_wrap(design ~ model, 
                   # scales = "free_x", 
                   ncol = 2,
                   labeller = function(labs){label_both(labs, multi_line = FALSE)}) + 
        # theme(axis.line=element_line()) +
        scale_x_log10() +
        # xlim(-3, 2) + 
        xlab("estimator squared error") +
        geom_vline(data = avg_summary, aes(xintercept = mse), lwd = 1, alpha = 0.5))

print(ggplot(res) +
        geom_density(aes(x = log10(imbalance)), position = "identity", fill = "blue", alpha = 0.2) +
        facet_wrap(design ~ ., 
                   # scales = "free_x", 
                   ncol = 3,
                   labeller = function(labs){label_both(labs, multi_line = FALSE)}) + 
        # theme(axis.line=element_line()) +
        # scale_x_log10() +
        # xlim(-3, 2) + 
        xlab("log10 balance") +
        geom_vline(data = avg_summary, aes(xintercept = avg_log10_imbalance), lwd = 1, alpha = 0.5))

# save(res, file = "demo_mse_linear_nonlinear2.RData")


mod_linear = aov(squared_error ~ design, res[model == "linear"])
# summary(mod_linear)
TukeyHSD(x = mod_linear, which = "design", conf.level = 0.95)


mod_nonlinear = aov(squared_error ~ design, res[model == "nonlinear"])
# summary(mod_nonlinear)
TukeyHSD(x = mod_nonlinear, which = "design", conf.level = 0.95)

avg_summary[order(avg_summary$model), ]
