options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign

source("illustration_params.R")
ns = c(32, 100, 1000)
ps = c(1, 2, 5, 10)
nR = 200
rerand_threshold = 0.01
nX = 50
nC = 3
objective = "mahal_dist"
res = data.frame(imbalance = numeric(), method = character(), n = numeric(), p = numeric(), nx = numeric())


set.seed(1)
for (nx in 1 : nX){
  for (n in ns){
    for (p in ps){
      cat("n", n, "p", p, "nx", nx, "\n")
      X = matrix(runif(n * p, -sqrt(3), sqrt(3)), nrow = n)
      inv_cov_X = solve(var(X))
      
      
      # X = data.matrix(MASS::Pima.tr[1 : n, 1 : p])
      #stdize
      # X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
      
      ###D = BCRD
      w_bcrds = complete_randomization_with_forced_balanced(n, nR)
      bcrd_imbalances = apply(w_bcrds, 1, function(w){compute_objective_val(X, w, objective = objective, inv_cov_X = inv_cov_X)})
      
      ###D = R
      w_rs = complete_randomization_with_forced_balanced(n, nR * 1 / rerand_threshold)
      bcrds_rerand_balances = apply(w_rs, 1, function(w){compute_objective_val(X, w, objective = objective, inv_cov_X = inv_cov_X)})
      r_imbalances = bcrds_rerand_balances[order(bcrds_rerand_balances)[1 : nR]]
      
      #D = M
      bmeds = binaryMatchExperimentalDesignSearch(X)
      w_ms = resultsBinaryMatchSearch(bmeds, num_vectors = nR)$indicTs
      m_imbalances = apply(w_ms, 1, function(w){compute_objective_val(X, w, objective = objective, inv_cov_X = inv_cov_X)})
      # dim(bmeds_res$indicTs)
      
      #D = G
      ged = initGreedyExperimentalDesignObject(X, objective = objective, max_designs = nR, wait = TRUE, num_cores = nC)
      w_gs = resultsGreedySearch(ged, max_vectors = nR)$ending_indicTs
      g_imbalances = apply(w_gs, 1, function(w){compute_objective_val(X, w, objective = objective, inv_cov_X = inv_cov_X)})
      
      #D = MG
      bmfged = binaryMatchFollowedByGreedyExperimentalDesignSearch(X, objective = objective, max_designs = nR, wait = TRUE, num_cores = nC)
      w_mgs = resultsBinaryMatchThenGreedySearch(bmfged)$indicTs
      mg_imbalances = apply(w_mgs, 1, function(w){compute_objective_val(X, w, objective = objective, inv_cov_X = inv_cov_X)})

      #D = MR
      w_mrs = resultsBinaryMatchSearch(bmeds, num_vectors = nR * 1 / rerand_threshold, objective = objective)$indicTs
      bmeds_rerand_balances = apply(w_mrs, 1, function(w){compute_objective_val(X, w, objective = objective, inv_cov_X = inv_cov_X)})
      mr_imbalances = bmeds_rerand_balances[order(bmeds_rerand_balances)[1 : nR]]
      
      #save
      res = rbind(res,
        data.frame(imbalance = bcrd_imbalances, method = "BCRD", n = n, p = p, nx = nx),
        data.frame(imbalance = r_imbalances,    method = "R",    n = n, p = p, nx = nx),
        data.frame(imbalance = m_imbalances,    method = "M",    n = n, p = p, nx = nx),
        data.frame(imbalance = g_imbalances,    method = "G",    n = n, p = p, nx = nx),
        data.frame(imbalance = mg_imbalances,   method = "MG",   n = n, p = p, nx = nx),
        data.frame(imbalance = mr_imbalances,   method = "MR",   n = n, p = p, nx = nx)
      )
      save(res, file = "imbalance_experiments.RData")
      
      res$method = factor(res$method, design_levels_properly_ordered)
      res_summary = res %>% group_by(n, p, method) %>% summarize(avg_log10_imbalance = mean(log10(imbalance)))
      print(ggplot(res_summary) +
              geom_line(aes(x = log10(n), y = avg_log10_imbalance, col = method)) +
              facet_wrap(p ~ ., scales = "free_x", ncol = 4) +
              scale_color_manual(values = alpha(manual_colors, 1)) +
              xlab("log10 n") +
              ylab("log10 d_M covariate imbalance"))
    }
  }
}

res$method = factor(res$method, design_levels_properly_ordered)
avg_summary = res %>% group_by(method, n, p) %>% summarize(med = median(log10(imbalance)))
print(ggplot(res) +
        geom_density(aes(x = log10(imbalance), col = method, fill = method), position = "identity", alpha = 0.2) +
        facet_wrap(n ~ p, 
                   scales = "free",
                   ncol = length(ps),
                   labeller = function(labs){label_both(labs, multi_line = FALSE)}) +
        scale_y_continuous(breaks = NULL) +
        scale_color_manual(values = alpha(manual_colors, 0)) +
        scale_fill_manual(values = manual_colors) +
        xlab("log10 Mahal distance imbalance in x") +
        geom_vline(data = avg_summary, aes(xintercept = med, col = method), lwd = 2, alpha = 0.5, linetype = "dashed"))


res_summary_nx = res %>% group_by(n, p, method, nx) %>% summarize(avg_log10_imbalance = mean(log10(imbalance)))
ggplot_obj = ggplot(res_summary_nx) +
  aes(x = log10(n), y = avg_log10_imbalance, col = method) +
  xlab("log10 n") +
  ylab("log10 d_M covariate imbalance") +
  scale_x_continuous(breaks = c(2, 3)) +
  scale_color_manual(values = alpha(manual_colors, 1)) +
  facet_wrap(p ~ ., scales = "free_x", ncol = 4 ,labeller = function(labs){label_both(labs, multi_line = FALSE)})  

for (nx_i in 1 : nX){
  ggplot_obj = ggplot_obj +
    geom_line(data = res_summary_nx %>% filter(nx == nx_i)) 
}
ggplot_obj


# res_summary = res %>% group_by(n, p, method) %>% summarize(avg_log10_imbalance = mean(log10(imbalance)))
# print(ggplot(res_summary) +
#         geom_line(aes(x = log10(n), y = avg_log10_imbalance, col = method)) +
#         facet_wrap(p ~ ., scales = "free_x", ncol = 4) +
#         scale_color_manual(values = alpha(manual_colors, 1)) +
#         xlab("log10 n") +
#         ylab("log10 d_M covariate imbalance"))

# res %>% group_by(method, n, p) %>% summarize(avg = mean(log10(imbalance)))
# 
# res = data.table(res)
# res = res[imbalance != 0]
# res[method == "pair-matching-then-greedy" & n == 500, imbalance]
# sort(res[method == "pair-matching-then-greedy" & n == 500, imbalance])
# 
# res[, .(avg_log10_imbalance = mean(log10(imbalance))), by = c("n", "p", "method")]
# 
# t.test(res[n == 1000 & method == "greedy", log10(imbalance)],
#        res[n == 1000 & method == "pair-matching-then-greedy", log10(imbalance)])

res_summary_nx %<>% filter(!is.infinite(avg_log10_imbalance))
power_res = data.frame(p = factor(), method = factor(), slope = numeric())
power_res_wide = matrix(NA, nrow = length(levels(res$method)), ncol = length(ps))
rownames(power_res_wide) = levels(res$method)
colnames(power_res_wide) = ps
for (i_p in 1 : length(ps)){
  coefs = coef(lm(avg_log10_imbalance ~ log10(n) * method, (res_summary_nx[res_summary_nx$p == ps[i_p], ])))[c(2,8:12)]
  coefs = coefs + coefs[1]
  power_res = rbind(power_res,
    data.frame(
      p = p, 
      method = levels(res$method), 
      slope = coefs
    )               
  )
  power_res_wide[, i_p] = coefs
}
power_res_wide = rbind(power_res_wide, power_res_wide[3,] + power_res_wide[4,])
xtable(-power_res_wide, digits = 3)


power_res$method = factor(power_res$method, 
                    design_levels_properly_ordered)

ggplot(power_res) + 
  geom_bar(aes(x = method, y = -slope, fill = method), position = "dodge2", stat = 'identity')  +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(breaks = 1:9) +
  scale_fill_manual(values = alpha(manual_colors, 1)) +
  xlab("") + ylab("r : d_M = O(n^-r)") +
  facet_wrap(p ~ ., ncol = 4, labeller = function(labs){label_both(labs, multi_line = FALSE)})

# summary(lm(avg_log10_imbalance ~ log10(n) * as.factor(p) * method, res_summary))

res_summary_nx %>% 
  filter(n == 100) %>% 
  group_by(p, method) %>% 
  summarize(avg_log10_imbalance = mean(avg_log10_imbalance)) %>%
  filter(method %in% c("G", "MG")) %>%
  data.frame
