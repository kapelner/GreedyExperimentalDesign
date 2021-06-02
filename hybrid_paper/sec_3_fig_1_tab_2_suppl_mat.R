options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table, viridis, RColorBrewer, ggsci, lmtest, sandwich)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign

set.seed(1)
ns = c(32, 100, 132, 200)
ps = 2
nX = 50
nR = 500
rerand_threshold = 0.01
nC = 10
objective = "abs_sum_diff"
res = data.table(n = factor(), p = factor(), nx = factor(), nsim = factor(), model = factor(), design = factor(), squared_error = numeric())

all_betas_and_correlations = list()
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 
  betas = c(0, 0, 0, 0, 0)))) #ZERO EFFECTS
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 
  betas = c(3, 3, 0, 0, 0)))) #LINEAR EVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 
  betas = c(3, 3, 1, 0, 0)))) #QUAD MORE UNEVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 
  betas = c(3, 3, 1, 1, 1)))) #QUAD MORE UNEVEN
all_betas_and_correlations = c(all_betas_and_correlations, list(setting = list(rho = 0, 
  betas = c(0, 0, 1, 1, 1)))) #COMPLETELY NONLINEAR
betaT = 1
sigma_e = 0.5

for (n in ns){
  for (p in ps){
    for (nx in 1 : nX){
      cat("n", n, "p", p, "nx", nx, "\n")
      X = matrix(runif(n * p, -sqrt(3), sqrt(3)), nrow = n)
      # X = matrix(rnorm(n * p, 0, 1), nrow = n)
      # X = data.matrix(MASS::Pima.tr[1 : n, 1 : p])
      #stdize
      # X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
      epsilons = rnorm(n, 0, sigma_e)
      
      ##ADD BCRD HERE
      bcrds = complete_randomization_with_forced_balanced(n, nR)
      
      bmeds = binaryMatchExperimentalDesignSearch(X)
      bmeds_res = resultsBinaryMatchSearch(bmeds, num_vectors = min(nR, 2^(n/2)), objective = objective)
      # dim(bmeds_res$indicTs)
      
      bmfged = binaryMatchFollowedByGreedyExperimentalDesignSearch(X, objective = objective, max_designs = nR, wait = TRUE, num_cores = nC)
      # bmfged$binary_match_design
      # bmfged$greedy_design
      bmfged_res = resultsBinaryMatchThenGreedySearch(bmfged)
      # dim(bmfged_res$indicTs)
      
      #rerand
      bcrds_rerand = complete_randomization_with_forced_balanced(n, nR * 1 / rerand_threshold)
      bcrds_rerand_balances = apply(bcrds_rerand, 1, function(w){compute_objective_val(X, w)})
      rerands = bcrds_rerand[order(bcrds_rerand_balances)[1 : nR], ]
      
      #matching then rerand
      bmeds_rerand = resultsBinaryMatchSearch(bmeds, num_vectors = nR * 1 / rerand_threshold, objective = objective)$indicTs
      bmeds_rerand_balances = apply(bmeds_rerand, 1, function(w){compute_objective_val(X, w)})
      matching_then_rerands = bmeds_rerand[order(bmeds_rerand_balances)[1 : nR], ]
      
      ged = initGreedyExperimentalDesignObject(X, objective = objective, wait = TRUE, num_cores = nC)
      ged_res = resultsGreedySearch(ged, max_vectors = nR)
      # dim(ged_res$ending_indicTs)
      
      
      for (nsim in 1 : nR){
        for (model_setting in all_betas_and_correlations){
          betas = model_setting$betas
          
          y = betas[1] * X[, 1] + 
            betas[2] * X[, 2] + 
            betas[3] * X[, 1]^2 + 
            betas[4] * X[, 2]^2 + 
            betas[5] * X[, 1] * X[, 2] + 
            epsilons
          
          #bcrd's
          w = bcrds[nsim, ]
          y_bcrd = y + betaT * w
          betaT_hat_cbrd = mean(y_bcrd[w == 1]) - mean(y_bcrd[w == 0])
          #rerand
          w = rerands[nsim, ]
          y_rerand = y + betaT * w
          betaT_hat_rerand = mean(y_rerand[w == 1]) - mean(y_rerand[w == 0])        
          #matching
          w = bmeds_res$indicTs[nsim,]
          y_matching = y + betaT * w
          betaT_hat_matching = mean(y_matching[w == 1]) - mean(y_matching[w == 0])
          #matching then greedy
          w = bmfged_res$indicTs[nsim,]
          y_matching_then_greedy = y + betaT * w
          betaT_hat_matching_then_greedy = mean(y_matching_then_greedy[w == 1]) - mean(y_matching_then_greedy[w == 0])
          #matching then rerand
          w = matching_then_rerands[nsim,]
          y_matching_then_rerand = y + betaT * w
          betaT_hat_matching_then_rerand = mean(y_matching_then_rerand[w == 1]) - mean(y_matching_then_rerand[w == 0])
          
          #just greedy
          w = ged_res$ending_indicTs[nsim,]
          y_greedy = y + betaT * w
          betaT_greedy = mean(y_greedy[w == 1]) - mean(y_greedy[w == 0])
  
          res = rbind(res, data.frame(design = "bcrd",                 squared_error = (betaT_hat_cbrd - betaT)^2,                 n = n, p = p, nx = nx, nsim = nsim, model = paste0(betas, collapse = "")))
          res = rbind(res, data.frame(design = "rerand",               squared_error = (betaT_hat_rerand - betaT)^2,               n = n, p = p, nx = nx, nsim = nsim, model = paste0(betas, collapse = "")))
          res = rbind(res, data.frame(design = "matching",             squared_error = (betaT_hat_matching - betaT)^2,             n = n, p = p, nx = nx, nsim = nsim, model = paste0(betas, collapse = "")))
          res = rbind(res, data.frame(design = "matching_then_greedy", squared_error = (betaT_hat_matching_then_greedy - betaT)^2, n = n, p = p, nx = nx, nsim = nsim, model = paste0(betas, collapse = "")))
          res = rbind(res, data.frame(design = "matching_then_rerand", squared_error = (betaT_hat_matching_then_rerand - betaT)^2, n = n, p = p, nx = nx, nsim = nsim, model = paste0(betas, collapse = "")))
          res = rbind(res, data.frame(design = "greedy",               squared_error = (betaT_greedy - betaT)^2,                   n = n, p = p, nx = nx, nsim = nsim, model = paste0(betas, collapse = "")))
        }     
      }
    }
    
    xL = quantile(res$squared_error, 0.20)
    xU = quantile(res$squared_error, 0.995)
    res_sub = res#[res$design %in% c("bcrd", "matching", "greedy", "matching_then_greedy"), ]
    avg_summary_nx = res %>% group_by(model, n, p, design, nx) %>% summarize(squared_error_avg_nx = mean(squared_error))
    avg_summary = avg_summary_nx %>% group_by(model, n, p, design) %>% summarize(squared_error_avg = mean(squared_error_avg_nx))
    
    print(ggplot(res_sub) +
      geom_density(aes(x = squared_error, col = design, fill = design), position = "identity", alpha = 0.2) +
      facet_wrap(n ~ model, 
                 scales = "free_x", 
                 ncol = length(all_betas_and_correlations),
                 labeller = function(labs){label_both(labs, multi_line = FALSE)}) +
      scale_x_log10(limits=c(min(avg_summary$squared_error_avg), max(avg_summary$squared_error_avg))) + 
      #scale_color_npg() + 
      xlab("squared_error") +
      geom_vline(data = avg_summary, 
        aes(xintercept = squared_error_avg, col = design), 
        lwd = 2, alpha = 0.5, linetype = "dashed"))

    save(res, file = "matching_then_greedy_squared_error_sims_unif.RData")
  }
}







#Fig 1
res$design = relevel(res$design, "matching_then_greedy")
manual_colors = c(
  "black", "purple", "red", "green", "blue", "orange"
)
res_sub = res
levels(res_sub$design) = c("MG", "BCRD", "R", "M", "MR", "G")
model_levels_properly_ordered = c("Z", "L", "LsNL", "LNL", "NL")
levels(res_sub$model)  = model_levels_properly_ordered
design_levels_properly_ordered = c("BCRD", "R", "G", "M", "MR", "MG")
res_sub$design = factor(res_sub$design, levels = design_levels_properly_ordered)

avg_summary_nx = res %>% group_by(model, n, p, design, nx) %>% summarize(squared_error_avg_nx = mean(squared_error))
avg_summary = avg_summary_nx %>% group_by(model, n, p, design) %>% summarize(squared_error_avg = mean(squared_error_avg_nx))

avg_summary_100 = data.table(avg_summary)[n == 100]
levels(avg_summary_100$design) = c("MG", "BCRD", "R", "M", "MR", "G")
levels(avg_summary_100$model)  = c("Z", "L", "LsNL", "LNL", "NL")

# res_sub[, design_col := manual_colors[as.numeric(design)]]
#manual jitter
jit = 0.00030
avg_summary_100[1 : 6, "squared_error_avg"] = 0.01 + jit * (0:5)
# 
# manual jitter (uniform sim)
avg_summary_100[11 : 12, "squared_error_avg"] = 0.01 + jit * (0:2)

avg_summary_100[16 : 17, "squared_error_avg"] = 0.0126 + jit * (0:1)
avg_summary_100[c(19, 23), "squared_error_avg"] = 0.0176 + jit * 1.5 * (0:1)
avg_summary_100[c(25, 28, 29), "squared_error_avg"] = 0.0170 + jit * 1.5 * (0:2)

#manual jitter (normal sim)
# avg_summary_100[10 : 12, "squared_error_avg"] = 0.01 + jit * (0:2)
# avg_summary_100[16 : 17, "squared_error_avg"] = 0.01 + jit * 1 * (0:1)
# avg_summary_100[22 : 23, "squared_error_avg"] = 0.054 + jit * 5 * (0:1)
# avg_summary_100[25 : 26, "squared_error_avg"] = 0.21 + jit * 20 * (0:1)
# avg_summary_100[28 : 29, "squared_error_avg"] = 0.054 + jit * 5 * (0:1)

ggplot(res_sub[n == 100]) +
  aes(x = squared_error, col = design, fill = design) +
  geom_density(aes(col = design), alpha = 0.1) +
  facet_grid(rows = "model", 
             scales = "free_x", #labeller = function(labs){label_both(labs, multi_line = FALSE)}
             ) +
  scale_x_log10(limits = c( #breaks = NULL, 
    quantile(res_sub[n == 100, squared_error], 0.15), #0.1
    quantile(res_sub[n == 100, squared_error], 0.98) #0.99
  )) + 
  scale_y_continuous(breaks = NULL) +
  # scale_color_lancet() + scale_fill_lancet() +
  xlab("squared error") +
  scale_color_manual(values = alpha(manual_colors, 0)) +
  scale_fill_manual(values = manual_colors) +
  geom_vline(data = avg_summary_100,
             aes(xintercept = squared_error_avg, color = design),
             lwd = 1, alpha = 0.5, linetype = "solid") +
  guides(color = FALSE, design = TRUE, fill = guide_legend(override.aes = list(alpha = 1)))


#table 2
pacman::p_load(xtable)
res_copy = res_sub
res_copy[n == 100, .(squared_error_avg = mean(squared_error)), by = c("model", "design")]
res_copy[n == 100 & model == "Z", .(squared_error_avg = mean(squared_error)), by = design]

TK_results = list()
for (model_name in model_levels_properly_ordered){
  res_mod = aov(squared_error ~ 0 + design, data = res_copy[n == 100 & model == model_name])
  coef(res_mod)
  TK_results[[model_name]] = TukeyHSD(x = res_mod, which = "design", conf.level = 0.95)$design
}

all_pval_res = matrix("", nrow = length(design_levels_properly_ordered), ncol = 0)
rownames(all_pval_res) = design_levels_properly_ordered

for (model in model_levels_properly_ordered){
  TK_res = TK_results[[model]]
  pval_res_model = matrix("", nrow = length(design_levels_properly_ordered), ncol = length(design_levels_properly_ordered) - 1)
  rownames(pval_res_model) = design_levels_properly_ordered
  colnames(pval_res_model) = design_levels_properly_ordered[-1]
  for (design1 in design_levels_properly_ordered){
    if (design1 == "BCRD")
      next
    for (design2 in design_levels_properly_ordered){
      # if (design2 == "MG")
      #   next
      comparison = paste0(design1, "-", design2)
      
      if (comparison %in% rownames(TK_res)){
        #do Bonferroni corrections
        pval = TK_res[comparison, 4]# * 5
        if (pval > 0.05){
          stars = ""
        } else if (pval < 0.001){
          stars = "***"
        } else if (pval < 0.01){
          stars = "**"
        } else { 
          stars = "*"
        }
        if (comparison == "R-M" || comparison == "G-M" ){
          pval_res_model[design1, design2] = stars
        } else {
          pval_res_model[design2, design1] = stars
        }
      }

    }
  }
  cat("\n\n\nmodel:", model, "\n")
  print(xtable(pval_res_model))
  all_pval_res = cbind(all_pval_res, pval_res_model)
}

#table A1-A4
res_mod = aov(squared_error ~ design * model, data = res_copy[n == 32])
summary(res_mod)
xtable(res_mod)
coeftest(res_mod, vcov = vcovHC(res_mod, type="HC1"))

res_mod = lm(squared_error ~ design * model, data = res_copy[n == 100])
summary(res_mod)
xtable(res_mod)
coeftest(res_mod, vcov = vcovHC(res_mod, type="HC1"))

res_mod = lm(squared_error ~ design * model, data = res_copy[n == 132])
summary(res_mod)
xtable(res_mod)
coeftest(res_mod, vcov = vcovHC(res_mod, type="HC1"))

res_mod = lm(squared_error ~ design * model, data = res_copy[n == 200])
summary(res_mod)
xtable(res_mod)
coeftest(res_mod, vcov = vcovHC(res_mod, type="HC1"))
res_mod = lm(squared_error ~ design, data = res[n == 100 & model == "33000"])
summary(res_mod)

# t.test(
#   res[n == 100 & model == "33000" & design == "matching_then_greedy", "squared_error"],
#   res[n == 100 & model == "33000" & design == "matching", "squared_error"]
# )
# t.test(
#   res[n == 100 & model == "33000" & design == "greedy", "squared_error"],
#   res[n == 100 & model == "33000" & design == "rerand", "squared_error"]
# )
# 
