



res$design = relevel(res$design, "matching_then_greedy")

res_sub = res
levels(res_sub$design) = design_codes
levels(res_sub$model)  = model_codes
res_sub$design = factor(res_sub$design, levels = design_levels_properly_ordered)



######
manual_colors = c("black", "purple", "red", "green", "blue", "orange")
design_codes = c("MG", "BCRD", "R", "M", "MR", "G")

model_codes = c("Z", "L", "LsNL", "LNL", "NL", "LH", "LNLH", "NLH")


avg_summary_nx = avg_summary_nx_neps %>% group_by(model, n, p, design, nx) %>% summarize(squared_error_avg_nx = mean(squared_error_avg_nx_neps))
avg_summary = avg_summary_nx %>% group_by(model, n, p, design) %>% summarize(squared_error_avg = mean(squared_error_avg_nx))


avg_summary_100 = data.table(avg_summary)[n == 100 & p == 2]
levels(avg_summary_100$design) = design_codes
levels(avg_summary_100$model)  = model_codes

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

ggplot(avg_summary_nx %>% filter(p == 10)) +
  aes(x = squared_error_avg_nx, col = design, fill = design) +
  geom_density(aes(col = design), alpha = 0.1) +
  facet_grid(rows = "model", 
             scales = "free_x", #labeller = function(labs){label_both(labs, multi_line = FALSE)}
  ) +
  scale_x_log10(limits = c( #breaks = NULL, 
    quantile(avg_summary_nx$squared_error_avg_nx, 0.01), #0.1
    quantile(avg_summary_nx$squared_error_avg_nx, 0.99) #0.99
  )) + 
  ylim(0, 10) + 
  # scale_y_continuous(breaks = NULL) +
  # scale_color_lancet() + scale_fill_lancet() +
  xlab("squared error") +
  scale_color_manual(values = alpha(manual_colors, 0)) +
  scale_fill_manual(values = manual_colors) +
  geom_vline(data = avg_summary %>% filter(p == 10),
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



























