pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table, viridis, RColorBrewer, ggsci, lmtest, sandwich, scales, xtable)
pacman::p_load_gh("zeehio/facetscales")

###run the simulation first
source("illustration_params.R")

avg_summary_nx_neps = res %>% 
  group_by(model, n, p, design, nx, neps) %>% 
  summarize(squared_error_avg_nx_neps = mean(squared_error))

###########################
avg_summary_nx_neps %<>% ungroup() %>% select(-n)
avg_summary_nx_neps %<>%
  mutate(model = factor(model), design = factor(design), p = factor(p), nx = factor(nx), neps = factor(neps))


levels(avg_summary_nx_neps$model) = c("L", "HL", "LsNL", "NL", "HNL", "LNL", "HLNL", "Z")
avg_summary_nx_neps$model = factor(avg_summary_nx_neps$model, model_codes_properly_ordered)
levels(avg_summary_nx_neps$design) = c("BCRD", "G", "M", "MG", "MR", "R")
avg_summary_nx_neps$design = factor(avg_summary_nx_neps$design, design_levels_properly_ordered)

avg_summary_nx_neps %<>%
  mutate(log10_squared_error_avg_nx_neps = log10(squared_error_avg_nx_neps))

summary_pop_frame = avg_summary_nx_neps %>%
  group_by(model, design, p) %>%
  summarize(squared_error_avg = mean(squared_error_avg_nx_neps), squared_error_log10_avg = mean(log10_squared_error_avg_nx_neps))

#Fig S1

x_custom_scales = list(
  Z = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(avg_summary_nx_neps %>% filter(model == "Z") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.15), #0.1
      quantile(avg_summary_nx_neps %>% filter(model == "Z") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.95) #0.1
      #0.99
    )
  ),
  L = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(avg_summary_nx_neps %>% filter(model == "L") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.05), #0.1
      quantile(avg_summary_nx_neps %>% filter(model == "L") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.99) #0.1
      #0.99
    )
  ),
  LsNL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(avg_summary_nx_neps %>% filter(model == "LsNL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.10), #0.1
      quantile(avg_summary_nx_neps %>% filter(model == "LsNL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.99) #0.1
      #0.99
    )
  ),
  LNL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(avg_summary_nx_neps %>% filter(model == "LNL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.10), #0.1
      quantile(avg_summary_nx_neps %>% filter(model == "LNL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.97) #0.1
      #0.99
    )
  ),
  NL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(avg_summary_nx_neps %>% filter(model == "NL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.01), #0.1
      quantile(avg_summary_nx_neps %>% filter(model == "NL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.999) #0.1
      #0.99
    )
  ),
  HL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(avg_summary_nx_neps %>% filter(model == "HL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.1), #0.1
      quantile(avg_summary_nx_neps %>% filter(model == "HL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.98) #0.1
      #0.99
    )
  ),
  HLNL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(avg_summary_nx_neps %>% filter(model == "HLNL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.10), #0.1
      quantile(avg_summary_nx_neps %>% filter(model == "HLNL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.97) #0.1
      #0.99
    )
  ),
  HNL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(avg_summary_nx_neps %>% filter(model == "HNL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.05), #0.1
      quantile(avg_summary_nx_neps %>% filter(model == "HNL") %>% select(log10_squared_error_avg_nx_neps) %>% as.matrix, 
               0.97) #0.1
      #0.99
    )
  )
)


ggplot(avg_summary_nx_neps %>% sample_n(200000)) +
  aes(x = log10_squared_error_avg_nx_neps, col = design, fill = design) +
  geom_density(aes(col = design), alpha = 0.1, trim = FALSE) +
  facet_grid_sc(
    rows = p ~ model, 
    scales = list(x = x_custom_scales, y = "free_y"), 
    labeller = labeller(p = label_both, model = label_both)
  ) +
  scale_y_continuous(breaks = NULL) +
  # scale_color_lancet() + scale_fill_lancet() +
  xlab("log10 squared error") +
  scale_color_manual(values = alpha(manual_colors, 0)) +
  scale_fill_manual(values = manual_colors) +
  geom_vline(data = summary_pop_frame,
             aes(xintercept = squared_error_log10_avg, color = design),
             lwd = 0.75, alpha = 0.7, linetype = "solid") + #2 * summary_neyman_frame$se_log10_squared_error
  guides(color = FALSE, design = TRUE, fill = guide_legend(override.aes = list(alpha = 1)))



############### order performance by avg in bucket

summary_pop_frame %<>%
  arrange(model, p, squared_error_avg)

summary_pop_frame_comps2 = summary_pop_frame %>%
  group_by(p, model) %>%
  mutate(placement = 1 : n()) %>%
  select(p, model, design, placement, squared_error_avg) %>%
  mutate(vs_1 = squared_error_avg / nth(squared_error_avg, n = 1)) %>%
  mutate(vs_2 = squared_error_avg / nth(squared_error_avg, n = 2)) %>%
  mutate(vs_3 = squared_error_avg / nth(squared_error_avg, n = 3)) %>%
  mutate(vs_4 = squared_error_avg / nth(squared_error_avg, n = 4)) %>%
  mutate(vs_5 = squared_error_avg / nth(squared_error_avg, n = 5)) %>%
  mutate(vs_1 = ifelse(vs_1 <= 0, NA, vs_1), vs_2 = ifelse(vs_2 <= 0, NA, vs_2), vs_3 = ifelse(vs_3 <= 0, NA, vs_3), vs_4 = ifelse(vs_4 <= 0, NA, vs_4), vs_5 = ifelse(vs_5 <= 0, NA, vs_5))

summary_pop_frame_comps1 = summary_pop_frame %>%
  group_by(p, model) %>%
  mutate(placement = 1 : n()) %>%
  select(p, model, design, placement, squared_error_avg) %>%
  mutate(vs_6 = 100 - squared_error_avg / nth(squared_error_avg, n = 6) * 100) %>%
  mutate(vs_5 = 100 - squared_error_avg / nth(squared_error_avg, n = 5) * 100) %>%
  mutate(vs_4 = 100 - squared_error_avg / nth(squared_error_avg, n = 4) * 100) %>%
  mutate(vs_3 = 100 - squared_error_avg / nth(squared_error_avg, n = 3) * 100) %>%
  mutate(vs_2 = 100 - squared_error_avg / nth(squared_error_avg, n = 2) * 100) %>%
  mutate(vs_6 = ifelse(vs_6 <= 0, NA, vs_6), vs_2 = ifelse(vs_2 <= 0, NA, vs_2), vs_3 = ifelse(vs_3 <= 0, NA, vs_3), vs_4 = ifelse(vs_4 <= 0, NA, vs_4), vs_5 = ifelse(vs_5 <= 0, NA, vs_5))

n = 100
sigma_e = 0.5
ideal_sq_err = 2 * sigma_e^2 / (n / 2)

summary_pop_frame_comps3 = summary_pop_frame %>%
  group_by(p, model) %>%
  mutate(placement = 1 : n()) %>%
  select(p, model, design, placement, squared_error_avg) %>%
  mutate(mult_ideal = squared_error_avg / ideal_sq_err) %>%
  ungroup() %>%
  mutate(pct_better_than_worst = summary_pop_frame_comps1$vs_6) %>%
  mutate(pct_worse_than_best = summary_pop_frame_comps2$vs_1)

#########Tables S2, S3, S4

p = 10
for (model in model_codes_properly_ordered){
  summary_pop_frame_p_mod = summary_pop_frame_comps3[summary_pop_frame_comps3$model == model & summary_pop_frame_comps3$p == p, ]
  
  cat("p", p, "model", model, "\n")
  print(xtable(summary_pop_frame_p_mod %>% select(design, squared_error_avg, mult_ideal, pct_better_than_worst, pct_worse_than_best), digits = 3), include.rownames = FALSE)
}





########### check if performance of the designs in model Z is statistically different
TK_results = list()
for (model_name in model_codes_properly_ordered){
  res_mod = aov(squared_error_avg_nx_neps ~ 0 + design, data = avg_summary_nx_neps[avg_summary_nx_neps$model == model_name, ])
  coef(res_mod)
  TK_results[[model_name]] = TukeyHSD(x = res_mod, which = "design", conf.level = 0.95)$design
}

all_pval_res = matrix("", nrow = length(design_levels_properly_ordered), ncol = 0)
rownames(all_pval_res) = design_levels_properly_ordered

for (model in model_codes_properly_ordered){
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
######conclusion: Z has no statistically significant differences
































