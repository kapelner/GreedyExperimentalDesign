pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table, viridis, RColorBrewer, ggsci, lmtest, sandwich, scales, xtable)
pacman::p_load_gh("zeehio/facetscales")

###run the simulation first
options(dplyr.tibble_min = 100)

source("illustration_params.R")

levels(res$model) = c("Z", "L", "LsNL", "LNL", "NL", "HL",  "HLNL", "HNL")
res$model = factor(res$model, model_codes_properly_ordered)

levels(res$design) = c("BCRD", "R", "M", "MG", "MR", "G")
res$design = factor(res$design, design_levels_properly_ordered)

res %<>% 
  mutate(log10_squared_error = log10(squared_error))
res_neyman_frame = res %>% filter(nx == 1 & neps == 1) 

summary_neyman_frame = res_neyman_frame %>% 
  group_by(model, p, design) %>% 
  summarize(
    avg_squared_error = mean(squared_error), 
    sd_squared_error = sd(squared_error), 
    se_squared_error = sd(squared_error) / sqrt(n()),
    avg_log10_squared_error = mean(log10_squared_error), 
    sd_log10_squared_error = sd(log10_squared_error), 
    se_log10_squared_error = sd(log10_squared_error) / sqrt(n())
  )


####Fig 1

# res_sub[, design_col := manual_colors[as.numeric(design)]]
#manual jitter
# jit = 0.00030
# avg_summary_100[1 : 6, "squared_error_avg"] = 0.01 + jit * (0:5)
# # 
# # manual jitter (uniform sim)
# avg_summary_100[11 : 12, "squared_error_avg"] = 0.01 + jit * (0:2)
# 
# avg_summary_100[16 : 17, "squared_error_avg"] = 0.0126 + jit * (0:1)
# avg_summary_100[c(19, 23), "squared_error_avg"] = 0.0176 + jit * 1.5 * (0:1)
# avg_summary_100[c(25, 28, 29), "squared_error_avg"] = 0.0170 + jit * 1.5 * (0:2)

#manual jitter (normal sim)
# avg_summary_100[10 : 12, "squared_error_avg"] = 0.01 + jit * (0:2)
# avg_summary_100[16 : 17, "squared_error_avg"] = 0.01 + jit * 1 * (0:1)
# avg_summary_100[22 : 23, "squared_error_avg"] = 0.054 + jit * 5 * (0:1)
# avg_summary_100[25 : 26, "squared_error_avg"] = 0.21 + jit * 20 * (0:1)
# avg_summary_100[28 : 29, "squared_error_avg"] = 0.054 + jit * 5 * (0:1)

x_custom_scales = list(
  Z = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(res_neyman_frame %>% filter(model == "Z") %>% select(log10_squared_error) %>% as.matrix, 
        0.15), #0.1
      quantile(res_neyman_frame %>% filter(model == "Z") %>% select(log10_squared_error) %>% as.matrix, 
        0.95) #0.1
      #0.99
    )
  ),
  L = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(res_neyman_frame %>% filter(model == "L") %>% select(log10_squared_error) %>% as.matrix, 
        0.15), #0.1
      quantile(res_neyman_frame %>% filter(model == "L") %>% select(log10_squared_error) %>% as.matrix, 
        0.99) #0.1
      #0.99
    )
  ),
  LsNL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(res_neyman_frame %>% filter(model == "LsNL") %>% select(log10_squared_error) %>% as.matrix, 
        0.15), #0.1
      quantile(res_neyman_frame %>% filter(model == "LsNL") %>% select(log10_squared_error) %>% as.matrix, 
        0.99) #0.1
      #0.99
    )
  ),
  LNL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(res_neyman_frame %>% filter(model == "LNL") %>% select(log10_squared_error) %>% as.matrix, 
        0.15), #0.1
      quantile(res_neyman_frame %>% filter(model == "LNL") %>% select(log10_squared_error) %>% as.matrix, 
        0.97) #0.1
      #0.99
    )
  ),
  NL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(res_neyman_frame %>% filter(model == "NL") %>% select(log10_squared_error) %>% as.matrix, 
        0.15), #0.1
      quantile(res_neyman_frame %>% filter(model == "NL") %>% select(log10_squared_error) %>% as.matrix, 
        0.93) #0.1
      #0.99
    )
  ),
  HL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(res_neyman_frame %>% filter(model == "HL") %>% select(log10_squared_error) %>% as.matrix, 
        0.23), #0.1
      quantile(res_neyman_frame %>% filter(model == "HL") %>% select(log10_squared_error) %>% as.matrix, 
        0.97) #0.1
      #0.99
    )
  ),
  HLNL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(res_neyman_frame %>% filter(model == "HLNL") %>% select(log10_squared_error) %>% as.matrix, 
        0.23), #0.1
      quantile(res_neyman_frame %>% filter(model == "HLNL") %>% select(log10_squared_error) %>% as.matrix, 
        0.97) #0.1
      #0.99
    )
  ),
  HNL = scale_x_continuous(
    breaks = -4 : 4,
    limits = c( #breaks = NULL,
      quantile(res_neyman_frame %>% filter(model == "HNL") %>% select(log10_squared_error) %>% as.matrix, 
        0.10), #0.1
      quantile(res_neyman_frame %>% filter(model == "HNL") %>% select(log10_squared_error) %>% as.matrix, 
        0.97) #0.1
      #0.99
    )
  )
)

ggplot(res_neyman_frame) +
  aes(x = log10_squared_error, col = design, fill = design) +
  geom_density(aes(col = design), alpha = 0.1, trim = FALSE) +
  facet_grid_sc(
    rows = p ~ model, 
    scales = list(x = x_custom_scales, y = "free_y"), 
    labeller = labeller(p = label_both, model = label_both)
  ) +
  scale_y_continuous(breaks = NULL) +
  xlab("log10 squared error") +
  scale_color_manual(values = alpha(manual_colors, 0)) +
  scale_fill_manual(values = manual_colors) +
  geom_vline(data = summary_neyman_frame,
             aes(xintercept = avg_log10_squared_error, color = design),
             lwd = 0.75, alpha = 0.7, linetype = "solid") + #2 * summary_neyman_frame$se_log10_squared_error
  guides(color = FALSE, design = TRUE, fill = guide_legend(override.aes = list(alpha = 1)))






############### order performance by avg in bucket

summary_neyman_frame %<>%
  arrange(model, p, avg_squared_error)

summary_neyman_frame_comps2 = summary_neyman_frame %>%
  group_by(p, model) %>%
  mutate(placement = 1 : n()) %>%
  select(p, model, design, placement, avg_squared_error) %>%
  mutate(vs_1 = avg_squared_error / nth(avg_squared_error, n = 1)) %>%
  mutate(vs_2 = avg_squared_error / nth(avg_squared_error, n = 2)) %>%
  mutate(vs_3 = avg_squared_error / nth(avg_squared_error, n = 3)) %>%
  mutate(vs_4 = avg_squared_error / nth(avg_squared_error, n = 4)) %>%
  mutate(vs_5 = avg_squared_error / nth(avg_squared_error, n = 5)) %>%
  mutate(vs_1 = ifelse(vs_1 <= 0, NA, vs_1), vs_2 = ifelse(vs_2 <= 0, NA, vs_2), vs_3 = ifelse(vs_3 <= 0, NA, vs_3), vs_4 = ifelse(vs_4 <= 0, NA, vs_4), vs_5 = ifelse(vs_5 <= 0, NA, vs_5))
summary_neyman_frame_comps2 %>% filter(design == "MG") %>% arrange(placement)


summary_neyman_frame_comps1 = summary_neyman_frame %>%
  group_by(p, model) %>%
  mutate(placement = 1 : n()) %>%
  select(p, model, design, placement, avg_squared_error) %>%
  mutate(vs_6 = 100 - avg_squared_error / nth(avg_squared_error, n = 6) * 100) %>%
  mutate(vs_5 = 100 - avg_squared_error / nth(avg_squared_error, n = 5) * 100) %>%
  mutate(vs_4 = 100 - avg_squared_error / nth(avg_squared_error, n = 4) * 100) %>%
  mutate(vs_3 = 100 - avg_squared_error / nth(avg_squared_error, n = 3) * 100) %>%
  mutate(vs_2 = 100 - avg_squared_error / nth(avg_squared_error, n = 2) * 100) %>%
  mutate(vs_6 = ifelse(vs_6 <= 0, NA, vs_6), vs_2 = ifelse(vs_2 <= 0, NA, vs_2), vs_3 = ifelse(vs_3 <= 0, NA, vs_3), vs_4 = ifelse(vs_4 <= 0, NA, vs_4), vs_5 = ifelse(vs_5 <= 0, NA, vs_5))
summary_neyman_frame_comps1
summary_neyman_frame_comps1 %>% filter(design == "MG") %>% arrange(placement)


summary_neyman_frame_comps3 = summary_neyman_frame %>%
  group_by(p, model) %>%
  mutate(placement = 1 : n()) %>%
  select(p, model, design, placement, avg_squared_error) %>%
  mutate(mult_ideal = avg_squared_error / ideal_sq_err) %>%
  ungroup() %>%
  mutate(pct_better_than_worst = summary_neyman_frame_comps1$vs_6) %>%
  mutate(pct_worse_than_best = summary_neyman_frame_comps2$vs_1)


#########Tables 2,3, S1

p = 10
for (model in model_codes_properly_ordered){
  summary_neyman_frame_p_mod = summary_neyman_frame_comps3[summary_neyman_frame_comps3$model == model & summary_neyman_frame_comps3$p == p, ]
  
  cat("p", p, "model", model, "\n")
  print(xtable(summary_neyman_frame_p_mod %>% select(design, avg_squared_error, mult_ideal, pct_better_than_worst, pct_worse_than_best), digits = 3), include.rownames = FALSE)
} 


