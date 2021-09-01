options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, doParallel, tidyverse, magrittr, data.table)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign

####population model sims
nX = 1
nEPS = 5
nR = 500
n = 100
ps = c(2, 5, 10)
rerand_threshold = 0.01
betaT = 1
sigma_e = 0.5
nC = 5
covariate_distributions = c("uniform", "normal", "exponential")
kernel_names = c("mahalanobis", "poly_2", "exponential", "gaussian", "laplacian", "inv_mult_quad")

kernel_objective_functions = c(
  kernel_names, 
 "mahalanobis_plus_gaussian", 
 "mahalanobis_plus_exponential", 
 "mahalanobis_plus_laplacian", 
 "mahalanobis_plus_inv_mult_quad",
 "mahalanobis_plus_gaussian_plus_inv_mult_quad"
)

set.seed(1)

source("common_sim_params.R")

exp_settings = expand.grid(
  neps = 1 : nEPS,
  nx = 1 : nX,
  p = ps,
  covariate_distribution = covariate_distributions
)

filename = "sec_S3"

source("common_sims.R")



res_summary = res %>%
  group_by(design, kernel, n, p, model, covariate_distribution, estimator) %>%
  summarize(ase = mean(squared_error)) %>%
  ungroup() %>%
  arrange(p, model, covariate_distribution) %>%
  group_by(p, model, covariate_distribution) %>%
  arrange(ase) %>%
  mutate(mult_of_best = ase / first(ase)) %>%
  mutate(prop_of_worst = ase / last(ase)) %>%
  ungroup() %>%
  arrange(p, model, covariate_distribution)
write_csv(res_summary, file = "basic_results_covariate_distribution.csv")

table(res_summary$model)
res_summary %>% filter(model == "linear_uniform") %>% arrange(desc(estimator)) %>% data.frame
res_summary %>% filter(model == "one_quadratic") %>% data.frame
res_summary %>% filter(model == "quadratics_and_interaction") %>% data.frame
res_summary %>% filter(model == "purely_nonlinear") %>% data.frame
res_summary %>% filter(model == "purely_nonlinear_sin") %>% data.frame

res_summary = res %>%
  group_by(design, kernel, n, p, model, estimator) %>%
  summarize(ase = mean(squared_error)) %>%
  ungroup() %>%
  arrange(p, model) %>%
  group_by(p, model) %>%
  arrange(model, estimator, ase) %>%
  mutate(mult_of_best = ase / first(ase)) %>%
  mutate(prop_of_worst = ase / last(ase))
write_csv(res_summary, file = "basic_results.csv")
