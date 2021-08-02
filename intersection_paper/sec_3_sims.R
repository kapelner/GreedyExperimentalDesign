####population model sims
nX = 5
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

kernel_objective_functions = c(kernel_names, "mahalanobis_plus_gaussian", "mahalanobis_times_gaussian")

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
  group_by(design, kernel, n, p, model, covariate_distribution) %>%
  summarize(ase = mean(squared_error))
write_csv(res_summary, file = "basic_results_covariate_distribution.csv")

res_summary = res %>%
  group_by(design, kernel, n, p, model) %>%
  summarize(ase = mean(squared_error))
write_csv(res_summary, file = "basic_results.csv")