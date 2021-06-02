ns = 100 #c(32, 100, 132, 200)
ps = c(2, 5, 10)
nX = 50
nEPS = 50
nR = 500
nR0 = 5000
rerand_threshold = 0.01
nC = 10
betaT = 1
sigma_e = 0.5

all_betas_and_correlations = list()
all_betas_and_correlations$zero_effects = list(
  beta_x1 = 0,
  beta_x2 = 0,
  beta_x1sq = 0,
  beta_x2sq = 0,
  beta_x1_x2 = 0,
  beta_x1_T = 0,
  beta_x2_T = 0
)
all_betas_and_correlations$linear_uniform = list(
  beta_x1 = 3,
  beta_x2 = 3,
  beta_x1sq = 0,
  beta_x2sq = 0,
  beta_x1_x2 = 0,
  beta_x1_T = 0,
  beta_x2_T = 0
)
all_betas_and_correlations$one_quadratic = list(
  beta_x1 = 3,
  beta_x2 = 3,
  beta_x1sq = 1,
  beta_x2sq = 0,
  beta_x1_x2 = 0,
  beta_x1_T = 0,
  beta_x2_T = 0
)
all_betas_and_correlations$quadratics_and_interaction = list(
  beta_x1 = 3,
  beta_x2 = 3,
  beta_x1sq = 1,
  beta_x2sq = 1,
  beta_x1_x2 = 1,
  beta_x1_T = 0,
  beta_x2_T = 0
)
all_betas_and_correlations$purely_nonlinear = list(
  beta_x1 = 0,
  beta_x2 = 0,
  beta_x1sq = 1,
  beta_x2sq = 1,
  beta_x1_x2 = 1,
  beta_x1_T = 0,
  beta_x2_T = 0
)
all_betas_and_correlations$linear_uniform_plus_hetero_tx_effect = list(
  beta_x1 = 3,
  beta_x2 = 3,
  beta_x1sq = 0,
  beta_x2sq = 0,
  beta_x1_x2 = 0,
  beta_x1_T = 1,
  beta_x2_T = 1
)
all_betas_and_correlations$quadratics_and_interaction_plus_hetero_tx_effect = list(
  beta_x1 = 3,
  beta_x2 = 3,
  beta_x1sq = 1,
  beta_x2sq = 1,
  beta_x1_x2 = 1,
  beta_x1_T = 1,
  beta_x2_T = 1
)
all_betas_and_correlations$purely_nonlinear_plus_hetero_tx_effect = list(
  beta_x1 = 0,
  beta_x2 = 0,
  beta_x1sq = 1,
  beta_x2sq = 1,
  beta_x1_x2 = 1,
  beta_x1_T = 1,
  beta_x2_T = 1
)


objective = "abs_sum_diff"