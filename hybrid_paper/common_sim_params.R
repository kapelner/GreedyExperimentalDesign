#https://console.aws.amazon.com/cost-management/home?#/custom?reportType=CostUsage&chartStyle=Stack&groupBy=InstanceType&forecastTimeRangeOption=None&hasBlended=false&hasAmortized=false&excludeDiscounts=true&usageAs=usageQuantity&excludeCategorizedResources=false&excludeTaggedResources=false&excludeForecast=false&reportName=&timeRangeOption=Last6Months&startDate=2020-12-01&endDate=2021-05-31&granularity=Monthly&filter=%5B%5D&isTemplate=true


n = 100
ps = c(2, 5, 10)
rerand_threshold = 0.01
betaT = 1
sigma_e = 0.5

nC = 6

set.seed(1)

Xalls = list()
for (nx in 1 : nX){
  Xalls[[nx]] = matrix(runif(n * max(ps), -sqrt(3), sqrt(3)), nrow = n)
  # Xtemp = matrix(rnorm(n * p, 0, 1), nrow = n)
  # Xtemp = data.matrix(MASS::Pima.tr[1 : n, 1 : p])
  # Xalls[[nx]] = apply(Xtemp, 2, function(xj){(xj - mean(xj)) / sd(xj)}) #stdize
  # rm(Xtemp)
}

epsilon_alls = list()
for (neps in 1 : nEPS){
  epsilon_alls[[neps]] = rnorm(n, 0, sigma_e)
}

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